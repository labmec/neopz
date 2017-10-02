#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzbfilestream.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "pzpoisson3d.h"
//#include "pzhybridpoisson.h"
#include "pzpoisson3dreferred.h"
#include "mixedpoisson.h"
#include "pzelasmat.h"
#include "pzelasthybrid.h"
#include "pzmat1dlin.h"
#include "TPZVecL2.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZLagrangeMultiplier.h"


#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "TPZCompMeshTools.h"
#include "pzcondensedcompel.h"
#include "pzfunction.h"
#include "pzgraphmesh.h"
#include "pzfmatrix.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"
#include "pzvisualmatrix.h"
#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "pzcheckgeom.h"

#include "TPZMHMeshControl.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZMHMixedHybridMeshControl.h"

#include "meshgen.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>



using namespace std;

struct TRunConfig;

TPZGeoMesh *SimpleMesh(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, TPZVec<long> &coarseindices, int ndiv);

/// Insert material objects for the MHM Mesh solution
void InsertMaterialObjects(TPZMHMeshControl &control);
/// Insert material objects for the MHM-H(div) solution
void InsertMaterialObjects(TPZMHMixedMeshControl &control);

/// Print the elements with geometric information and connect values
void PrintElements(TPZCompMesh *cmesh, std::ostream &out);

/// copy the solution between one computation mesh to the other assuming the geometric elements match
void CopySolution(TPZCompMesh *from, TPZCompMesh *to);

/// unwrap de TPZCondensedCompel and TPZElementGroup elements
void UnwrapMesh(TPZCompMesh *cmesh);

/// function that randomly refines some elements
void RandomRefine(TPZGeoMesh *gmesh, TPZVec<long> &coarseindices, int nref);

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mainskeleton"));
#endif


const int matInterno = 1;
const int matCoarse = 2;
const int skeleton = 4;
const int secondskeleton = 3;
const int matpressure = 6;

const int dirichlet = 0;
const int neumann = 1;

int const bc1=-1;
int const bc2=-2;
int const bc3=-3;
int const bc4=-4;

static void DirichletValidacao(const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &gradres){
    result[0] = loc[0];
}

TAnalyticSolution *example = 0;

int main(int argc, char *argv[])
{
    TExceptionManager except;
    
#ifdef _AUTODIFF
    //    example = new TLaplaceExampleSmooth;
#endif
    
//    example = new TLaplaceExample1;
    
    TRunConfig Configuration;
    
    /// computation type :
    // (0) - compute reference mesh
    // (1) - compute MHM H1 mesh and compute MHM(div) mesh
    int ComputationType = 1;
    /// numhdiv - number of h-refinements
    Configuration.numHDivisions = 1;
    /// PolynomialOrder - p-order
    Configuration.pOrderInternal = 1;
    Configuration.pOrderSkeleton = 1;
    Configuration.numDivSkeleton = 0;
    Configuration.nelxcoarse = 2;
    Configuration.nelycoarse = 2;
    Configuration.Hybridize = 0;
    Configuration.Condensed = 1;

    // to avoid singular internal matrices
    if (Configuration.numHDivisions == 0 && Configuration.pOrderInternal <= Configuration.pOrderSkeleton) {
        Configuration.pOrderInternal = Configuration.pOrderSkeleton+1;
    }


    HDivPiola = 1;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    // tototo
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
//    gRefDBase.InitializeUniformRefPattern(ECube);
    
    TPZGeoMesh *gmesh = 0;
    TPZManVector<long> coarseindices;
    int ndiv = Configuration.numHDivisions;
    TPZManVector<REAL,3> x0(2,0.),x1(2,0.);
    // for using the aligned mesh
    x0.Fill(0.);
    x1.Fill(1.);
    gmesh = SimpleMesh(Configuration.nelxcoarse, Configuration.nelycoarse, x0, x1, coarseindices, ndiv);
    
    TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
    TPZAutoPointer<TPZMHMeshControl> MHM;
    
    std::stringstream MHMPref;
    MHMPref << "MHM";

    /// setting up the datastructure
    {
        TPZMHMeshControl *mhm = new TPZMHMeshControl(gmeshauto);
        mhm->DefinePartitionbyCoarseIndices(coarseindices);
        MHM = mhm;
    }
    MHM->SwitchLagrangeMultiplierSign(true);
    if (Configuration.LagrangeMult) {
        MHM->SetLagrangeAveragePressure(true);
    }
    
    InsertMaterialObjects(MHM);

    MHM->SetInternalPOrder(Configuration.pOrderInternal);
    MHM->SetSkeletonPOrder(Configuration.pOrderSkeleton);
    MHM->DivideSkeletonElements(Configuration.numDivSkeleton);
    if (Configuration.Hybridize)
    {
        MHM->Hybridize();
    }
    
    bool substructure = (bool) Configuration.Condensed;
    // this is where the mesh is constructed
    
    MHM->BuildComputationalMesh(substructure);
#ifdef PZDEBUG
    if(1)
    {
        std::ofstream file("MHM.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(MHM->GMesh().operator->(), file,true);
    }
#endif
#ifdef PZDEBUG
    if(0)
    {
        std::ofstream out("MHMMeshControl.txt");
        MHM->Print(out);
    }
#endif

    std::cout << "MHM Computational meshes created\n";
#ifdef PZDEBUG
    if(1)
    {

        std::ofstream out_mhm("MHM_hybrid.txt");
        MHM->CMesh()->Print(out_mhm);

    }
#endif
    std::cout << "Number of equations MHM equals " << MHM->CMesh()->NEquations() << std::endl;
    
    std::string configuration;
    
    {
        std::stringstream sout;
        sout << "H" << Configuration.numHDivisions << "-P" << Configuration.pOrderInternal;
        configuration = sout.str();
    }
    if(Configuration.LagrangeMult)
    {
        MHMPref << "_Lagr";
    }
    if (Configuration.Hybridize) {
        MHMPref << "_Hybr";
    }
    // compute the MHM solution
    Configuration.fGlobalSystemWithLocalCondensationSize = MHM->fGlobalSystemWithLocalCondensationSize;
    Configuration.fGlobalSystemSize = MHM->fGlobalSystemSize;
    Configuration.fNumeq = MHM->fNumeq;
    SolveProblem(MHM->CMesh(), MHM->GetMeshes(), example, MHMPref.str(), Configuration);
    
    
    return 0;
}



void InsertMaterialObjects(TPZMHMeshControl &control)
{
    control.fMaterialIds.insert(matInterno);
    control.fMaterialBCIds.insert(bc1);
    control.fMaterialBCIds.insert(bc2);
    control.fMaterialBCIds.insert(bc3);
    control.fMaterialBCIds.insert(bc4);

    TPZCompMesh &cmesh = control.CMesh();
	/// criar materiais
	int dim = cmesh.Dimension();
    TPZMatLaplacianHybrid *material1 = new TPZMatLaplacianHybrid(matInterno,dim);
    
    
    
    material1->SetParameters(1., 0.);
    if(example)
    {
        material1->SetPermeabilityFunction(example->ConstitutiveLawFunction());
        material1->SetForcingFunction(example->ForcingFunction());
    }
    
    
	cmesh.InsertMaterialObject(material1);
	
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,1.), val2(2,1,0.);
	
    //BC -1
    TPZBndCond * BCondD1 = dynamic_cast<TPZBndCond *>( material1->CreateBC(material1, bc1,neumann, val1, val2));
    if (example) {
        BCondD1->SetType(dirichlet);
        BCondD1->TPZDiscontinuousGalerkin::SetForcingFunction(example->ValueFunction());
    }
    //TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet1 = new TPZDummyFunction<REAL>(DirichletValidacao);
    //BCondD1->SetForcingFunction(bcmatDirichlet1);
    cmesh.InsertMaterialObject(BCondD1);
    
    //BC -2
	TPZMaterial * BCondD2 = material1->CreateBC(material1, bc2,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet2 = new TPZDummyFunction<STATE>(DirichletValidacao);
    BCondD2->SetForcingFunction(bcmatDirichlet2);
    if (example) {
        BCondD2->SetForcingFunction(example->ValueFunction());
    }
    cmesh.InsertMaterialObject(BCondD2);
    
    //BC -3
	TPZBndCond* BCondD3 = material1->CreateBC(material1, bc3,neumann, val1, val2);
//    TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet3 = new TPZDummyFunction<STATE>(DirichletValidacao);
//    BCondD3->SetForcingFunction(bcmatDirichlet3);
    if (example) {
        BCondD3->SetType(dirichlet);
        BCondD3->TPZDiscontinuousGalerkin::SetForcingFunction(example->ValueFunction());
    }
    cmesh.InsertMaterialObject(BCondD3);
    
    //BC -4
	TPZMaterial * BCondD4 = material1->CreateBC(material1, bc4,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet4 = new TPZDummyFunction<STATE>(DirichletValidacao);
    BCondD4->SetForcingFunction(bcmatDirichlet4);
    if (example) {
        BCondD4->SetForcingFunction(example->ValueFunction());
    }
    cmesh.InsertMaterialObject(BCondD4);
    
}

void InsertMaterialObjects(TPZMHMixedMeshControl &control)
{
    TPZCompMesh &cmesh = control.CMesh();

    TPZGeoMesh &gmesh = control.GMesh();
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(1,1,0.), val2Flux(1,1,0.), val2Pressure(1,1,10.);
    val2Pressure(0,0) = 1000.;

    int dim = gmesh.Dimension();
    cmesh.SetDimModel(dim);
    
    TPZCompMesh *MixedFluxPressureCmesh = &cmesh;
    
    // Material medio poroso
    TPZMixedPoisson * mat = new TPZMixedPoisson(1,dim);
    mat->SetSymmetric();
    mat->SetPermeability(1.);
    if(!example)
    {
    } else
    {
        mat->SetPermeabilityFunction(example->ConstitutiveLawFunction());
        mat->SetForcingFunction(example->ForcingFunction());
    }
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat);
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, bc1, typePressure, val1, val2Flux);
    TPZAutoPointer<TPZFunction<STATE> > force = new TPZDummyFunction<STATE>(DirichletValidacao);
    bcN->SetForcingFunction(0, force);
    if (example) {
        bcN->SetType(typePressure);
        bcN->TPZMaterial::SetForcingFunction(example->ValueFunction());
    }
    //    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    bcN = mat->CreateBC(mat, bc3, typePressure, val1, val2Flux);
    bcN->SetForcingFunction(0, force);
    if (example) {
        bcN->SetType(typePressure);
        bcN->TPZDiscontinuousGalerkin::SetForcingFunction(example->ValueFunction());
    }
    //    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, bc2, typePressure, val1, val2Pressure);
    bcS->SetForcingFunction(0, force);
    if (example) {
        bcS->TPZMaterial::SetForcingFunction(example->ValueFunction());
    }

    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    bcS = mat->CreateBC(mat, bc4, typePressure, val1, val2Pressure);
    bcS->SetForcingFunction(0, force);
    if (example) {
        bcS->TPZMaterial::SetForcingFunction(example->ValueFunction());
    }
    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    
    TPZMatLaplacian *matdim = new TPZMatLaplacian(1);
    matdim->SetDimension(gmesh.Dimension());
    control.PressureMesh()->InsertMaterialObject(matdim);

    
    control.InsertPeriferalMaterialObjects();
    

}

/// Insert material objects for the MHM-H(div) solution
void InsertMaterialObjects(TPZMHMixedHybridMeshControl &control)
{
    TPZGeoMesh &gmesh = control.GMesh();
    
    int meshdim = gmesh.Dimension();
    control.fFractureFlowDim1MatId.insert(10);
    // Material medio poroso
    TPZMixedPoisson * mat = new TPZMixedPoisson(10,meshdim-1);
    mat->SetSymmetric();
    mat->SetPermeability(1.e-3);
    TPZFNMatrix<9,REAL> K(3,3,0.),KInv(3,3,0.);
    K(0,0) = 1.;
    K(1,1) = 1.e-3;
    KInv(0,0) = 1.;
    KInv(1,1) = 1000.;
    mat->SetPermeabilityTensor(K, KInv);
    
    
    control.CMesh()->InsertMaterialObject(mat);
    
    InsertMaterialObjects((TPZMHMixedMeshControl &) control);
    
}



TPZCompMesh * CreateHDivMHMMesh(TPZGeoMesh * gmesh, int porder)
{
    int meshdim = gmesh->Dimension();
    TPZCompMesh * cmeshHDiv = new TPZCompMesh(gmesh);
    cmeshHDiv->SetDimModel(meshdim);
    cmeshHDiv->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    cmeshHDiv->SetDefaultOrder(porder);
    TPZVecL2 *matl2 = new TPZVecL2(1);
    matl2->SetDimension(2);
    cmeshHDiv->InsertMaterialObject(matl2);
    for (int matid = 2; matid<10; matid++) {
        TPZVecL2 *matl2 = new TPZVecL2(matid);
        matl2->SetDimension(2);
        cmeshHDiv->InsertMaterialObject(matl2);
    }
    TPZFNMatrix<1,STATE> val1(1,1,0.),val2(1,1,0.);
    TPZBndCond *bc = matl2->CreateBC(matl2, -1, 0, val1, val2);
    cmeshHDiv->InsertMaterialObject(bc);
    bc = matl2->CreateBC(matl2, -2, 0, val1, val2);
    cmeshHDiv->InsertMaterialObject(bc);
    cmeshHDiv->AutoBuild();
    
#ifdef PZDEBUG
    {
        std::ofstream outmesh("BigHDivMesh.txt");
        cmeshHDiv->Print(outmesh);
    }
#endif
    return cmeshHDiv;
}


TPZCompMesh * CreatePressureMHMMesh(TPZGeoMesh * gmesh, int porder, int dimension)
{
    TPZCompMesh * cmeshPressure = new TPZCompMesh(gmesh);
    cmeshPressure->SetDimModel(dimension);
    cmeshPressure->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmeshPressure->ApproxSpace().CreateDisconnectedElements(true);
    cmeshPressure->SetDefaultOrder(porder);
    TPZMatLaplacian *matl2 = new TPZMatLaplacian(1);
    matl2->SetDimension(dimension);
    cmeshPressure->InsertMaterialObject(matl2);
    for (int matid = 2; matid<10; matid++) {
        TPZMatLaplacian *matl2 = new TPZMatLaplacian(matid);
        matl2->SetDimension(dimension);
        cmeshPressure->InsertMaterialObject(matl2);
    }

    cmeshPressure->AutoBuild();
    long nc = cmeshPressure->NConnects();
    for (long ic=0; ic<nc; ic++) {
        cmeshPressure->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    return cmeshPressure;
}

TPZCompMesh * CreateHDivPressureMHMMesh(TPZVec<TPZCompMesh * > & cmeshes)
{
    TPZGeoMesh *gmesh = cmeshes[0]->Reference();
    if(!gmesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    int dim = gmesh->Dimension();
    
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(1,1,0.), val2Flux(1,1,0.), val2Pressure(1,1,10.);
    val2Pressure(0,0) = 1000.;
    
    // Malha computacional
    TPZCompMesh * MixedFluxPressureCmesh = new TPZCompMesh(gmesh);
    
    // Material medio poroso
    TPZMixedPoisson * mat = new TPZMixedPoisson(1,dim);
    mat->SetSymmetric();
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat);

    for (int matid = 2; matid<10; matid++)
    {
        // Material medio poroso
        TPZMixedPoisson * mat = new TPZMixedPoisson(matid,dim);
        mat->SetSymmetric();
        //    mat->SetForcingFunction(One);
        MixedFluxPressureCmesh->InsertMaterialObject(mat);
        

    }
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, -1, typeFlux, val1, val2Flux);
    TPZAutoPointer<TPZFunction<STATE> > force = new TPZDummyFunction<STATE>(DirichletValidacao);
    //    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, -2, typePressure, val1, val2Pressure);
    bcS->SetForcingFunction(0, force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    
    
    
    
    
    MixedFluxPressureCmesh->SetDimModel(dim);
    MixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElem();
    MixedFluxPressureCmesh->AutoBuild();
    
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    
    
    meshvector[0] = cmeshes[0];
    meshvector[1] = cmeshes[1];
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, MixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, MixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, MixedFluxPressureCmesh);
    
    return MixedFluxPressureCmesh;

}


TPZAutoPointer<TPZRefPattern> DivideQuadbyTriangles(const std::string refpatname)
{
    TPZGeoMesh gmesh;
    gmesh.NodeVec().Resize(5);
    REAL nodeco[][3] =
    {
        {-1,-1,0},
        {1,-1,0},
        {1,1,0},
        {-1,1,0},
        {0,0,0}
    };
    long nodeindexes[][3] = {
        {0,1,4},
        {1,2,4},
        {2,3,4},
        {3,0,4}
    };
    for (int i=0; i<5; i++) {
        TPZManVector<REAL,3> coord(3);
        for (int c=0; c<3; c++) {
            coord[c] = nodeco[i][c];
        }
        gmesh.NodeVec()[i].Initialize(coord, gmesh);
    }
    TPZManVector<long> corners(4);
    for (long i=0; i<4; i++) {
        corners[i] = i;
    }
    long elindex;
    gmesh.CreateGeoElement(EQuadrilateral, corners, 1, elindex);
    
    long fatherindex = elindex;
    
    for (int is=0; is<4; is++)
    {
        for (long i=0; i<3; i++) {
            corners[i] = nodeindexes[is][i];
        }
        gmesh.CreateGeoElement(ETriangle, corners, 1, elindex);
        gmesh.Element(elindex)->SetFather(fatherindex);
    }
    gmesh.BuildConnectivity();
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(gmesh);
    refpat->SetName(refpatname);
    TPZAutoPointer<TPZRefPattern> found = gRefDBase.FindRefPattern(refpat);
    if(!found)
    {
        gRefDBase.InsertRefPattern(refpat);
        refpat->InsertPermuted();
    }
    else
    {
        refpat = found;
    }
    return refpat;
}


TPZAutoPointer<TPZRefPattern> DivideTriangleby9Triangles(const std::string refpatname)
{
    TPZGeoMesh gmesh;
    gmesh.NodeVec().Resize(10);
    REAL nodeco[][3] =
    {
        {0,0,0}, //0
        {1,0,0}, //1
        {2,0,0},  //2
        {3,0,0},  //3
        {0,1,0},  //4
        {1,1,0},  //5
        {2,1,0},  //6
        {0,2,0},  //7
        {1,2,0},  //8
        {0,3,0} //9
    };
    long nodeindexes[][3] = {
        {0,3,9},
        {0,1,4},
        {1,5,4},
        {1,2,5},
        {2,6,5},
        {2,3,6},
        {4,5,7},
        {5,8,7},
        {5,6,8},
        {7,8,9}
    };
    for (int i=0; i<10; i++) {
        TPZManVector<REAL,3> coord(3);
        for (int c=0; c<3; c++) {
            coord[c] = nodeco[i][c];
        }
        gmesh.NodeVec()[i].Initialize(coord, gmesh);
    }
    TPZManVector<long> corners(3);
    for (long i=0; i<3; i++) {
        corners[i] = nodeindexes[0][i];
    }
    long elindex;
    gmesh.CreateGeoElement(ETriangle, corners, 1, elindex);
    
    long fatherindex = elindex;
    
    for (int is=1; is<10; is++)
    {
        for (long i=0; i<3; i++) {
            corners[i] = nodeindexes[is][i];
        }
        gmesh.CreateGeoElement(ETriangle, corners, 1, elindex);
        gmesh.Element(elindex)->SetFather(fatherindex);
    }
    gmesh.BuildConnectivity();
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(gmesh);
    refpat->SetName(refpatname);
    if(!gRefDBase.FindRefPattern(refpat))
    {
        gRefDBase.InsertRefPattern(refpat);
        refpat->InsertPermuted();
    }
    return refpat;

}

TPZGeoMesh *SimpleMesh(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, TPZVec<long> &coarseindices, int ndiv)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    int dimension = 2;
    gmesh->SetDimension(dimension);
    TPZManVector<int,2> nx(2,3);
    nx[0] = nelx;
    nx[1] = nely;
    TPZGenGrid gengrid(nx, x0, x1);
    gengrid.SetRefpatternElements(true);
    gengrid.Read(gmesh, 1);
    gengrid.SetBC(gmesh, 4, -1);
    gengrid.SetBC(gmesh, 5, -2);
    gengrid.SetBC(gmesh, 6, -3);
    gengrid.SetBC(gmesh, 7, -4);
    
    long nel = gmesh->NElements();
    
    coarseindices.resize(nel);
    long elcount = 0;
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement() ||  gel->Dimension() != dimension) {
            continue;
        }
        coarseindices[elcount] = el;
        elcount++;
    }
    coarseindices.resize(elcount);
#ifdef PZDEBUG
    {
        std::ofstream file("SimpleGMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }
#endif
    
    TPZCheckGeom geom(gmesh);
    geom.UniformRefine(ndiv);
    //    InsertInterfaceElements(gmesh,1,2);
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
#ifdef PZDEBUG
    {
        std::ofstream file("SimpleGMeshDivided.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }
#endif
    
    return gmesh;
}



/// Print the elements with geometric information and connect values
void PrintElements(TPZCompMesh *cmesh, std::ostream &out)
{
    long nelem = cmesh->NElements();
    for (long el = 0; el < nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        if (gel->Dimension() == cmesh->Dimension()) {
            DebugStop();
        }
        TPZManVector<REAL,3> co1(3),co2(3,-100.);
        gel->Node(0).GetCoordinates(co1);
        if(gel->NCornerNodes() > 1) gel->Node(1).GetCoordinates(co2);
        out << "gel index " << gel->Index() << " node loc " << co1 << " and " << co2 << std::endl;
        int nc = cel->NConnects();
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if (c.NShape()) {
                c.Print(*cmesh,out);
            }
        }
    }
}


/// function that randomly refines some elements
void RandomRefine(TPZGeoMesh *gmesh, TPZVec<long> &coarseindices, int nref)
{
    long nel = coarseindices.size();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(coarseindices[el]);
        while (gel->HasSubElement()) {
            int nsub = gel->NSubElements();
            int isub = rand()%nsub;
            gel = gel->SubElement(isub);
        }
        for (int iref = 0; iref<nref; iref++)
        {
            TPZManVector<TPZGeoEl *,10> gelsub;
            gel->Divide(gelsub);
            int nsub = gel->NSubElements();
            int isub = rand()%nsub;
            gel = gel->SubElement(isub);
        }
    }
}
