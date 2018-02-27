#ifdef HAVE_CONFIG_H
#include <pz_config.h>
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
#include "TPZMixedPoissonParabolic.h"



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

TPZGeoMesh *MalhaGeomFred(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, const std::string quad, const std::string triangle, TPZVec<int64_t> &coarseindices, int ndiv);

/// Create a Refinement Pattern that divides a quadrilateral by two triangles
TPZAutoPointer<TPZRefPattern> DivideQuadbyTriangles(const std::string refpatname);

/// Create a Refinement Patterns that divides a triangle into nine triangles
TPZAutoPointer<TPZRefPattern> DivideTriangleby9Triangles(const std::string refpatname);

/// Insert material objects for the MHM Mesh solution
void InsertMaterialObjects(TPZMHMeshControl &control);
/// Insert material objects for the MHM-H(div) solution
void InsertMaterialObjects(TPZMHMixedMeshControl &control);
/// Insert material objects for the MHM-H(div) solution
void InsertMaterialObjects(TPZMHMixedHybridMeshControl &control);

/// Print the elements with geometric information and connect values
void PrintElements(TPZCompMesh *cmesh, std::ostream &out);

/// unwrap de TPZCondensedCompel and TPZElementGroup elements
void UnwrapMesh(TPZCompMesh *cmesh);

/// function that randomly refines some elements
void RandomRefine(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices, int nref);

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
int const bc5=-5;

static void DirichletValidacao(const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &gradres){
    result[0] = loc[0];
}


int main(int argc, char *argv[])
{
    TExceptionManager except;
    
#ifdef _AUTODIFF
    //    example = new TLaplaceExampleSmooth;
#endif
    
    TRunConfig Configuration;
    
    /// computation type :
    // (0) - compute reference mesh
    // (1) - compute MHM H1 mesh and compute MHM(div) mesh
    int ComputationType = 1;
    /// numhdiv - number of h-refinements
    Configuration.numHDivisions = 1;
    /// PolynomialOrder - p-order
    Configuration.pOrderInternal = 2;
    
    
    Configuration.pOrderSkeleton = 1;
    Configuration.numDivSkeleton = 0;
    TPZManVector<REAL,3> x0(2,0.),x1(2,0.);
    // for using the aligned mesh
    x0[0] = 1.;
    int nelxref = 2;
    int nelyref = 2;
    Configuration.nelxcoarse = nelxref;
    Configuration.nelycoarse = nelyref;
    Configuration.Hybridize = 1;
    Configuration.Condensed = 0;

    if (argc == 8)
    {
        std::cout << "Executing using command line arguments\n";
        Configuration.nelxcoarse = atoi(argv[1]);
        Configuration.nelycoarse = atoi(argv[2]);
        Configuration.numHDivisions = atoi(argv[3]);
        Configuration.pOrderInternal = atoi(argv[4]);
        Configuration.numDivSkeleton = atoi(argv[5]);
        Configuration.pOrderSkeleton = atoi(argv[6]);
        Configuration.newline = atoi(argv[7]);
    }
    
    
    x0.Fill(0.);
    x1.Fill(1.);

    HDivPiola = 1;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    // tototo
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
//    gRefDBase.InitializeUniformRefPattern(ECube);

//    gRefDBase.InitializeRefPatterns(2);
    
    TPZGeoMesh *gmesh = 0;
    TPZManVector<int64_t> coarseindices;
    // verifying differences between the MHM-original and MHM with mixed approximations
    int nelx = Configuration.nelxcoarse;
    int nely = Configuration.nelycoarse;
    {
        std::ofstream out("DiffResults.nb",std::ios::app);
        out << "(* Running quadrilateral mesh with numsubdomains " << nelx << ", " << nely << " *)\n";
    }
    /// Analise the regularity of the subdomain problems
    TPZManVector<int,3> nelvec(2),nsub(2);
    nelvec[0] = Configuration.nelxcoarse;
    nelvec[1] = Configuration.nelycoarse;
    nsub[0] = nelx;
    nsub[1] = nely;
    int ndiv = Configuration.numHDivisions;
    gmesh = MalhaGeomFredQuadrada(nelx, nely, x0, x1, coarseindices, ndiv);
    RandomRefine(gmesh, coarseindices,1);
    
    TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
    TPZAutoPointer<TPZMHMixedHybridMeshControl> MHM;
    TPZAutoPointer<TPZMHMixedMeshControl> MHMixed;
    
    std::stringstream MHMPref;


    if(1)
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmesh);
        TPZMHMixedHybridMeshControl *mhm = new TPZMHMixedHybridMeshControl(gmeshauto);
        mhm->DefinePartitionbyCoarseIndices(coarseindices);
        MHMPref << "MHMixedHybrid";
        MHM = mhm;
        TPZMHMeshControl &meshcontrol = *mhm;
        MHM->SwitchLagrangeMultiplierSign(true);

        if (Configuration.LagrangeMult) {
            meshcontrol.SetLagrangeAveragePressure(true);
        }
        
        InsertMaterialObjects(*mhm);

        meshcontrol.SetInternalPOrder(Configuration.pOrderInternal);
        meshcontrol.SetSkeletonPOrder(Configuration.pOrderSkeleton);
        
        meshcontrol.DivideSkeletonElements(Configuration.numDivSkeleton);
        if (Configuration.Hybridize)
        {
            meshcontrol.Hybridize();
        }
        
        bool substructure = (bool) Configuration.Condensed;
        meshcontrol.BuildComputationalMesh(substructure);
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream file("GMeshControl.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file,true);
        }
#endif
#ifdef PZDEBUG
        if(0)
        {
            std::ofstream out("MHMMeshControl.txt");
            meshcontrol.Print(out);
        }
#endif

        std::cout << "MHM Computational meshes created\n";
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream gfile("geometry.txt");
            gmesh->Print(gfile);

            std::ofstream out_mhm("MHM_hybrid.txt");
            meshcontrol.CMesh()->Print(out_mhm);

        }
#endif
        std::cout << "Number of equations MHM equals " << MHM->CMesh()->NEquations() << std::endl;
    
    }
    
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
    SolveProblem(MHM->CMesh(), MHM->GetMeshes(), 0, MHMPref.str(), Configuration);
    

    return 0;
}



void InsertMaterialObjects(TPZMHMeshControl &control)
{
    TPZCompMesh &cmesh = control.CMesh();
	/// criar materiais
	int dim = cmesh.Dimension();
    TPZMatLaplacianHybrid *material1 = new TPZMatLaplacianHybrid(matInterno,dim);
    
    material1->SetParameters(1., 0.);
    
	TPZMaterial * mat1(material1);
    
    TPZMat1dLin *materialCoarse = new TPZMat1dLin(matCoarse);
    TPZFNMatrix<1,STATE> xk(1,1,0.),xb(1,1,0.),xc(1,1,0.),xf(1,1,0.);
    materialCoarse->SetMaterial(xk, xc, xb, xf);
    
    cmesh.InsertMaterialObject(materialCoarse);
    materialCoarse = new TPZMat1dLin(skeleton);
    materialCoarse->SetMaterial(xk, xc, xb, xf);
    cmesh.InsertMaterialObject(materialCoarse);
    materialCoarse = new TPZMat1dLin(secondskeleton);
    materialCoarse->SetMaterial(xk, xc, xb, xf);
    cmesh.InsertMaterialObject(materialCoarse);
    materialCoarse = new TPZMat1dLin(matpressure);
    materialCoarse->SetMaterial(xk, xc, xb, xf);
    cmesh.InsertMaterialObject(materialCoarse);
    
    
    
    //    REAL diff = -1.;
    //	REAL conv = 0.;
    //	TPZVec<REAL> convdir(3,0.);
    //	REAL flux = 8.;
    //
    //	material1->SetParameters(diff, conv, convdir);
    //	material1->SetInternalFlux( flux);
    
	cmesh.InsertMaterialObject(mat1);
	
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,1.), val2(2,1,0.);
	
    //BC -1
    TPZBndCond * BCondD1 = dynamic_cast<TPZBndCond *>( material1->CreateBC(mat1, bc1,neumann, val1, val2));
    //TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet1 = new TPZDummyFunction<REAL>(DirichletValidacao);
    //BCondD1->SetForcingFunction(bcmatDirichlet1);
    cmesh.InsertMaterialObject(BCondD1);
    
    //BC -2
	TPZMaterial * BCondD2 = material1->CreateBC(mat1, bc2,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet2 = new TPZDummyFunction<STATE>(DirichletValidacao);
    BCondD2->SetForcingFunction(bcmatDirichlet2);
    cmesh.InsertMaterialObject(BCondD2);
    
    //BC -3
	TPZBndCond* BCondD3 = material1->CreateBC(mat1, bc3,neumann, val1, val2);
//    TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet3 = new TPZDummyFunction<STATE>(DirichletValidacao);
//    BCondD3->SetForcingFunction(bcmatDirichlet3);
    cmesh.InsertMaterialObject(BCondD3);
    
    //BC -4
	TPZMaterial * BCondD4 = material1->CreateBC(mat1, bc4,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet4 = new TPZDummyFunction<STATE>(DirichletValidacao);
    BCondD4->SetForcingFunction(bcmatDirichlet4);
    cmesh.InsertMaterialObject(BCondD4);
    
    //BC -5: dirichlet nulo
    TPZMaterial * BCondD5 = material1->CreateBC(mat1, bc5,dirichlet, val1, val2);
    cmesh.InsertMaterialObject(BCondD5);
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
    int matid = 1;
    int dimension = 2;
    //    TPZMixedPoisson * mat = new TPZMixedPoisson(1,dim);
    TPZMixedPoissonParabolic *mat = new TPZMixedPoissonParabolic(matid,dimension);
    mat->SetDeltaT(0.0001);
    mat->SetSymmetric();
    mat->SetPermeability(1.);
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat);
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, -1, typeFlux, val1, val2Flux);
    //    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    bcN = mat->CreateBC(mat, -3, typeFlux, val1, val2Flux);
    //    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, -2, typePressure, val1, val2Pressure);

    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    bcS = mat->CreateBC(mat, -4, typeFlux, val1, val2Pressure);
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
    control.InsertFractureFlowMaterial(10);
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




/// Create a reference geometric mesh starting with nelx by nely domains
TPZGeoMesh *CreateReferenceGMesh(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, int numref)
{
    TPZManVector<int,3> nx(2);
    nx[0] = nelx;
    nx[1] = nely;
    TPZGenGrid gengrid(nx,x0,x1);
    gengrid.SetRefpatternElements(true);
    TPZGeoMesh *result = new TPZGeoMesh;
    int matid = 1;
    gengrid.Read(result, matid);
    gengrid.SetBC(result, 4, -1);
    gengrid.SetBC(result, 5, -2);
    gengrid.SetBC(result, 6, -1);
    gengrid.SetBC(result, 7, -2);
    int matidpoint = 10;
    
    int64_t firstnode = nelx+2;
    int64_t numnodes = nelx-1;
    for (int64_t ynode = 1; ynode < nely; ynode++)
    {
        for (int64_t node = firstnode; node < firstnode+numnodes; node++) {
            TPZManVector<int64_t,2> nodeindices(1);
            nodeindices[0] = node;
            int64_t index;
            result->CreateGeoElement(EPoint, nodeindices, matidpoint, index);
        }
        firstnode += nelx+1;
    }
    result->BuildConnectivity();
    // refina a malha uma vez uniformemente
    int numuni = 1;
    for (int uni=0; uni<numuni; uni++)
    {
        int64_t nelem = result->NElements();
        for (int64_t el=0; el<nelem; el++) {
            TPZGeoEl *gel = result->Element(el);
            if (gel->Dimension() == 0) {
                continue;
            }
            TPZManVector<TPZGeoEl *,8> subs;
            gel->Divide(subs);
        }
    }
    // refina a malha na direcao dos elementos ponto
    std::set<int> matids;
    matids.insert(matidpoint);
    for (int cycle = 0; cycle < numref; cycle++) {
        int64_t nelem = result->NElements();
        for (int64_t el=0; el<nelem; el++) {
            TPZGeoEl *gel = result->Element(el);
            int targetmatid = cycle+2;
            TPZRefPatternTools::RefineDirectional(gel, matids, targetmatid);
        }
    }
    
    {
        std::ofstream out("../ReferenceGMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(result, out);
    }
    return result;
}


void UnwrapMesh(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    bool change = true;
    while(change)
    {
        change = false;
        for (int64_t el=0; el<nel; el++) {
            
            TPZCompEl *cel = cmesh->Element(el);
            TPZCondensedCompEl *condense = dynamic_cast<TPZCondensedCompEl *>(cel);
            if (condense) {
                condense->Unwrap();
                change = true;
            }
            cel = cmesh->Element(el);
            TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
            if (elgr) {
                elgr->Unwrap();
                change = true;
            }
        }
    }
}


/// Print the elements with geometric information and connect values
void PrintElements(TPZCompMesh *cmesh, std::ostream &out)
{
    int64_t nelem = cmesh->NElements();
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        if (gel->Dimension() != 1) {
            DebugStop();
        }
        TPZManVector<REAL,3> co1(3),co2(3);
        gel->Node(0).GetCoordinates(co1);
        gel->Node(1).GetCoordinates(co2);
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
void RandomRefine(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices, int nref)
{
    int64_t nel = coarseindices.size();
    for (int64_t el=0; el<nel; el++) {
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
