#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "TPZBFileStream.h"

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

#include "meshgen.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>



using namespace std;

struct TRunConfig;

TPZGeoMesh *MalhaGeomFred(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, const std::string quad, const std::string triangle, TPZVec<long> &coarseindices, int ndiv);

/// Create a Refinement Pattern that divides a quadrilateral by two triangles
TPZAutoPointer<TPZRefPattern> DivideQuadbyTriangles(const std::string refpatname);

/// Create a Refinement Patterns that divides a triangle into nine triangles
TPZAutoPointer<TPZRefPattern> DivideTriangleby9Triangles(const std::string refpatname);

/// Insert material objects for the MHM Mesh solution
void InsertMaterialObjects(TPZMHMeshControl &control);
/// Insert material objects for the MHM-H(div) solution
void InsertMaterialObjects(TPZMHMixedMeshControl &control);

/// Compute the differences at the submesh level
void ComputeDifferencesBySubmesh(TRunConfig &config, TPZMHMeshControl &MHM, TPZMHMixedMeshControl &MHMixed, const std::string &filename);

/// Create a reference geometric mesh starting with nelx by nely domains
// called in the main program
TPZGeoMesh *CreateReferenceGMesh(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, int numref);

/// compute the reference solution and return created mesh
// call CreateReferenceCMesh
TPZCompMesh *ComputeReferenceSolution(TPZGeoMesh *gmesh, int porder, TPZVec<TPZCompMesh *> &meshvec);

/// create the computational mesh of the reference solution
// call CreateHDivMHMMesh and CreatePressureMHMMesh
TPZCompMesh *CreateReferenceCMesh(TPZGeoMesh *gmesh, TPZVec<TPZCompMesh *> &meshvec, int porder);

/// Create an HDiv mesh used as a reference mesh
TPZCompMesh * CreateHDivMHMMesh(TPZGeoMesh * gmesh, int porder);
/// Create a pressure mesh for the reference mesh
TPZCompMesh * CreatePressureMHMMesh(TPZGeoMesh * gmesh, int porder, int dimension);

/// Create a multiphysics mesh from an H(div) and Pressure mesh
/// Called by the method which creates the reference computacional mesh
TPZCompMesh * CreateHDivPressureMHMMesh(TPZVec<TPZCompMesh * > &cmesh);


/// Analise the regularity of the subdomain problems
void AnalyseRegularity(const TPZVec<int> &pos0,const TPZVec<int> &nelx, TPZVec<int> &nsub, TPZFMatrix<REAL> &lowestexp);

/// Print the elements with geometric information and connect values
void PrintElements(TPZCompMesh *cmesh, std::ostream &out);

/// copy the solution between one computation mesh to the other assuming the geometric elements match
void CopySolution(TPZCompMesh *from, TPZCompMesh *to);

/// unwrap de TPZCondensedCompel and TPZElementGroup elements
void UnwrapMesh(TPZCompMesh *cmesh);

/// function that returns the permeability for a given coordinate
void Permeability(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &diff);

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mainskeleton"));
#endif


const int matInterno = 1;
const int matCoarse = 2;
const int skeleton = 4;
const int secondskeleton = 3;
const int matInterface = 5;
const int matpressure = 6;

const int dirichlet = 0;
const int neumann = 1;
const int mixed = 2;

int const bc1=-1;
int const bc2=-2;
int const bc3=-3;
int const bc4=-4;
int const bc5=-5;

static void DirichletValidacao(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    result[0] = loc[0];
}

TPZFMatrix<REAL> gPorous(500,100,0.);

TAnalyticSolution *example = 0;

int main(int argc, char *argv[])
{
    TExceptionManager except;
    
#ifdef _AUTODIFF
    example = new TLaplaceExample1;
#endif
    
    TRunConfig Configuration;
    
    /// computation type :
    // (0) - compute reference mesh
    // (1) - compute MHM H1 mesh and compute MHM(div) mesh
    int ComputationType = 1;
    /// numhdiv - number of h-refinements
    Configuration.numHDivisions = 3;
    /// PolynomialOrder - p-order
    Configuration.pOrderInternal = 2;
    
    
    Configuration.pOrderSkeleton = 1;
    Configuration.numDivSkeleton = 0;
    TPZManVector<REAL,3> x0(2,0.),x1(2,0.);
    // for using the aligned mesh
    x0[0] = 1.;
    if (!example)
    {
        int nelxref = 64;
        int nelyref = 16;
        Configuration.nelxcoarse = nelxref;
        Configuration.nelycoarse = nelyref;
    }
    else
    {
        Configuration.nelxcoarse = 8;
        Configuration.nelycoarse = 8;
    }
    Configuration.Hybridize = 0;
    if (argc == 7)
    {
        std::cout << "Executing using command line arguments\n";
        Configuration.nelxcoarse = atoi(argv[1]);
        Configuration.nelycoarse = atoi(argv[2]);
        Configuration.numHDivisions = atoi(argv[3]);
        Configuration.pOrderInternal = atoi(argv[4]);
        Configuration.numDivSkeleton = atoi(argv[5]);
        Configuration.pOrderSkeleton = atoi(argv[6]);
    }
    
    // to avoid singular internal matrices
    if (Configuration.numHDivisions == 0 && Configuration.pOrderInternal <= Configuration.pOrderSkeleton) {
        Configuration.pOrderInternal = Configuration.pOrderSkeleton+1;
    }
    x1[0] = x0[0]+0.01*Configuration.nelxcoarse;
    x1[1] = x0[1]+0.01*Configuration.nelycoarse;
    
    if(example)
    {
        x0.Fill(0.);
        x1.Fill(1.);
    }

    HDivPiola = 1;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    {
        std::ifstream pores("../porous.txt");
        for (int j=0; j<100; j++) {
            for (int i=0; i<500; i++) {
                pores >> gPorous(i,j);
                if (!pores) {
                    DebugStop();
                }
            }
        }
    }
    // tototo
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
//    gRefDBase.InitializeUniformRefPattern(ECube);
    
    TPZManVector<TPZCompMesh *,2> ReferenceMeshVec(2,0);
    TPZGeoMesh *ReferenceGMesh = 0;
    TPZCompMesh *ReferenceCMesh = 0;

//    gRefDBase.InitializeRefPatterns();
    
    if(ComputationType == 0)
    {
        // generate the reference solution, save it on disk and exit
        int nelx = Configuration.nelxcoarse, nely = Configuration.nelycoarse;
        int numref = 1;
        TPZGeoMesh *gmesh = CreateReferenceGMesh(nelx, nely, x0, x1, numref);
        TPZManVector<TPZCompMesh *,2> meshvec(2);
        int porder = 1;
        TPZCompMesh *cmesh = ComputeReferenceSolution(gmesh,porder,meshvec);
        TPZBFileStream meshfile;
        meshfile.OpenWrite("Ref.bin");
        gmesh->ResetReference();
        gmesh->Write(meshfile, false);
        meshvec[0]->Write(meshfile, false);
        meshvec[1]->Write(meshfile, false);
        //       cmesh->Write(meshfile, false);
        if(0)
        {
            ofstream out("gmesh1.txt");
            gmesh->Print(out);
        }
        if(0)
        {
            ofstream out1("cmeshwrite0.txt");
            meshvec[0]->Print(out1);
            ofstream out2("cmeshwrite1.txt");
            meshvec[1]->Print(out2);
        }
        UnwrapMesh(cmesh);

        if(0)
        {
            ofstream out1("mfmeshwrite.txt");
            cmesh->Print(out1);
        }

        delete cmesh;
        delete meshvec[0];
        delete meshvec[1];
        delete gmesh;
        exit(0);
    }
    if(0)
    {
        // read the reference solution from the file
        TPZBFileStream meshfile;
        meshfile.OpenRead("Ref.bin");
        TPZGeoMesh *gmesh = new TPZGeoMesh;
        gmesh->Read(meshfile, 0);
        ReferenceGMesh = gmesh;
        if(0)
        {
            ofstream out("gmesh2.txt");
            gmesh->Print(out);
        }
        TPZManVector<TPZCompMesh *,2> meshvec(2);
        meshvec[0] = new TPZCompMesh;
        meshvec[1] = new TPZCompMesh;
        meshvec[0]->Read(meshfile, gmesh);
        meshvec[1]->Read(meshfile, gmesh);
        ReferenceMeshVec = meshvec;
        if(0)
        {
            ofstream out1("cmeshread0.txt");
            meshvec[0]->Print(out1);
            ofstream out2("cmeshread1.txt");
            meshvec[1]->Print(out2);
        }
        TPZCompMesh *cmesh = CreateHDivPressureMHMMesh(meshvec);
        for (int i=0; i<10; i++) {
            TPZMaterial *mat = cmesh->FindMaterial(i);
            if (mat) {
                TPZMixedPoisson *mixed = dynamic_cast<TPZMixedPoisson *>(mat);
                if (!mixed) {
                    DebugStop();
                }
                TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Permeability);
                dummy->SetPolynomialOrder(0);
                TPZAutoPointer<TPZFunction<STATE> > func(dummy);
                mixed->SetPermeabilityFunction(func);
            }
        }

        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, cmesh);
        ReferenceCMesh = cmesh;
        if(0)
        {
            std::string plotfile("referencesolutionRead.vtk");
            TPZStack<std::string> scalnames,vecnames;
            scalnames.Push("Pressure");
            scalnames.Push("Permeability");
            vecnames.Push("Derivative");
            vecnames.Push("Flux");
            TPZAnalysis an(cmesh,false);
            if(0)
            {
                ofstream out1("mfmeshread.txt");
                cmesh->Print(out1);
            }
            an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotfile);
            int resolution = 0;
            an.PostProcess(resolution,cmesh->Dimension());
        }
    }
    TPZGeoMesh *gmesh = 0;
    TPZVec<long> coarseindices;
    if(0)
    {
        // original research paper - the mesh was not aligned with the heterogeneities
        std::string quad = "QuadByTriangles";
        std::string triangle = "TriangleBy9Triangles";
        TPZAutoPointer<TPZRefPattern> refpatquad = DivideQuadbyTriangles(quad);
        quad = refpatquad->Name();
        TPZAutoPointer<TPZRefPattern> refpattriangle = DivideTriangleby9Triangles(triangle);
        int nelx = 15;
        int nely = 5;
        Configuration.nelxcoarse = nelx;
        Configuration.nelycoarse = nely;
        int ndiv = Configuration.numHDivisions;
        gmesh = MalhaGeomFred(nelx, nely, x0, x1, quad, triangle, coarseindices, ndiv);
        {
            std::ofstream out("DiffResults.nb",std::ios::app);
            out << "(* Running triangular mesh with subdomains " << nelx << " " << nely << " *)\n";
        }
    }
    else if(!example)
    {
        // verifying differences between the MHM-original and MHM with mixed approximations
        int nelx = 16;
        int nely = 4;
        Configuration.nelxcoarse = nelx;
        Configuration.nelycoarse = nely;
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
        TPZFMatrix<REAL> lowestexp;
        TPZManVector<int,2> pos0(2,0);
        pos0[0] = 100;

        AnalyseRegularity(pos0, nelvec,  nsub,  lowestexp);

        VisualMatrixVTK(lowestexp, "regularity.vtk");
        {
            std::ofstream out("regularity.nb");
            lowestexp.Print("Regularity=",out,EMathematicaInput);
        }

        int ndiv = Configuration.numHDivisions;
        gmesh = MalhaGeomFredQuadrada(nelx, nely, x0, x1, coarseindices, ndiv);
    }
    else
    {
        {
            std::ofstream out("DiffResults.nb",std::ios::app);
            out << "(* Running quadrilateral mesh with Config { ";
            Configuration.MathematicaInlinePrint(out);
            out << "} *)\n";
        }
        int ndiv = Configuration.numHDivisions;
        gmesh = MalhaGeomFredQuadrada(Configuration.nelxcoarse, Configuration.nelycoarse, x0, x1, coarseindices, ndiv);
        
    }
    
    TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
    TPZAutoPointer<TPZMHMeshControl> MHM;
    TPZAutoPointer<TPZMHMixedMeshControl> MHMixed;

    if (1)
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmesh);
        TPZMHMeshControl *mhm = new TPZMHMeshControl(gmeshauto,coarseindices);
        MHM = mhm;
        TPZMHMeshControl &meshcontrol = *mhm;
        MHM->SwitchLagrangeMultiplierSign(true);

        meshcontrol.SetLagrangeAveragePressure(true);
        
        InsertMaterialObjects(meshcontrol);

        meshcontrol.SetInternalPOrder(Configuration.pOrderInternal);
        meshcontrol.SetSkeletonPOrder(Configuration.pOrderSkeleton);
        
        meshcontrol.CreateSkeletonElements(skeleton);
        
        meshcontrol.DivideSkeletonElements(Configuration.numDivSkeleton);
        if (Configuration.Hybridize)
        {
            meshcontrol.Hybridize(secondskeleton, matpressure);
        }
        
        bool substructure = true;
        meshcontrol.BuildComputationalMesh(substructure);
#ifdef PZDEBUG
        if(0)
        {
            std::ofstream file("GMeshControl.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file);
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
        if(0)
        {
            std::ofstream gfile("geometry.txt");
            gmesh->Print(gfile);

            std::ofstream out_mhm("MHM_h1.txt");
            meshcontrol.CMesh()->Print(out_mhm);

        }
#endif
        std::cout << "Number of equations MHM equals " << MHM->CMesh()->NEquations() << std::endl;
    
    }
    
    if(1)
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmesh);
        TPZMHMixedMeshControl *mhm = new TPZMHMixedMeshControl(gmeshauto,coarseindices);
        MHMixed = mhm;
        TPZMHMixedMeshControl &meshcontrol = *mhm;
        
        
        InsertMaterialObjects(meshcontrol);
        
        meshcontrol.SetInternalPOrder(Configuration.pOrderInternal);
        meshcontrol.SetSkeletonPOrder(Configuration.pOrderSkeleton);
        
        // criam-se apenas elementos geometricos
        int matskeleton = skeleton;
        meshcontrol.CreateSkeletonElements(matskeleton);
        meshcontrol.DivideSkeletonElements(Configuration.numDivSkeleton);

        if (Configuration.Hybridize)
        {
            meshcontrol.TPZMHMeshControl::Hybridize(secondskeleton, matpressure);
            
        }
        
        bool substructure = true;
        meshcontrol.BuildComputationalMesh(substructure);
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream file("GMeshControlHDiv.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file);
        }
#endif
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream out("MixedMeshControlHDiv.txt");
            meshcontrol.Print(out);
        }
#endif
        
        std::cout << "MHM Hdiv Computational meshes created\n";
#ifdef PZDEBUG
        if(0)
        {
            std::ofstream gfile("geometryMHMHdiv.txt");
            gmeshauto->Print(gfile);
            std::ofstream out_mhm("MHM_hdiv.txt");
            meshcontrol.CMesh()->Print(out_mhm);
            
        }
#endif
        
        std::cout << "Number of equations MHMixed " << MHMixed->CMesh()->NEquations() << std::endl;

    }
    else
    {
        DebugStop();
    }
    std::string configuration;
    
    {
        std::stringstream sout;
        sout << "H" << Configuration.numHDivisions << "-P" << Configuration.pOrderInternal;
        configuration = sout.str();
    }
    // compute the MHM solution
    SolveProblem(MHM->CMesh(), MHM->GetMeshes(), example, "MHM", Configuration);
    
    // compute the MHM H(div) solution
    SolveProblem(MHMixed->CMesh(), MHMixed->GetMeshes(), example, "MHMixed", Configuration);
    
//    CopySolution(MHMixed->CMesh().operator->(), MHM->CMesh().operator->());
    
    {
        std::string filename = "MHMixed_" + configuration + ".txt";
        std::ofstream out(filename);
        PrintElements(MHMixed->CMesh().operator->(), out);
    }
    {
        std::string filename = "MHM_" + configuration + ".txt";
        std::ofstream out(filename);
        PrintElements(MHM->CMesh().operator->(), out);
    }
    
//    ComputeDifferencesBySubmesh(Configuration, MHM, MHMixed, "DiffResults.nb");
    if(0)
    {
        TPZManVector<STATE,10> square_errors(3,0.);
        TPZCompMeshTools::ComputeDifferenceNorm(MHMixed->CMesh().operator->(), MHM->CMesh().operator->(), square_errors);
        std::cout << "Difference between both formulations " << square_errors << std::endl;
        {
            std::ofstream out("DiffResults.nb",std::ios::app);
            out << "(* domain size " << Configuration.nelxcoarse << " " << Configuration.nelycoarse << " num subdomains " << MHM->Coarse_to_Submesh().size() << " *)\n";
            out << "AppendTo[results, {";
            out << " ";
            Configuration.MathematicaInlinePrint(out);
            out << " ,";
            out << " {";
            out << square_errors;
            out << " } }];\n";
        }
    }

    return 0;
}



void InsertMaterialObjects(TPZMHMeshControl &control)
{
    TPZCompMesh &cmesh = control.CMesh();
	/// criar materiais
	int dim = cmesh.Dimension();
    TPZMatLaplacianHybrid *material1 = new TPZMatLaplacianHybrid(matInterno,dim);
    
    material1->SetParameters(1., 0.);
    if(!example)
    {
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Permeability);
        dummy->SetPolynomialOrder(0);
        TPZAutoPointer<TPZFunction<STATE> > func(dummy);
        material1->SetPermeabilityFunction(func);
    }
    else
    {
        material1->SetPermeabilityFunction(example->ConstitutiveLawFunction());
        material1->SetForcingFunction(example->ForcingFunction());
    }
    
    
    
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
    if (example) {
        BCondD1->SetType(dirichlet);
        BCondD1->TPZDiscontinuousGalerkin::SetForcingFunction(example->ValueFunction());
    }
    //TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet1 = new TPZDummyFunction<REAL>(DirichletValidacao);
    //BCondD1->SetForcingFunction(bcmatDirichlet1);
    cmesh.InsertMaterialObject(BCondD1);
    
    //BC -2
	TPZMaterial * BCondD2 = material1->CreateBC(mat1, bc2,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet2 = new TPZDummyFunction<STATE>(DirichletValidacao);
    BCondD2->SetForcingFunction(bcmatDirichlet2);
    if (example) {
        BCondD2->SetForcingFunction(example->ValueFunction());
    }
    cmesh.InsertMaterialObject(BCondD2);
    
    //BC -3
	TPZBndCond* BCondD3 = material1->CreateBC(mat1, bc3,neumann, val1, val2);
//    TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet3 = new TPZDummyFunction<STATE>(DirichletValidacao);
//    BCondD3->SetForcingFunction(bcmatDirichlet3);
    if (example) {
        BCondD3->SetType(dirichlet);
        BCondD3->TPZDiscontinuousGalerkin::SetForcingFunction(example->ValueFunction());
    }
    cmesh.InsertMaterialObject(BCondD3);
    
    //BC -4
	TPZMaterial * BCondD4 = material1->CreateBC(mat1, bc4,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet4 = new TPZDummyFunction<STATE>(DirichletValidacao);
    BCondD4->SetForcingFunction(bcmatDirichlet4);
    if (example) {
        BCondD4->SetForcingFunction(example->ValueFunction());
    }
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
    TPZMixedPoisson * mat = new TPZMixedPoisson(1,dim);
    mat->SetSymmetric();
    mat->SetPermeability(1.);
    if(!example)
    {
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Permeability);
        dummy->SetPolynomialOrder(0);
        TPZAutoPointer<TPZFunction<STATE> > func(dummy);
        mat->SetPermeabilityFunction(func);
    } else
    {
        mat->SetPermeabilityFunction(example->ConstitutiveLawFunction());
        mat->SetForcingFunction(example->ForcingFunction());
    }
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat);
    int LagrangeMatIdLeft = 50;
    int LagrangeMatIdRight = 51;
    int nstate = 1;
    TPZLagrangeMultiplier *matleft = new TPZLagrangeMultiplier(LagrangeMatIdLeft,dim,nstate);
    TPZLagrangeMultiplier *matright = new TPZLagrangeMultiplier(LagrangeMatIdRight,dim,nstate);

    MixedFluxPressureCmesh->InsertMaterialObject(matleft);
    MixedFluxPressureCmesh->InsertMaterialObject(matright);
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, -1, typeFlux, val1, val2Flux);
    TPZAutoPointer<TPZFunction<STATE> > force = new TPZDummyFunction<STATE>(DirichletValidacao);
    if (example) {
        bcN->SetType(typePressure);
        bcN->TPZMaterial::SetForcingFunction(example->ValueFunction());
    }
    //    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    bcN = mat->CreateBC(mat, -3, typeFlux, val1, val2Flux);
    if (example) {
        bcN->SetType(typePressure);
        bcN->TPZDiscontinuousGalerkin::SetForcingFunction(example->ValueFunction());
    }
    //    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
    TPZBndCond * bcPressure = mat->CreateBC(mat, matpressure, typePressure, val1, val2Flux);
    //    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcPressure);
    
    TPZBndCond * bcFlux = mat->CreateBC(mat, skeleton, typePressure, val1, val2Flux);
    //    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcFlux);
    
    TPZBndCond * bcSecondFlux = mat->CreateBC(mat, secondskeleton, typePressure, val1, val2Flux);
    //    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcSecondFlux);
    
    

    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, -2, typePressure, val1, val2Pressure);
    bcS->SetForcingFunction(0, force);
    if (example) {
        bcS->TPZMaterial::SetForcingFunction(example->ValueFunction());
    }

    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    bcS = mat->CreateBC(mat, -4, typePressure, val1, val2Pressure);
    bcS->SetForcingFunction(0, force);
    if (example) {
        bcS->TPZMaterial::SetForcingFunction(example->ValueFunction());
    }
    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    
    // Bc N
    val2Pressure.Zero();
    TPZBndCond * skelmat = mat->CreateBC(mat, skeleton, typePressure, val1, val2Pressure);
//    TPZAutoPointer<TPZFunction<STATE> > force = new TPZDummyFunction<STATE>(DirichletValidacao);
    //    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(skelmat);
    
    

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

TPZGeoMesh *MalhaGeomFred(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, const std::string quad, const std::string triangle, TPZVec<long> &coarseindices, int ndiv)
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
    gengrid.SetBC(gmesh, 6, -1);
    gengrid.SetBC(gmesh, 7, -2);
    
    TPZAutoPointer<TPZRefPattern> refquad,reftriangle;
    refquad = gRefDBase.FindRefPattern(quad);
    reftriangle = gRefDBase.FindRefPattern(triangle);
    if (!refquad || ! reftriangle) {
        DebugStop();
    }
    long nel = gmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->Type() == EQuadrilateral   ) {
            gel->SetRefPattern(refquad);
            TPZManVector<TPZGeoEl *,4> subs;
            gel->Divide(subs);
        }
    }
    nel = gmesh->NElements();

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
        std::ofstream file("GMeshFredCoarse.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }
#endif

    if(1)
    {
    
        for (long el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel->HasSubElement() &&  gel->Type() == ETriangle) {
                gel->SetRefPattern(reftriangle);
                TPZManVector<TPZGeoEl *,12> subs;
                gel->Divide(subs);
            }
        }
        nel = gmesh->NElements();
        
        for (long el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel->HasSubElement() &&  gel->Type() == EOned) {
                TPZAutoPointer<TPZRefPattern> refpat = TPZRefPatternTools::PerfectMatchRefPattern(gel);
                if (!refpat) {
                    DebugStop();
                }
                gel->SetRefPattern(refpat);
                TPZManVector<TPZGeoEl *,12> subs;
                gel->Divide(subs);
            }
        }
    }
    
#ifdef PZDEBUG
    {
        std::ofstream file("GMeshFredInit.vtk");
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
        std::ofstream file("GMeshFred.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }
#endif
    
    return gmesh;
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
    
    long firstnode = nelx+2;
    long numnodes = nelx-1;
    for (long ynode = 1; ynode < nely; ynode++)
    {
        for (long node = firstnode; node < firstnode+numnodes; node++) {
            TPZManVector<long,2> nodeindices(1);
            nodeindices[0] = node;
            long index;
            result->CreateGeoElement(EPoint, nodeindices, matidpoint, index);
        }
        firstnode += nelx+1;
    }
    result->BuildConnectivity();
    // refina a malha uma vez uniformemente
    int numuni = 1;
    for (int uni=0; uni<numuni; uni++)
    {
        long nelem = result->NElements();
        for (long el=0; el<nelem; el++) {
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
        long nelem = result->NElements();
        for (long el=0; el<nelem; el++) {
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

void Permeability(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &diff)
{
    long ix = x[0]*100;
    long iy = x[1]*100;
    static int count = 0;
    if((fabs(ix-x[0]*100) < 1.e-6 || fabs(ix-x[1]*100) < 1.e-6) && count < 20)
    {
        count++;
        std::cout << "probing for a permeability at the interface of two regions\n";
        std::cout << "x = " << x << std::endl;
    }
    if (IsZero(x[1]-1.)) {
        iy = 99;
    }
    if (IsZero(x[0]-5.)) {
        ix = 499;
    }
    REAL valporous = gPorous(ix,iy);
    // totototototo
    //    valporous = 1.+0.3*sin(x[0]*50)*cos(x[1]*50.);
    for (int i=0; i<2; i++) {
        diff(i,i) = valporous;
        diff(i,1-i)=0.;
        diff(2+i,i) = 1./valporous;
        diff(2+i,1-i) = 0.;
    }
//    for (int i=0; i<2; i++) {
//        diff(i,i) = 1.;
//        diff(2+i,i) = 1.;
//    }
}

/// compute the reference solution and return created mesh
TPZCompMesh *ComputeReferenceSolution(TPZGeoMesh *gmesh, int porder, TPZVec<TPZCompMesh *> &meshvec)
{
    TPZCompMesh *CHDivPressureMesh = CreateReferenceCMesh(gmesh, meshvec, porder);
    //calculo solution
    TPZAnalysis an(CHDivPressureMesh);
    std::cout << "Assembling and Solving " << CHDivPressureMesh->NEquations() << " equations\n";
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(CHDivPressureMesh);
    strmat.SetNumThreads(0);
    an.SetStructuralMatrix(strmat);

#else
    TPZSkylineStructMatrix strmat(CHDivPressureMesh);
#endif
#ifndef PZDEBUG
    //    skyl.SetNumThreads(16);
#endif
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();
    if(0)
    {
        std::ofstream global("Global.nb");
        TPZAutoPointer<TPZStructMatrix> strmat = an.StructMatrix();
        an.Solver().Matrix()->Print("Glob = ",global,EMathematicaInput);
        an.Rhs().Print("Rhs = ",global,EMathematicaInput);
    }
    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
#ifdef PZDEBUG
    {
        std::ofstream out("MeshWithSol.txt");
        CHDivPressureMesh->Print(out);
    }
#endif
    an.LoadSolution(); // compute internal dofs
                       //    an.Solution().Print("sol = ");
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, an.Mesh());
    //    TPZBuildMultiphysicsMesh::TransferFromMeshes(cmeshes, an.Mesh());
    //    for (int i=0; i<cmeshes.size(); i++) {
    //        cmeshes[i]->Solution().Print("sol = ");
    //    }
    //    cmeshes[0]->Solution().Print("solq = ");
    //    cmeshes[1]->Solution().Print("solp = ");
    UnwrapMesh(an.Mesh());
    std::string plotfile("referencesolution.vtk");
    TPZStack<std::string> scalnames,vecnames;
    scalnames.Push("Pressure");
    scalnames.Push("Permeability");
    vecnames.Push("Derivative");
    vecnames.Push("Flux");
    an.DefineGraphMesh(CHDivPressureMesh->Dimension(), scalnames, vecnames, plotfile);
    int resolution = 0;
    an.PostProcess(resolution,CHDivPressureMesh->Dimension());

    return CHDivPressureMesh;
}

/// create the computational mesh of the reference solution
TPZCompMesh *CreateReferenceCMesh(TPZGeoMesh *gmesh, TPZVec<TPZCompMesh *> &meshvec, int porder)
{
    meshvec.resize(2);
    int dimension = 2;
    meshvec[0] = CreateHDivMHMMesh(gmesh,porder);
    meshvec[1] = CreatePressureMHMMesh(gmesh, porder, dimension);
    int hdivplusplus=1;
    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshvec[0], hdivplusplus);
    TPZCompMeshTools::SetPressureOrders(meshvec[0], meshvec[1]);
    TPZCompMesh *cmesh = CreateHDivPressureMHMMesh(meshvec);
    for (int i=0; i<10; i++) {
        TPZMaterial *mat = cmesh->FindMaterial(i);
        if (mat) {
            TPZMixedPoisson *mixed = dynamic_cast<TPZMixedPoisson *>(mat);
            if (!mixed) {
                DebugStop();
            }
            TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Permeability);
            dummy->SetPolynomialOrder(0);
            TPZAutoPointer<TPZFunction<STATE> > func(dummy);
            mixed->SetPermeabilityFunction(func);
        }
    }
    TPZCompMeshTools::GroupElements(cmesh);
    TPZCompMeshTools::CreatedCondensedElements(cmesh, true);
    return cmesh;
}

void UnwrapMesh(TPZCompMesh *cmesh)
{
    long nel = cmesh->NElements();
    bool change = true;
    while(change)
    {
        change = false;
        for (long el=0; el<nel; el++) {
            
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

REAL objectivefunc(REAL K1, REAL K2, REAL K3, REAL K4, REAL Lambda)
{
    if (abs(Lambda-1.) < 1.e-6)
    {
        REAL val = 8*K2*K3*(K1*K1*(K3*(-K3 + K4) + K2*(K3 + K4)) +
                             K2*K4*(K2*(K3 - K4) + K3*(K3 + K4)) +
                             K1*(K2*K2*(K3 + K4) + K3*K4*(K3 + K4) +
                                 K2*(K3*K3 + 6*K3*K4 + K4*K4)) +
                             (K1 + K2)*(K2 + K3)*(K1 + K4)*(K3 + K4)*cos(M_PI*Lambda));
        return val;
    }
    REAL val = (8*K2*K3*(K1*K1*(K3*(-K3 + K4) + K2*(K3 + K4)) +
              K2*K4*(K2*(K3 - K4) + K3*(K3 + K4)) +
              K1*(K2*K2*(K3 + K4) + K3*K4*(K3 + K4) +
                  K2*(K3*K3 + 6*K3*K4 + K4*K4)) +
              (K1 + K2)*(K2 + K3)*(K1 + K4)*(K3 + K4)*cos(M_PI*Lambda))*tan((M_PI*Lambda)/2.))/
    (K4*(K1*(K2*(K3 - K4) - K3*(K3 + K4)) + K2*(K2*(K3 - K4) + K3*(K3 + K4)) +
         (K1 + K2)*(K2 + K3)*(K3 + K4)*cos(M_PI*Lambda)));
    return val;
}

REAL Power(REAL val, int expon)
{
    if (expon != 2) {
        DebugStop();
    }
    return val*val;
}

REAL Sec(REAL val)
{
    return 1./cos(val);
}

REAL Dobjectivefunc(REAL K1, REAL K2, REAL K3, REAL K4, REAL Lambda)
{
    const double Pi = M_PI;
    REAL nom = 4*K2*K3*Pi*(2*Power(K1 + K2,2)*Power(K2 + K3,2)*(K1 + K4)*Power(K3 + K4,2)*cos(Pi*Lambda) +
                           2*(K1 + K2)*(K2 + K3)*(K3 + K4)*
                           (K1*K3*(K1*(K2 - 3*K3) + K2*(K2 + K3)) +
                            (Power(K1,2)*(K2 + K3) + K2*K3*(K2 + K3) +
                             K1*(Power(K2,2) + 10*K2*K3 + Power(K3,2)))*K4 +
                            (K2*(-3*K2 + K3) + K1*(K2 + K3))*Power(K4,2) -
                            2*K1*(K2 + K3)*K4*(K1 + K2 + K3 + K4)*cos(Pi*Lambda)) +
                           4*Power(K1*K3 - K2*K4,2)*(Power(K2,2)*K4 + K1*(Power(K3,2) + (K2 + K3)*K4))*
                           Power(Sec((Pi*Lambda)/2.),2));
//    4*K2*K3*M_PI*(
//    2*(K1 + K2)*(K1+K2)*(K2 + K3)*(K2+K3)*(K1 + K4)*(K3 + K4)*(K3+K4)*cos(M_PI*Lambda) +
//                            2*(K1 + K2)*(K2 + K3)*(K3 + K4)*(
//                                K1*K3*(K1*(K2 - 3*K3) + K2*(K2 + K3)) +
//                                    ((K1*K1)*(K2 + K3) + K2*K3*(K2 + K3) + K1*((K2*K2) + 10*K2*K3 + (K3*K3)))*K4 +
//                             (K2*(-3*K2 + K3) + K1*(K2 + K3))*(K4*K4) -
//                             2*K1*(K2 + K3)*K4*(K1 + K2 + K3 + K4)*cos(M_PI*Lambda)
//                                                             ) +
//                            4*(K1*K3 - K2*K4)*(K1*K3 - K2*K4)*((K2*K2)*K4 + K1*((K3*K3) + (K2 + K3)*K4))/((cos((M_PI*Lambda)/2.)*cos((M_PI*Lambda)/2.)))
//                             );
    REAL denom =
    (K4*
     (K1*(K2 - K3)*K3 - K1*(K2 + K3)*K4 + K2*(K2*(K3 - K4) + K3*(K3 + K4)) + (K1 + K2)*(K2 + K3)*(K3 + K4)*cos(M_PI*Lambda))
     *(K1*(K2 - K3)*K3 - K1*(K2 + K3)*K4 + K2*(K2*(K3 - K4) + K3*(K3 + K4)) + (K1 + K2)*(K2 + K3)*(K3 + K4)*cos(M_PI*Lambda))
     );
    REAL val = nom/denom;
    return val;
}

REAL StartGuess(REAL K1, REAL K2, REAL K3, REAL K4)
{
    REAL inc = 0.1;
    REAL start = 0.60;
    REAL val = objectivefunc(K1,K2,K3,K4,start);
    if(val > 0.) DebugStop();
    while(abs(inc) > 1.e-2)
    {
        REAL val2 = objectivefunc(K1,K2,K3,K4,start+inc);
        if (val2 < 0) {
            start += inc;
            if (start > 1.-1.e-6) {
                start -=inc;
                inc /=2.;
            }
        }
        else
        {
            inc /= 2.;
        }
    }
    return start;
}

REAL ConvergeNewton(REAL K1, REAL K2, REAL K3, REAL K4, REAL start)
{
    REAL val = objectivefunc(K1, K2, K3, K4, start);
    int nstep = 0;
    while(abs(val) > 1.e-6 && nstep < 15)
    {
        REAL deriv = Dobjectivefunc(K1, K2, K3, K4, start);
        start -= val/deriv;
        val = objectivefunc(K1, K2, K3, K4, start);
        nstep++;
    }
    if (nstep == 15) {
        std::cout << "val = " << val << " lambda = " << start << " nsteps " << nstep;
    }
    return start;
}

REAL SolutionRegularity(TPZFMatrix<REAL> &perms)
{
    REAL K1,K2,K3,K4;
    K1 = perms(0,0);
    K2 = perms(1,0);
    K3 = perms(1,1);
    K4 = perms(0,1);
    REAL lambda = StartGuess(K1, K2, K3, K4);
    lambda = ConvergeNewton(K1, K2, K3, K4, lambda);
    return lambda;
}

/// Analise the regularity of the subdomain problems
void AnalyseRegularity(const TPZVec<int> &pos0,const TPZVec<int> &nelx, TPZVec<int> &nsub, TPZFMatrix<REAL> &lowestexp)
{
#ifdef PZDEBUG
    if(nelx[0]%nsub[0] || nelx[1]%nsub[1])
    {
        DebugStop();
    }
#endif
    lowestexp.Redim(nelx[0]-1, nelx[1]-1);
    for (int ir=0; ir<nelx[0]-1; ir++) {
        for (int ic=0; ic<nelx[1]-1; ic++) {
            TPZFNMatrix<4,REAL> perm(2,2);
            perm(0,0) = gPorous(pos0[0]+ir,pos0[1]+ic);
            perm(1,0) = gPorous(pos0[0]+ir+1,pos0[1]+ic);
            perm(0,1) = gPorous(pos0[0]+ir,pos0[1]+ic+1);
            perm(1,1) = gPorous(pos0[0]+ir+1,pos0[1]+ic+1);
            lowestexp(ir,ic) = SolutionRegularity(perm);
            if (lowestexp(ir,ic) < 0 || lowestexp(ir,ic) > 1. || isnan(lowestexp(ir,ic))) {
                DebugStop();
            }
        }
    }
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

/// copy the solution between one computation mesh to the other assuming the geometric elements match
void CopySolution(TPZCompMesh *from, TPZCompMesh *to)
{
    long nelem = from->NElements();
    TPZGeoMesh *gfrom = from->Reference();
    TPZGeoMesh *gto = to->Reference();
    for (long el = 0; el < nelem; el++) {
        TPZCompEl *celfrom = from->Element(el);
        if(!celfrom) continue;
        TPZGeoEl *gelfrom = celfrom->Reference();
        if(!gelfrom) continue;
        if (gelfrom->Dimension() != 1) {
            DebugStop();
        }
        long index = gelfrom->Index();
        TPZGeoEl *gelto = gto->Element(index);
        TPZCompEl *celto = gelto->Reference();
        int nc = celfrom->NConnects();
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &cfrom = celfrom->Connect(ic);
            TPZConnect &cto = celto->Connect(ic);
            long seqfrom = cfrom.SequenceNumber();
            long seqto = cto.SequenceNumber();
            long size = from->Block().Size(seqfrom);
            for (int i=0; i<size; i++) {
                to->Block()(seqto,0,i,0) = from->Block()(seqfrom,0,i,0);
            }
        }
    }
    to->LoadSolution(to->Solution());
}

/// Compute the differences at the submesh level
void ComputeDifferencesBySubmesh(TRunConfig &config, TPZMHMeshControl &MHM, TPZMHMixedMeshControl &MHMixed, const std::string &filename)
{
    long nsubmeshes = MHM.Coarse_to_Submesh().size();
    TPZManVector<STATE> difference(nsubmeshes,0.);
    TPZManVector<STATE,10> square_errors(3,0.);
    std::map<long,long>::iterator it;
    long count = 0;
    for (it = MHM.Coarse_to_Submesh().begin(); it != MHM.Coarse_to_Submesh().end(); it++)
    {
        square_errors[0] = 0;
        square_errors[1] = 0;
        square_errors[2] = 0;
        long skelindex = it->first;
        long MHM_index = it->second;
        if (MHMixed.Coarse_to_Submesh().find(skelindex) == MHMixed.Coarse_to_Submesh().end()) {
            DebugStop();
        }
        long MHMixed_index = MHMixed.Coarse_to_Submesh()[skelindex];
        if (MHM_index == -1 || MHMixed_index == -1) {
            continue;
        }
        TPZSubCompMesh *submhm = dynamic_cast<TPZSubCompMesh *>(MHM.CMesh()->Element(MHM_index));
        TPZSubCompMesh *submhmixed = dynamic_cast<TPZSubCompMesh *>(MHMixed.CMesh()->Element(MHMixed_index));
        TPZCompMeshTools::ComputeDifferenceNorm(submhmixed, submhm, square_errors);
        difference[count] = square_errors[0];
        //        count++;
        //        TPZCompMeshTools::ComputeDifferenceNorm(submhm, submhmixed, square_errors);
        //        difference[count] = square_errors[0];
        count++;
    }
    //    ComputationalMesh->Reference()->ResetReference();
    //    ComputationalMesh->LoadReferences();
    //    return 0;
    {
        std::ofstream out("DiffResults.nb",std::ios::app);
        out << "(* domain size " << config.nelxcoarse << " " << config.nelycoarse << " num subdomains " << count << " *)\n";
        out << "AppendTo[results, {";
        out << " "  << config.numHDivisions << " , " << config.pOrderInternal << " ,";
        out << " {";
        for (long el=0; el<difference.size(); el++) {
            out << difference[el];
            if (el != difference.size()-1) {
                out << " , ";
            }
        }
        out << " } }];\n";
    }
    
}

