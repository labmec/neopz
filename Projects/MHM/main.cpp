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
#include "tpzhierarquicalgrid.h"
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

/// function that randomly refines some elements
void RandomRefine(TPZGeoMesh *gmesh, TRunConfig &config, int nref);

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

TPZFMatrix<REAL> gPorous(500,100,0.);

TAnalyticSolution *example = 0;

int main33(int argc, char *argv[])
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
    Configuration.numHDivisions = 1;
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
        Configuration.nelxcoarse = 4;
        Configuration.nelycoarse = 4;
    }
    Configuration.Hybridize = 1;
    Configuration.Condensed = 1;

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
#ifdef MACOSX
        std::ifstream pores("../porous.txt");
#else
        std::ifstream pores("porous.txt");
#endif
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
    TPZGeoMesh * gmesh;
    bool UseGenGridQ = false;
    REAL Lx = 1000.,Ly = 100., Lz = 10;
    int nref = 1;
    
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
    std::string quad = "QuadByTriangles";
    std::string triangle = "TriangleBy9Triangles";
    TPZAutoPointer<TPZRefPattern> refpatquad = DivideQuadbyTriangles(quad);
    TPZAutoPointer<TPZRefPattern> refpattriangle = DivideTriangleby9Triangles(triangle);
    
    int nelx = 10;
    int nely = 4;
    TPZVec<long> coarseindices;
    gmesh = MalhaGeomFred(nelx, nely, quad, triangle, coarseindices);
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
//        RandomRefine(gmesh, Configuration,1);
        
    }
    
    TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
    TPZAutoPointer<TPZMHMixedHybridMeshControl> MHM;
    TPZAutoPointer<TPZMHMixedMeshControl> MHMixed;
    std::stringstream MHMPref, MHMMixedPref;

    if(1)
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmesh);
        TPZMHMixedHybridMeshControl *mhm = new TPZMHMixedHybridMeshControl(gmeshauto,coarseindices);
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
        
        meshcontrol.CreateSkeletonElements(skeleton);
        
        meshcontrol.DivideSkeletonElements(Configuration.numDivSkeleton);
        if (Configuration.Hybridize)
        {
            meshcontrol.Hybridize(secondskeleton, matpressure);
        }
        
        bool substructure = (bool) Configuration.Condensed;
        meshcontrol.BuildComputationalMesh(substructure);
#ifdef PZDEBUG
        if(0)
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
        if(0)
        {
            std::ofstream gfile("geometry.txt");
            gmesh->Print(gfile);

            std::ofstream out_mhm("MHM_h1.txt");
            meshcontrol.CMesh()->Print(out_mhm);

        }
#endif
    
    TPZCompMesh * CHDivPressureMesh = meshcontrol.CMesh().operator->();

    std::cout << "Number of equations " << CHDivPressureMesh->NEquations() << std::endl;
    
    
//    bool KeepOneLagrangian = true;

    std::cout << "Reduced number of equations " << CHDivPressureMesh->NEquations() << std::endl;
    
    //calculo solution
    TPZAnalysis an(CHDivPressureMesh);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(CHDivPressureMesh);
    strmat.SetNumThreads(16);
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
        std::ofstream out("MHM_hdiv_sol.txt");
        an.Mesh()->Print(out);
    }
#endif
    an.LoadSolution(); // compute internal dofs
//    an.Solution().Print("sol = ");
    
    TPZManVector<TPZCompMesh *,5> cmeshes;
    meshcontrol.GetMeshVec(cmeshes);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmeshes, an.Mesh());
    
//    TPZBuildMultiphysicsMesh::TransferFromMeshes(cmeshes, an.Mesh());
//    for (int i=0; i<cmeshes.size(); i++) {
//        cmeshes[i]->Solution().Print("sol = ");
//    }
//    cmeshes[0]->Solution().Print("solq = ");
//    cmeshes[1]->Solution().Print("solp = ");
    std::string plotfile("mixed_solution.vtk");
    TPZStack<std::string> scalnames,vecnames;
    scalnames.Push("Pressure");
    scalnames.Push("Permeability");
    vecnames.Push("Derivative");
    vecnames.Push("Flux");
    an.DefineGraphMesh(CHDivPressureMesh->Dimension(), scalnames, vecnames, plotfile);
    int resolution = 0;
    an.PostProcess(resolution,CHDivPressureMesh->Dimension());
    TPZManVector<STATE,10> square_errors(3,0.);
    CHDivPressureMesh->Reference()->ResetReference();
    CHDivPressureMesh->LoadReferences();
    TPZCompMeshTools::ComputeDifferenceNorm(ReferenceCMesh, CHDivPressureMesh, square_errors);
    
    for (int i=0; i<square_errors.size(); i++) {
        square_errors[i] = sqrt(square_errors[i]);
    }
    std::cout << "The error norms of the differences are " << square_errors << std::endl;
    return 0;
}


int main_not_used(int argc, char *argv[])
{
    HDivPiola = 1;
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);

    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    if(problemasuave || problemaarctan){
        gmesh= MalhaGeom2(1, 1);}
        std::cout << "Number of equations MHM equals " << MHM->CMesh()->NEquations() << std::endl;
    }
    
    if(1)
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmesh);
        TPZMHMixedHybridMeshControl *mhm = new TPZMHMixedHybridMeshControl(gmeshauto,coarseindices);
        MHMMixedPref << "MHMixed";
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
        
        bool substructure = (bool) Configuration.Condensed;
        meshcontrol.BuildComputationalMesh(substructure);
#ifdef PZDEBUG
        if(0)
        {
            std::ofstream file("GMeshControlHDiv.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file);
        }
#endif
#ifdef PZDEBUG
        if(0)
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
    if(Configuration.LagrangeMult)
    {
        MHMPref << "_Lagr";
        MHMMixedPref << "_Lagr";
    }
    if (Configuration.Hybridize) {
        MHMPref << "_Hybr";
        MHMMixedPref << "_Hybr";
    }
    // compute the MHM solution
    Configuration.fGlobalSystemWithLocalCondensationSize = MHM->fGlobalSystemWithLocalCondensationSize;
    Configuration.fGlobalSystemSize = MHM->fGlobalSystemSize;
    Configuration.fNumeq = MHM->fNumeq;
    SolveProblem(MHM->CMesh(), MHM->GetMeshes(), example, MHMPref.str(), Configuration);
    
    // compute the MHM H(div) solution
    Configuration.fGlobalSystemWithLocalCondensationSize = MHMixed->fGlobalSystemWithLocalCondensationSize;
    Configuration.fGlobalSystemSize = MHMixed->fGlobalSystemSize;
    Configuration.fNumeq = MHMixed->fNumeq;
    SolveProblem(MHMixed->CMesh(), MHMixed->GetMeshes(), example, MHMMixedPref.str(), Configuration);
    
//    CopySolution(MHMixed->CMesh().operator->(), MHM->CMesh().operator->());
    
    if (Configuration.Condensed)
    {
        std::string filename = "MHMixed_" + configuration + ".txt";
        std::ofstream out(filename);
        PrintElements(MHMixed->CMesh().operator->(), out);
    }
    if(Configuration.Condensed)
    {
        std::string filename = "MHM_" + configuration + ".txt";
        std::ofstream out(filename);
        PrintElements(MHM->CMesh().operator->(), out);
    }
    
//    ComputeDifferencesBySubmesh(Configuration, MHM, MHMixed, "DiffResults.nb");
    if(!example)
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
<<<<<<< HEAD
#endif
    
    //calculo solution
    TPZAnalysis an(mhm.CMesh());
    TPZSkylineStructMatrix skyl(mhm.CMesh());
    an.SetStructuralMatrix(skyl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();
    an.Solve();
    an.Solution().Print("solution = ");
    long neq = an.Solution().Rows();
    long numeq = MIN(10, neq);
    TPZManVector<long> equationindices(numeq);
    for (int i=0; i<numeq; i++) {
        equationindices[i] = i;
    }
    an.ShowShape("Shapes.vtk", equationindices);
    an.SetStep(0);
    
    if(problemasuave || problemaarctan){
        std::string plotfile("result.vtk");
        TPZStack<std::string> scalnames,vecnames;
        scalnames.Push("Solution");
        scalnames.Push("ExactSolution");
        scalnames.Push("PressureConstante");
        an.DefineGraphMesh(mhm.CMesh()->Dimension(), scalnames, vecnames, plotfile);
        an.PostProcess(0,2);
    }
    
//    an.SetExact(*SolExataSteklov);
//    TPZVec<REAL> erros(3);
//    an.PostProcessError(erros);
    
    //    //construir elementos 1D de interface
    //    TPZCompMesh * cmesh1 = MalhaCompTemporaria(gmesh);
    //    gmesh->ResetReference();
    //
    //    //mudar matId dos elementos 1D de interface
    //    ChangeIndex(gmesh, matCoarse);
    //    ofstream arg3("gmesh3.txt");
    //	gmesh->Print(arg3);
    //
    //    ofstream file1("malhageometricaCoarse.vtk");
    //    TPZVTKGeoMesh::PrintGMeshVTK(gmesh.operator->(), file1, true);
    //
    //
    //    if(coarseindex.find(6) != coarseindex.end())
    //    {
    //        std::cout << "\n\n\nNAO ACHEI O NUMERO DESEJADO\n\n\n";
    //    }
    //
    //
    //    std::set<long>::iterator it;
    //    std::cout << "coarse index: \n";
    //    for (it=coarseindex.begin(); it!=coarseindex.end(); ++it)
    //        std::cout << ' ' << *it;
    //    std::cout << '\n';
    //
    ////-------- malha mais fina -------
    //    //refinamento uniforme dos elementos 2D
    //    dims.Resize(1, 0);
    //    dims[0]=2;
    //    RefinamentoUniforme(gmesh, 1, dims);
    //    ofstream arg4("gmesh4.txt");
    //	gmesh->Print(arg4);
    //
    //    ofstream file2("malhageometricaFina.vtk");
    //    TPZVTKGeoMesh::PrintGMeshVTK(gmesh.operator->(), file2, true);
    //
    ////malha computacional
    //    TPZCompMesh * cmesh = MalhaComp2(gmesh,1,coarseindex);
    //    ofstream arg5("cmesh.txt");
    //	cmesh->Print(arg5);
    //
    //    InterfaceToCoarse(cmesh, matInterno, matCoarse, matInterface);
    //
    //    ofstream arg6("cmesh2.txt");
    //	cmesh->Print(arg6);
    
    return EXIT_SUCCESS;
}

TPZGeoMesh *MalhaGeom2(REAL Lx, REAL Ly)
{
    int Qnodes = 4;
	long dim = 2;
    
	TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <long> TopolQuad(4);
	TPZVec <long> TopolLine(2);
	
	//indice dos nos
	long id = 0;
	REAL valx;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = Lx - xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,Ly);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
    
	//indice dos elementos
	id = 0;
    
    //elementos internos
    TopolQuad[0] = 0;
	TopolQuad[1] = 1;
	TopolQuad[2] = 2;
	TopolQuad[3] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matInterno,*gmesh);
	id++;
    
    //elementos de contorno
	TopolLine[0] = 0;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
	id++;
	
	TopolLine[0] = 1;
	TopolLine[1] = 2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
	id++;
	
	TopolLine[0] = 2;
	TopolLine[1] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
	id++;
	
	TopolLine[0] = 3;
	TopolLine[1] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
	id++;
    
    //construir a malha
	gmesh->BuildConnectivity();
	
	return gmesh;
}


TPZGeoMesh *GMeshSteklov(bool triang_elements)
{
    TPZManVector<int,2> nx(2,2);
    REAL extent = 1;
    nx[1] =1;
    TPZManVector<REAL,3> x0(3,0.),x1(3,extent);
    x0[0] = -extent;
    TPZGenGrid gengrid(nx,x0,x1);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    if(triang_elements)
    {
        gengrid.SetElementType(ETriangle);
    }
    gengrid.Read(gmesh);
    
    //elementos de contorno
    TPZManVector<REAL,3> firstpoint(3,0.),secondpoint(3,0.);
    firstpoint[0] = extent;
    secondpoint[0] = extent;
    secondpoint[1] = extent;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc2);
    gengrid.SetBC(gmesh,6,bc3);
    gengrid.SetBC(gmesh,7,bc4);
    firstpoint[0] = -extent;
    secondpoint[0] = 0.;
    secondpoint[1] = 0.;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc5);
    firstpoint = secondpoint;
    secondpoint[0] = extent;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc1);
    
    return gmesh;
}

TPZCompMesh* MalhaCompTemporaria(TPZAutoPointer<TPZGeoMesh>  gmesh)
{
	/// criar materiais
	int dim = 2;
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matInterno,dim);
	TPZMaterial * mat1(material);
	material->NStateVariables();
    
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->InsertMaterialObject(mat1);
	
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	
    TPZMaterial * BCondD1 = material->CreateBC(mat1, bc2,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCondD1);
	
	TPZMaterial * BCondD2 = material->CreateBC(mat1, bc4,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCondD2);
    
	TPZMaterial * BCondN1 = material->CreateBC(mat1, bc1,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCondN1);
    
    TPZMaterial * BCondN2 = material->CreateBC(mat1, bc3,dirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCondN2);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    cmesh->AutoBuild();
    
    TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
    
    return cmesh;
}

void InsertMaterialObjects(TPZMHMeshControl &control)
{
    TPZCompMesh &cmesh = control.CMesh();
	/// criar materiais
	int dim = cmesh.Dimension();
    TPZMatLaplacianLagrange *material1 = new TPZMatLaplacianLagrange(matInterno,dim);
    
    material1->SetParameters(10., 0.);
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Porosity);
    dummy->SetPolynomialOrder(0);
    TPZAutoPointer<TPZFunction<STATE> > func(dummy);
    material1->SetPermeabilityFunction(func);

    
	TPZMaterial * mat1(material1);
    
    TPZMat1dLin *materialCoarse = new TPZMat1dLin(matCoarse);
    TPZFNMatrix<1,STATE> xk(1,1,0.),xb(1,1,0.),xc(1,1,0.),xf(1,1,0.);
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
    TPZMaterial * BCondD1 = material1->CreateBC(mat1, bc1,dirichlet, val1, val2);
    //TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet1 = new TPZDummyFunction<REAL>(DirichletValidacao);
    //BCondD1->SetForcingFunction(bcmatDirichlet1);
    cmesh.InsertMaterialObject(BCondD1);
    
    //BC -2
	TPZMaterial * BCondD2 = material1->CreateBC(mat1, bc2,neumann, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet2 = new TPZDummyFunction<REAL>(DirichletValidacaoMenos);
    BCondD2->SetForcingFunction(bcmatDirichlet2);
    cmesh.InsertMaterialObject(BCondD2);
    
    //BC -3
	TPZMaterial * BCondD3 = material1->CreateBC(mat1, bc3,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet3 = new TPZDummyFunction<REAL>(DirichletValidacao);
    BCondD3->SetForcingFunction(bcmatDirichlet3);
    cmesh.InsertMaterialObject(BCondD3);
    
    //BC -4
	TPZMaterial * BCondD4 = material1->CreateBC(mat1, bc4,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet4 = new TPZDummyFunction<REAL>(DirichletValidacao);
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
    TPZMixedPoisson * mat = new TPZMixedPoisson(1,dim);
    mat->SetSymmetric();
    mat->SetPermeability(10.);
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Porosity);
    dummy->SetPolynomialOrder(0);
    TPZAutoPointer<TPZFunction<STATE> > func(dummy);
    mat->SetPermeabilityFunction(func);
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat);
    
    
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, -1, typeFlux, val1, val2Flux);
    TPZAutoPointer<TPZFunction<STATE> > force = new TPZDummyFunction<STATE>(DirichletValidacao);
    //    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, -2, typePressure, val1, val2Pressure);
    bcS->SetForcingFunction(0, force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    
    
    

}

void InsertMaterialObjectsSuave(TPZCompMesh &cmesh)
{
	/// criar materiais
	int dim = cmesh.Dimension();
    TPZMatLaplacianLagrange *materialFiner = new TPZMatLaplacianLagrange(matInterno,dim);
    
    TPZAutoPointer<TPZFunction<REAL> > forcef;
    TPZAutoPointer<TPZFunction<STATE> > solExata;
    
    if(problemaarctan)
    {
        forcef = new TPZDummyFunction<REAL>(ForcingTang);
        solExata = new TPZDummyFunction<STATE>(SolArcTan);
    }else
    {
        forcef = new TPZDummyFunction<REAL>(ForceSuave);
        solExata = new TPZDummyFunction<STATE>(SolSuave);
    }
    
    materialFiner->SetForcingFunction(forcef);
    materialFiner->SetForcingFunctionExact(solExata);
    
	TPZMaterial * mat1(materialFiner);
    
    TPZMat1dLin *materialSkeleton = new TPZMat1dLin(matCoarse);
    TPZFNMatrix<1,STATE> xk(1,1,0.),xb(1,1,0.),xc(1,1,0.),xf(1,1,0.);
    materialSkeleton->SetMaterial(xk, xc, xb, xf);
    
	cmesh.InsertMaterialObject(mat1);
    cmesh.InsertMaterialObject(materialSkeleton);
	
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,1.), val2(2,1,0.);
    TPZMaterial * BCondD1;
    TPZMaterial * BCondD2;
    TPZMaterial * BCondD3;
    TPZMaterial * BCondD4;
    
    if(problemaarctan)
    {
        //u=0 no contorno de Omega
        BCondD1 = materialFiner->CreateBC(mat1, bc1,neumann, val1, val2);
        BCondD2 = materialFiner->CreateBC(mat1, bc2,neumann, val1, val2);
        BCondD3 = materialFiner->CreateBC(mat1, bc3,neumann, val1, val2);
        BCondD4 = materialFiner->CreateBC(mat1, bc4,neumann, val1, val2);
    }
    else
    {
        //BC -1
        BCondD1 = materialFiner->CreateBC(mat1, bc1,neumann, val1, val2);
        TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet1 = new TPZDummyFunction<REAL>(DirichletSuave);
        BCondD1->SetForcingFunction(bcmatDirichlet1);
        
        //BC -2
        BCondD2 = materialFiner->CreateBC(mat1, bc2,neumann, val1, val2);
        TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet2 = new TPZDummyFunction<REAL>(DirichletSuave);
        BCondD2->SetForcingFunction(bcmatDirichlet2);
        
        //BC -3
        BCondD3 = materialFiner->CreateBC(mat1, bc3,neumann, val1, val2);
        TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet3 = new TPZDummyFunction<REAL>(DirichletSuave);
        BCondD3->SetForcingFunction(bcmatDirichlet3);
        
        //BC -4
        BCondD4 = materialFiner->CreateBC(mat1, bc4,neumann, val1, val2);
        TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet4 = new TPZDummyFunction<REAL>(DirichletSuave);
        BCondD4->SetForcingFunction(bcmatDirichlet4);
    }
    
    cmesh.InsertMaterialObject(BCondD1);
    cmesh.InsertMaterialObject(BCondD2);
    cmesh.InsertMaterialObject(BCondD3);
    cmesh.InsertMaterialObject(BCondD4);
}


TPZCompMesh* MalhaComp2(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,std::set<long> coarseindex)
{
	/// criar materiais
	int dim = 2;
    TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
    
    DebugStop();
//    InsertMaterialObjects(*cmesh);
    TPZMatPoisson3d *material2 = new TPZMatPoisson3d(matCoarse,dim);
    cmesh->InsertMaterialObject(material2);
    
    //    REAL diff = -1.;
    //	REAL conv = 0.;
    //	TPZVec<REAL> convdir(3,0.);
    //	REAL flux = 8.;
    //
    //	material1->SetParameters(diff, conv, convdir);
    //	material1->SetInternalFlux( flux);
    
    cmesh->SetDefaultOrder(pOrder);
    
    cmesh->SetAllCreateFunctionsContinuous();
    
    //Criar elementos computacionais malha MHM
    int nel = gmesh->NElements();
    int matid, eldim;
    long index;
    int hassubel, nsubels;
    int iel, is;
    
    TPZGeoEl *gel = NULL;
    TPZGeoEl *gsubel = NULL;
    
    for(iel = 0; iel<nel; iel++)
    {
        gel = gmesh->ElementVec()[iel];
        if(!gel) DebugStop();
        
        eldim = gel->Dimension();
        matid = gel->MaterialId();
        
        //elementos de dimensao = dim (malha fina)
        if(eldim==dim)
        {
            index = gel->Index();
            if(coarseindex.find(index) != coarseindex.end())
            {
                nsubels = gel->NSubElements();
                for (is=0; is<nsubels; is++)
                {
                    gsubel = gel->SubElement(is);
                    hassubel = gsubel->HasSubElement();
                    if(!hassubel){
                        cmesh->CreateCompEl(gsubel,index);
                    }
                }
                for (is=0; is<nsubels; is++)
                {
                    gsubel = gel->SubElement(is);
                    hassubel = gsubel->HasSubElement();
                    if(!hassubel){
                        gsubel->ResetReference();
                    }
                }
            }
            continue;
        }
        
        //elementos de dimensao = dim-1
        
        //malha coarse
        if(matid==matCoarse)
        {
            cmesh->CreateCompEl(gel, index);
            gel->ResetReference();
            
            continue;
        }
        
        //elementos de contorno
        hassubel =gel->HasSubElement();
        if(!hassubel)
        {
            cmesh->CreateCompEl(gel, index);
            gel->ResetReference();
        }
    }
    
    cmesh->LoadReferences();
    cmesh->ExpandSolution();
    
    return cmesh;
}



void RefinamentoUniforme(TPZAutoPointer<TPZGeoMesh> gmesh, int nref,TPZVec<int> dims)
{
    RefinamentoUniforme(gmesh.operator->(), nref, dims);
}



void RefinamentoUniforme(TPZGeoMesh *gmesh, int nref,TPZVec<int> dims)
{
    
    int ir, iel, k;
    int nel=0, dim=0;
    int ndims = dims.size();
	for(ir = 0; ir < nref; ir++ )
    {
		TPZVec<TPZGeoEl *> filhos;
        nel = gmesh->NElements();
        
		for (iel = 0; iel < nel; iel++ )
        {
			TPZGeoEl * gel = gmesh->ElementVec()[iel];
            if(!gel) DebugStop();
            
            dim = gel->Dimension();
            
            for(k = 0; k<ndims; k++)
            {
                if(dim == dims[k])
                {
                    gel->Divide (filhos);
                    break;
                }
            }
		}
	}
    
}

void RefinamentoAdaptado(TPZAutoPointer<TPZGeoMesh> gmesh, TPZStack<TPZManVector<REAL,3> > coordcentro)
{
    int size = coordcentro.NElements();
    std::set<TPZGeoEl *> setgeo;
    
    TPZVec<REAL> qsi(3,0.);
    long iniEl = 0;
    
    TPZStack<TPZGeoEl *> vecgel;
    vecgel.Resize(0);
    int eldim = 2;
    for(int ix = 0; ix < size; ix++)
    {
        TPZGeoEl * gel = NULL;
        
        if(coordcentro[ix][0]== 0. || coordcentro[ix][1]==0. || coordcentro[ix][0]== 1. || coordcentro[ix][1]==1.)
        {
            eldim = 1;
            gel = gmesh->FindElement(coordcentro[ix], qsi, iniEl, eldim);
            eldim = 2;
        }else gel = gmesh->FindElement(coordcentro[ix], qsi, iniEl, eldim);
        if(!gel) DebugStop();
        setgeo.insert(gel);
    }
    
    
    int nel = gmesh->NElements();
    TPZVec<TPZGeoEl *> filhos;
    for(int iel = 0; iel<nel; iel++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        std::set<TPZGeoEl *>::iterator found = setgeo.find(gel);
        
        if(gel==(*found)){
            continue;
        }
        gel->Divide (filhos);
    }
}

void GetElIndexCoarseMesh(TPZGeoMesh *  gmesh, std::set<long> &coarseindex)
{
    int nel = gmesh->NElements();
    int iel;
    int hassubel=0;
    int dim = gmesh->Dimension();
    int eldim;
    for(iel = 0; iel<nel; iel++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        if(!gel) DebugStop();
        
        hassubel = gel->HasSubElement();
        eldim = gel->Dimension();
        if(!hassubel && eldim ==dim)
        {
            coarseindex.insert(gel->Index());
        }
    }
    
}

void ChangeIndex(TPZAutoPointer<TPZGeoMesh> gmesh, int matcoarse1D)
{
    int nel = gmesh->NElements();
    int iel;
    int hassubel, ninterf;
    
    for(iel = 0; iel<nel; iel++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        if(!gel) DebugStop();
        
        hassubel = gel->HasSubElement();
        if(!hassubel)
        {
            ninterf = gel->NumInterfaces();
            if(ninterf > 1) DebugStop();
            if (ninterf==1)
            {
                gel->SetMaterialId(matcoarse1D);
                gel->DecrementNumInterfaces();
            }
        }
    }
    
}

//TPZCompMesh *SkeletonCoarseCompMesh (TPZCompMesh *cmesh)
//{
//
//    //ponteiro para a malha geometrica de mesh
//    TPZGeoMesh *gmesh = cmesh->Reference();
//    if(!gmesh)
//    {
//        DebugStop();
//    }
//
//    //Reseta as referencias do ponteiro para a malha geometrica criada
//    //e criar uma nova malha computacional baseada nesta malha geometrica
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    //TPZCompMesh *newcmesh = new TPZCompMesh(gmesh);
//
//    int dim = cmesh->Dimension();
//    //newcmesh->SetDimModel(dim);
//
//    TPZStack<TPZCompElSide> neighequal,neighsmaller;
//    TPZCompElSide neighbigger;
//    int nneighs=0;
//
//    int nel = cmesh->ElementVec().NElements();
//    for(long el = 0; el < nel; el++)
//	{
//		TPZCompEl * compEl = cmesh->ElementVec()[el];
//		if(!compEl) continue;
//
//        int eldim = compEl->Reference()->Dimension();
//        int elmatId = compEl->Reference()->MaterialId();
//
//        //convencao PZ: elemento de contorno tem sempre Id negativo
//		if(elmatId < 0)
//		{
//			compEl->Reference()->ResetReference();
//		}
//        else if(eldim == dim)
//        {
//            compEl->Reference()->ResetReference();
//
//            int nsides = compEl->Reference()->NSides();
//            int ncorn = compEl->Reference()->NCornerNodes();
//            for(int side = ncorn; side < nsides; side++)
//            {
//                neighequal.Resize(0);
//                neighsmaller.Resize(0);
//
//                TPZCompElSide celside(compEl,side);
//                celside.EqualLevelElementList(neighequal, 0, 0);
//                neighbigger = celside.LowerLevelElementList(0);
//
//                if(neighbigger){
//                    neighequal.Push(neighbigger);
//                }
//                nneighs = neighequal.NElements();
//
//                //                celside.HigherLevelElementList(neighsmaller, 1, 1);
//                //                int nneighsmaller = neighsmaller.size();
//                //                if(nneighs && nneighsmaller)
//                //                {
//                //                    DebugStop();
//                //                }
//
//                if(nneighs != 0)
//                {
//                    //Loop on neighboring elements greater or equal level
//                    for(int i =0; i<nneighs; i++)
//                    {
//                        TPZGeoEl *geoside = neighequal[i].Reference().Element();
//
//                        //Do not assume neighbors by node
//                        if(neighequal[i].Side() < geoside->NCornerNodes()) continue;
//
//                        //verificando se eh elemento 1d
//                        if(geoside->Dimension() == dim-1 && geoside->MaterialId() > 0) continue;
//
//                        long index;
//                        TPZInterpolatedElement *newcompel;
//                        newcompel = dynamic_cast<TPZInterpolatedElement *> (cmesh->CreateCompEl(geoside, index));
//                        newcompel->LoadElementReference();
//                    }
//                }
//            }
//        }
//        else continue;
//	}
//
//    return cmesh;
//}

TPZCompMesh *SkeletonCoarseCompMesh (TPZCompMesh *cmesh, int matId)
{
    
    //ponteiro para a malha geometrica de mesh
    TPZGeoMesh *gmesh = cmesh->Reference();
    if(!gmesh)
    {
        DebugStop();
    }
    
    //Resetar as referencias do ponteiro para a malha geometrica criada
    //e criar uma nova malha computacional baseada nesta malha geometrica
    //gmesh->ResetReference();
    
    //    TPZCompMesh *newcmesh = new TPZCompMesh(gmesh);
    //    gmesh->ResetReference();
    //    newcmesh->LoadReferences();
    
    int dim = cmesh->Dimension();
    //newcmesh->SetDimModel(dim);
    
    TPZStack<TPZCompElSide> neighequal,neighsmaller;
    TPZCompElSide neighbigger;
    int nneighs=0;
    
    int nel = gmesh->ElementVec().NElements();
    for(long el = 0; el < nel; el++)
	{
		TPZGeoEl * gel = gmesh->ElementVec()[el];
		if(!gel) continue;
        if(gel->HasSubElement()) continue;
        
        int eldim = gel->Dimension();
        //int elmatId = gel->MaterialId();
        
        //        //convencao PZ: elemento de contorno tem sempre Id negativo
        //		if(elmatId < 0)
        //		{
        //			geo->Reference()->ResetReference();
        //		}
        if(eldim == dim)
        {
            //compEl->Reference()->ResetReference();
            
            int nsides = gel->NSides();
            int ncorn = gel->NCornerNodes();
            for(int side = ncorn; side < nsides; side++)
            {
                neighequal.Resize(0);
                neighsmaller.Resize(0);
                
                TPZGeoElSide gelside(gel,side);
                gelside.EqualLevelCompElementList(neighequal, 0, 0);
                neighbigger = gelside.LowerLevelCompElementList2(0);
                
                if(neighbigger){
                    neighequal.Push(neighbigger);
                }
                nneighs = neighequal.NElements();
                
                //                gelside.HigherLevelCompElementList2(neighsmaller, 1, 1);
                //                int nneighsmaller = neighsmaller.size();
                //                if(nneighs && nneighsmaller)
                //                {
                //                    DebugStop();
                //                }
                
                if(nneighs != 0)
                {
                    //Loop on neighboring elements greater or equal level
                    for(int i =0; i<nneighs; i++)
                    {
                        TPZGeoEl *geoside = neighequal[i].Reference().Element();
                        
                        //Do not assume neighbors by node
                        if(neighequal[i].Side() < geoside->NCornerNodes()) continue;
                        
                        //verificando se eh elemento 1d
                        if(geoside->Dimension() == dim-1 && geoside->MaterialId() > 0) continue;
                        
                        long index;
                        TPZInterpolatedElement *newcompel;
                        geoside->ResetReference();
                        newcompel = dynamic_cast<TPZInterpolatedElement *> (cmesh->CreateCompEl(geoside, index));
                        newcompel->Reference()->SetMaterialId(matId);
                        newcompel->LoadElementReference();
                        newcompel->Print();
                        cout<<"\n";
                        newcompel->Reference()->Print();
                    }
                }
            }
            gel->ResetReference();
        }
        else continue;
	}
    //cmesh->LoadReferences();
    cmesh->InitializeBlock();
    return cmesh;
}

void InterfaceToCoarse(TPZCompMesh *cmesh, int matvolume, int matskeleton, int matinterface)
{
    int nelem = cmesh->NElements();
    for (int el=0; el<nelem; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        int matid = cel->Reference()->MaterialId();
        if (matid != matskeleton) {
            continue;
        }
        // the element cel is a skeleton element
        TPZGeoEl *gel = cel->Reference();
        std::set<long> celindices;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZCompElSide celskeleton = gelside.Reference();
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            TPZStack<TPZCompElSide> celstack;
            gelside.HigherLevelCompElementList2(celstack, 1, 0);
            long nstack = celstack.NElements();
            for (long i=0; i<nstack; i++) {
                TPZCompElSide celside = celstack[i];
                long celindex = celside.Element()->Index();
                if (celindices.find(celindex) != celindices.end()) {
                    continue;
                }
                celindices.insert(celindex);
                int side = gel->NSides()-1;
                long index;
                TPZGeoEl *gelnew = gel->CreateBCGeoEl(side, matinterface);
                new TPZInterfaceElement(*cmesh, gelnew , index, celside, celskeleton);
                
            }
            neighbour = neighbour.Neighbour();
        }
    }
}


void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    const REAL n = 0;
    
    const REAL x = loc[0];
    const REAL y = loc[1];
    const REAL r = sqrt(x*x+y*y);
    const REAL t = atan2(y,x);
    const REAL sol = pow((REAL)2,0.25 + n/2.)*pow(r,0.5 + n)*cos((0.5 + n)*t);
    u[0] = sol;
    
    du(0,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(x*cos((0.5 + n)*atan2(y,x)) + y*sin((0.5 + n)*atan2(y,x)));
    du(1,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(y*cos((0.5 + n)*atan2(y,x)) - x*sin((0.5 + n)*atan2(y,x)));
    
}

void RefinamentoSingular(TPZAutoPointer<TPZGeoMesh> gmesh,int nref)
{
    RefinamentoSingular(gmesh.operator->(), nref);
}

void RefinamentoSingular(TPZGeoMesh *gmesh,int nref)
{
    long nnodes = gmesh->NNodes();
    long in;
    for (in=0; in<nnodes; in++) {
        TPZGeoNode *gno = &gmesh->NodeVec()[in];
        if (abs(gno->Coord(0))< 1.e-6 && abs(gno->Coord(1)) < 1.e-6) {
            break;
        }
    }
    if (in == nnodes) {
        DebugStop();
    }
    TPZGeoElSide gelside;
    long nelem = gmesh->NElements();
    for (long el = 0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        int ncorner = gel->NCornerNodes();
        for (int ic=0; ic<ncorner; ic++) {
            long nodeindex = gel->NodeIndex(ic);
            if (nodeindex == in) {
                gelside = TPZGeoElSide(gel, ic);
                break;
            }
        }
        if (gelside.Element()) {
            break;
        }
    }
    if (!gelside.Element()) {
        DebugStop();
    }
    for (int iref = 0; iref <nref; iref++) {
        TPZStack<TPZGeoElSide> gelstack;
        gelstack.Push(gelside);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            gelstack.Push(neighbour);
            neighbour = neighbour.Neighbour();
        }
        long nstack = gelstack.size();
        for (long ist=0; ist < nstack; ist++) {
            if (!gelstack[ist].Element()->HasSubElement()) {
                TPZVec<TPZGeoEl *> subel;
                gelstack[ist].Element()->Divide(subel);
            }
        }
    }
}

void SolSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    const REAL x = loc[0];
    const REAL y = loc[1];
    const REAL sol = sin(M_PI*x)*sin(M_PI*y);
    u[0] = sol;
    du.Resize(2, 1);
    du(0,0) = M_PI*cos(M_PI*x)*sin(M_PI*y);
    du(1,0) = M_PI*cos(M_PI*y)*sin(M_PI*x);
}

void ForceSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &force){
    const REAL x = loc[0];
    const REAL y = loc[1];
    const REAL f = 2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    force[0] = f;
}

void DirichletSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    TPZFMatrix<STATE> du(2,1);
    SolSuave(loc,result,du);
}

#define Power pow
#define ArcTan atan
#define Sqrt sqrt

void SolArcTan(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux){
    REAL x = pt[0];
    REAL y = pt[1];
    
    p[0]=0;
    flux(0,0)=0;
    flux(1,0)=0;

    
    
    
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

    long nel = gmesh->NElements();
    int dim = gmesh->Dimension();
    for (long el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->Level() != levelinterface || gel->Dimension() != dimension) {
            continue;
        }
        int nsides = gel->NSides();
        for (int is = gel->NCornerNodes(); is<nsides; is++) {
            if (gel->SideDimension(is) != dimension-1) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->Dimension() == dimension-1) {
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour != gelside) {
                continue;
            }
            TPZGeoElSide neighfather = gelside;
            for (int il = level; il< levelinterface; il++) {
                neighfather = neighfather.Father2();
            }
            if (!neighfather || neighfather.Dimension() != dimension-1) {
                continue;
            }
            
            if (neighbour == gelside) {
                TPZGeoElBC(gelside, 2);
            }
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
    
    control.InsertPeriferalMaterialObjects();
    
}


TPZCompMesh * CreateHDivMHMMesh(TPZGeoMesh * gmesh, int porder)
{

    int dim = gmesh->Dimension();
    TPZCompMesh * cmeshHDiv = new TPZCompMesh(gmesh);
    cmeshHDiv->SetDimModel(dim);
    cmeshHDiv->SetAllCreateFunctionsHDiv();

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
    
    cmeshHDiv->SetDefaultOrder(porder);
    cmeshHDiv->AutoBuild();
    
//    cmeshHDiv->SetDefaultOrder(porder);
//    std::set<int> mat_skeleton;
//    mat_skeleton.insert(2);
//    cmeshHDiv->AutoBuild(mat_skeleton);
//    
//    cmeshHDiv->SetDefaultOrder(porder);
//    std::set<int> mat_vol;
//    mat_vol.insert(1);
//    mat_vol.insert(-1);
//    mat_vol.insert(-2);
//    cmeshHDiv->AutoBuild(mat_vol);
    

#ifdef PZDEBUG
    {
        std::ofstream outmesh("CmeshFlux.txt");
        cmeshHDiv->Print(outmesh);
    }
#endif
    return cmeshHDiv;
}


void DuplicateNeighbouringConnects(TPZCompMesh * HDivMesh)
{
    TPZGeoMesh *gmesh = HDivMesh->Reference();
    int dimension = gmesh->Dimension();
    gmesh->ResetReference();
    HDivMesh->LoadReferences();
    HDivMesh->ComputeNodElCon();
    long nel = HDivMesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = HDivMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();

        if (!gel || gel->Dimension() != dimension) {
            continue;
        }
        int nc = cel->NConnects();
        for (int ic =0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if (c.HasDependency() && c.NElConnected() == 2)
            {
                // duplicate the connect
                long cindex = HDivMesh->AllocateNewConnect(c);
                TPZConnect &newc = HDivMesh->ConnectVec()[cindex];
                newc = c;
                c.DecrementElConnected();
                newc.DecrementElConnected();
                cel->SetConnectIndex(ic, cindex);
            }
        }
    }
    HDivMesh->ExpandSolution();
}

TPZCompMesh * CreatePressureMHMMesh(TPZGeoMesh * gmesh, int porder, int dimension)
{
    int dim = gmesh->Dimension();
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
    
#ifdef PZDEBUG
    {
        std::ofstream outmesh("CmeshPressure.txt");
        cmeshPressure->Print(outmesh);
    }
#endif
    
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
    TPZFMatrix<STATE> val1(1,1,0.), val2Flux(1,1,0.), val2Pressure(1,1,0.);
    
    // Malha computacional
    TPZCompMesh * MixedFluxPressureCmesh = new TPZCompMesh(gmesh);
    
    // Material medio poroso
    TPZMixedPoisson * mat = new TPZMixedPoisson(1,dim);
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
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, -2, typeFlux, val1, val2Flux);
//    bcN->SetBCForcingFunction(0, force);
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
    
#ifdef PZDEBUG
    {
        std::ofstream outmesh("CmeshMixed.txt");
        MixedFluxPressureCmesh->Print(outmesh);
    }
#endif
    
    return MixedFluxPressureCmesh;

}

<<<<<<< HEAD
void HideTheElements(TPZCompMesh * Multiphysics, bool KeepOneLagrangian, TPZVec<long> &coarseindices)
{
    typedef std::set<long> TCompIndexes;
    std::map<long, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = Multiphysics->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    Multiphysics->LoadReferences();
    long nelg = coarseindices.NElements();
    for (long iel=0; iel<nelg; iel++) {
        long el = coarseindices[iel];
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->Dimension() != dim && gel->MaterialId() > 0) {
            DebugStop();
        }
        // we took any neighbour of gel and identified a mapindex with it??
        TPZStack<TPZCompElSide> highlevel;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        gelside.HigherLevelCompElementList3(highlevel, 0, 0);
        long nelst = highlevel.size();
        for (long elst=0; elst<nelst; elst++) {
            ElementGroups[el].insert(highlevel[elst].Element()->Index());
        }
        if (gel->Reference()) {
            if (nelst) {
                DebugStop();
            }
            ElementGroups[el].insert(gel->Reference()->Index());
        }
    }
    std::cout << "Number of element groups " << ElementGroups.size() << std::endl;
    std::map<long,TCompIndexes>::iterator it;
    for (it=ElementGroups.begin(); it != ElementGroups.end(); it++) {
        std::cout << "Group " << it->first << " group size " << it->second.size() << std::endl;
        std::cout << " elements ";
        std::set<long>::iterator its;
        for (its = it->second.begin(); its != it->second.end(); its++) {
            std::cout << *its << " ";
        }
        std::cout << std::endl;
    }
    
    std::set<long> submeshindices;
    TPZCompMeshTools::PutinSubmeshes(Multiphysics, ElementGroups, submeshindices, KeepOneLagrangian);
    std::cout << "After putting in substructures\n";
    
    Multiphysics->ComputeNodElCon();
    Multiphysics->CleanUpUnconnectedNodes();
    for (std::set<long>::iterator it=submeshindices.begin(); it != submeshindices.end(); it++) {
        TPZCompEl *cel = Multiphysics->Element(*it);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            DebugStop();
        }
        TPZCompMeshTools::GroupElements(subcmesh);
        subcmesh->ComputeNodElCon();
        TPZCompMeshTools::CreatedCondensedElements(subcmesh, KeepOneLagrangian);
        subcmesh->SetAnalysisSkyline(16, 0, 0);
    }
    Multiphysics->ComputeNodElCon();
    Multiphysics->CleanUpUnconnectedNodes();
    std::cout << "Finished substructuring\n";
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

/// function that randomly refines some elements
void RandomRefine(TPZGeoMesh *gmesh, TRunConfig &config, int nref)
{
    int nel = config.nelxcoarse*config.nelycoarse;
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
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
