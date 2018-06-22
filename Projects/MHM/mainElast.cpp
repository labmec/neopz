
#include "pzlog.h"
#include "pzgmesh.h"
#include "pzmanvector.h"
#include "TPZMHMeshControl.h"
#include "TPZVTKGeoMesh.h"

#include "pzelasmat.h"
#include "TPZElasticity2DHybrid.h"
#include "pzmat1dlin.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "meshgen.h"

#ifndef USING_MKL
#include "pzskylstrmatrix.h"
#endif
/// Insert material objects for the MHM Mesh solution
void InsertMaterialObjects(TPZMHMeshControl &control);

/// Insert material objects for the MHM Mesh solution
void InsertMaterialObjects(TPZCompMesh &control);

/// Compute an approximation using an H1 approximation
TPZAutoPointer<TPZCompMesh> ComputeH1Approximation(int nelx, int nely, int porder, std::string prefix);

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mainskeleton"));
#endif


const int matInterno = 1;
const int matCoarse = 2;
const int skeleton = 4;
const int secondskeleton = 3;
const int matpressure = 6;

const int dirichlet = 0;

int const bc1=-1;
int const bc2=-2;
int const bc3=-3;
int const bc4=-4;
int const bc5=-5;

TAnalyticSolution *example;
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    TExceptionManager except;
    
#ifdef _AUTODIFF
    example = new TElasticityExample1;
#endif
    TRunConfig Configuration;
    
    /// numhdiv - number of h-refinements
    int j = 6;
    int n_div = 2<<j;
    int j_int = 7 - j;
    int n_div_internal = j_int;
    Configuration.numHDivisions = n_div_internal;
    /// PolynomialOrder - p-order
    Configuration.pOrderInternal = 3;
    Configuration.pOrderSkeleton = 2;
    Configuration.numDivSkeleton = 0;
    Configuration.nelxcoarse = n_div;
    Configuration.nelycoarse = n_div;
    Configuration.Hybridize = 0;
    Configuration.Condensed = 1;
    Configuration.LagrangeMult = 0;
    Configuration.n_threads = 12;

    if (argc == 3)
    {
        Configuration.nelxcoarse = atoi(argv[1]);
        Configuration.nelycoarse = Configuration.nelxcoarse;
        Configuration.pOrderSkeleton = atoi(argv[2]);
        Configuration.pOrderInternal = Configuration.pOrderSkeleton+1;
    }
    HDivPiola = 1;

    if(0)
    {
        /// Compute an approximation using an H1 approximation
        int nx = 40;
        int porder = 2;
        std::string prefix = "H1Aprox";
        
        TPZAutoPointer<TPZCompMesh> cmeshH1 = ComputeH1Approximation(nx, nx, porder, prefix);
    }

    
    // to avoid singular internal matrices
    if (Configuration.numDivSkeleton == Configuration.numHDivisions && Configuration.pOrderInternal <= Configuration.pOrderSkeleton) {
        Configuration.pOrderInternal = Configuration.pOrderSkeleton+1;
    }

    
    TPZGeoMesh *gmesh = 0;
    TPZVec<int64_t> coarseindices;
   
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.);
    x1[2] = 0.;
    int ndiv = Configuration.numHDivisions;
    gmesh = MalhaGeomFredQuadrada(Configuration.nelxcoarse, Configuration.nelycoarse, x0, x1, coarseindices, ndiv);

    TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
    TPZAutoPointer<TPZMHMeshControl> MHM;
    
    
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmesh);
        TPZMHMeshControl *mhm = new TPZMHMeshControl(gmeshauto);
        mhm->DefinePartitionbyCoarseIndices(coarseindices);
        MHM = mhm;
        TPZMHMeshControl &meshcontrol = *mhm;
        
        meshcontrol.SetLagrangeAveragePressure(Configuration.LagrangeMult);
        
        InsertMaterialObjects(meshcontrol);
        
        meshcontrol.SetInternalPOrder(Configuration.pOrderInternal);
        meshcontrol.SetSkeletonPOrder(Configuration.pOrderSkeleton);
        
        meshcontrol.DivideSkeletonElements(Configuration.numDivSkeleton);
        if(Configuration.Hybridize)
        {
            meshcontrol.SetHybridize(true);
        }
        
        bool substructure = true;
        if (Configuration.Condensed == 0) {
            substructure = false;
        }
        meshcontrol.BuildComputationalMesh(substructure);
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream file("GMeshControl.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file);
        }
#endif
#ifdef PZDEBUG
        if(1)
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
            
            std::ofstream out_mhm("MHM_h1.txt");
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

    // compute the MHM solution
    SolveProblem(MHM->CMesh(), MHM->GetMeshes(), example, "MHMElast", Configuration);

    return 0;
}

void InsertMaterialObjects(TPZMHMeshControl &control)
{
    TPZCompMesh &cmesh = control.CMesh();
    /// criar materiais
//    int dim = cmesh.Dimension();
    control.SetProblemType(TPZMHMeshControl::EElasticity2D);
    STATE Young = 1000., nu = 0.3, fx = 0., fy = 0.;
    TPZElasticity2DHybrid *material1 = new TPZElasticity2DHybrid(matInterno,Young,nu,fx,fy);
    material1->SetPlaneStrain();
    control.fMaterialIds.insert(matInterno);
    if(example) material1->SetForcingFunction(example->ForcingFunction());
    
    if(example) material1->SetElasticityFunction(example->ConstitutiveLawFunction());
    TPZMaterial * mat1(material1);
    
    TPZMat1dLin *materialCoarse = new TPZMat1dLin(matCoarse);
    TPZFNMatrix<1,STATE> xk(2,2,0.),xb(2,2,0.),xc(2,2,0.),xf(2,1,0.);
    materialCoarse->SetMaterial(xk, xc, xb, xf);
    
    cmesh.InsertMaterialObject(materialCoarse);
//    materialCoarse = new TPZMat1dLin(skeleton);
//    materialCoarse->SetMaterial(xk, xc, xb, xf);
//    cmesh.InsertMaterialObject(materialCoarse);
//    materialCoarse = new TPZMat1dLin(secondskeleton);
//    materialCoarse->SetMaterial(xk, xc, xb, xf);
//    cmesh.InsertMaterialObject(materialCoarse);
//    materialCoarse = new TPZMat1dLin(matpressure);
//    materialCoarse->SetMaterial(xk, xc, xb, xf);
//    cmesh.InsertMaterialObject(materialCoarse);
    
    
    
    //    REAL diff = -1.;
    //	REAL conv = 0.;
    //	TPZVec<REAL> convdir(3,0.);
    //	REAL flux = 8.;
    //
    //	material1->SetParameters(diff, conv, convdir);
    //	material1->SetInternalFlux( flux);
    
    cmesh.InsertMaterialObject(mat1);
    
    
    control.SwitchLagrangeMultiplierSign(true);
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    //BC -1
    val1(0,0) = 0.;
    val2.Zero();
    val1(0,0) = 0;
    val1(1,1) = 0;
    val2(0,0) = 10.;
    TPZMaterial * BCondD1 = material1->CreateBC(mat1, bc1,dirichlet, val1, val2);
    if(example) BCondD1->SetForcingFunction(example->ValueFunction());
    cmesh.InsertMaterialObject(BCondD1);
    control.fMaterialBCIds.insert(bc1);
    //BC -2
    val1.Zero();
    val2(0,0) = 10.;
    TPZMaterial * BCondD2 = material1->CreateBC(mat1, bc2,dirichlet, val1, val2);
    if(example) BCondD2->SetForcingFunction(example->ValueFunction());
    cmesh.InsertMaterialObject(BCondD2);
    control.fMaterialBCIds.insert(bc2);

    //BC -3
    val1.Zero();
    val2.Zero();
    val2(0,0) = 10.;
    TPZMaterial * BCondD3 = material1->CreateBC(mat1, bc3,dirichlet, val1, val2);
    if(example) BCondD3->SetForcingFunction(example->ValueFunction());
    cmesh.InsertMaterialObject(BCondD3);
    control.fMaterialBCIds.insert(bc3);

    //BC -4
    val1(0,0) = 0;
    val1(1,1) = 1.e9;
    val2(0,0) = 10.;
    TPZMaterial * BCondD4 = material1->CreateBC(mat1, bc4,dirichlet, val1, val2);
    if(example) BCondD4->SetForcingFunction(example->ValueFunction());
    cmesh.InsertMaterialObject(BCondD4);
    control.fMaterialBCIds.insert(bc4);

    //BC -5: dirichlet nulo
    TPZMaterial * BCondD5 = material1->CreateBC(mat1, bc5,dirichlet, val1, val2);
    cmesh.InsertMaterialObject(BCondD5);
    control.fMaterialBCIds.insert(bc5);

}

void InsertMaterialObjects(TPZCompMesh &cmesh)
{
    /// criar materiais
    //    int dim = cmesh.Dimension();
    STATE Young = 1000., nu = 0.3, fx = 0., fy = 0.;
    TPZElasticityMaterial *material1 = new TPZElasticityMaterial(matInterno,Young,nu,fx,fy);
    material1->SetPlaneStrain();
    
    if(example) material1->SetForcingFunction(example->ForcingFunction());
    if(example) material1->SetElasticityFunction(example->ConstitutiveLawFunction());

    TPZMaterial * mat1(material1);
    
    
    //    REAL diff = -1.;
    //	REAL conv = 0.;
    //	TPZVec<REAL> convdir(3,0.);
    //	REAL flux = 8.;
    //
    //	material1->SetParameters(diff, conv, convdir);
    //	material1->SetInternalFlux( flux);
    
    cmesh.InsertMaterialObject(mat1);
    
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    //BC -1
    val1(0,0) = 0.;
    val2.Zero();
    val1(0,0) = 0;
    val1(1,1) = 0;
    TPZMaterial * BCondD1 = material1->CreateBC(mat1, bc1,dirichlet, val1, val2);
    if(example) BCondD1->SetForcingFunction(example->ValueFunction());
    cmesh.InsertMaterialObject(BCondD1);
    
    //BC -2
    val1.Zero();
    val2(0,0) = 10.;
    TPZMaterial * BCondD2 = material1->CreateBC(mat1, bc2,dirichlet, val1, val2);
    if(example) BCondD2->SetForcingFunction(example->ValueFunction());
    cmesh.InsertMaterialObject(BCondD2);
    
    //BC -3
    val1.Zero();
    val2.Zero();
    TPZMaterial * BCondD3 = material1->CreateBC(mat1, bc3,dirichlet, val1, val2);
    if(example) BCondD3->SetForcingFunction(example->ValueFunction());
    cmesh.InsertMaterialObject(BCondD3);
    
    //BC -4
    val1(0,0) = 0;
    val1(1,1) = 1.e9;
    val2(0,0) = -1.;
    TPZMaterial * BCondD4 = material1->CreateBC(mat1, bc4,dirichlet, val1, val2);
    if(example) BCondD4->SetForcingFunction(example->ValueFunction());
    cmesh.InsertMaterialObject(BCondD4);
    
    //BC -5: dirichlet nulo
    TPZMaterial * BCondD5 = material1->CreateBC(mat1, bc5,dirichlet, val1, val2);
    cmesh.InsertMaterialObject(BCondD5);
}

/// Compute an approximation using an H1 approximation
TPZAutoPointer<TPZCompMesh> ComputeH1Approximation(int nelx, int nely, int porder, std::string prefix)
{
    TPZGeoMesh *gmesh = 0;
    TPZVec<int64_t> coarseindices;
    
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.);
    x1[2] = 0.;
    int ndiv = 0;
    gmesh = MalhaGeomFredQuadrada(nelx, nely, x0, x1, coarseindices, ndiv);
    
    gmesh->SetDimension(2);
    TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
    
    TPZAutoPointer<TPZCompMesh> cmeshauto = new TPZCompMesh(gmeshauto);
    cmeshauto->SetDimModel(2);
    InsertMaterialObjects(cmeshauto);
    
    cmeshauto->SetAllCreateFunctionsContinuous();
    
    cmeshauto->AutoBuild();
    
    
    //calculo solution
    bool shouldrenumber = true;
    TPZAnalysis an(cmeshauto,shouldrenumber);
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmeshauto.operator->());
    strmat.SetNumThreads(8);
    
#else
    TPZSkylineStructMatrix strmat(cmeshauto.operator->());
    strmat.SetNumThreads(0);
#endif
    
    
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();
    if(0)
    {
        std::string filename = prefix;
        filename += "_Global.nb";
        std::ofstream global(filename.c_str());
        TPZAutoPointer<TPZStructMatrix> strmat = an.StructMatrix();
        an.Solver().Matrix()->Print("Glob = ",global,EMathematicaInput);
        an.Rhs().Print("Rhs = ",global,EMathematicaInput);
    }
    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
    an.LoadSolution(); // compute internal dofs
                       //    an.Solution().Print("sol = ");
    
    
#ifdef PZDEBUG
    {
        std::ofstream out(prefix+"_MeshWithSol.txt");
        cmeshauto->Print(out);
    }
#endif
    
    std::string configuration;
    {
        std::stringstream sout;
        sout << nelx << "-" << nely;
        configuration = sout.str();
    }
    //    TPZBuildMultiphysicsMesh::TransferFromMeshes(cmeshes, an.Mesh());
    //    for (int i=0; i<cmeshes.size(); i++) {
    //        cmeshes[i]->Solution().Print("sol = ");
    //    }
    //    cmeshes[0]->Solution().Print("solq = ");
    //    cmeshes[1]->Solution().Print("solp = ");
    std::stringstream sout;
    sout << prefix << "Approx-";
    sout << configuration << ".vtk";
    std::string plotfile = sout.str();
    std::cout << "plotfile " << plotfile.c_str() << std::endl;
    TPZStack<std::string> scalnames,vecnames;
    TPZMaterial *mat = cmeshauto->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }
    if (mat->NStateVariables() == 2)
    {
        scalnames.Push("SigmaX");
        scalnames.Push("SigmaY");
        scalnames.Push("TauXY");
        vecnames.Push("Displacement");
    }
    else if(mat->NStateVariables() == 1)
    {
        scalnames.Push("Pressure");
        vecnames.Push("Flux");
        vecnames.Push("Derivative");
    }
    an.DefineGraphMesh(cmeshauto->Dimension(), scalnames, vecnames, plotfile);
    int resolution = 2;
    an.PostProcess(resolution,cmeshauto->Dimension());

#ifdef _AUTODIFF
    std::cout << "Computing errors\n";
    int64_t neq = cmeshauto->NEquations();
    an.SetExact(TElasticityExample1::GradU);
    TPZVec<REAL> errors(3,0.);
    an.PostProcessError(errors);
    std::cout << "Errors computed " << errors << std::endl;

    std::stringstream filename;
    filename << prefix << "Errors.txt";
    std::ofstream out (filename.str(),std::ios::app);
    out << "nelx " << nelx << " nely " << nely << " porder " << porder << " neq " << neq <<  " Energy " << errors[0] << " L2 " << errors[1] << " H1 " << errors[2] << std::endl;
#endif
    return cmeshauto;
    
}

