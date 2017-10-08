
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

TAnalyticSolution *example;

void RefinePightCorners(TPZGeoMesh *gmesh);

void RefineParticular(TPZGeoMesh *gmesh,int a,int b, int c);

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    TExceptionManager except;
    
#ifdef _AUTODIFF
    example = new TElasticityExample1;
#endif
    TRunConfig Configuration(4,4,4,2,1,2,0,1,0);  // Config: nelxcoarse, nelycoarse, nhdiv, porderint, nskeldiv, porderskel, hybridize, condensed, LagrangeMult.

    HDivPiola = 1;
    
    // to avoid singular internal matrices
    if (Configuration.numDivSkeleton == Configuration.numHDivisions && Configuration.pOrderInternal <= Configuration.pOrderSkeleton) {
        Configuration.pOrderInternal = Configuration.pOrderSkeleton+1;
    }

    TPZGeoMesh *gmesh = 0;
    TPZVec<long> coarseindices;
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.);
    x1[2] = 0.;
    gmesh = MalhaGeomFredQuadrada(Configuration, x0, x1, coarseindices);

	///Jorge
//    RefineParticular(gmesh,90,94,1800);
	//std::ofstream saida("gmesh_refine1.vtk");
	//TPZVTKGeoMesh::PrintGMeshVTK(gmesh,saida);

//    TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
    TPZAutoPointer<TPZMHMeshControl> MHM;
    
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmesh);
        TPZMHMeshControl *mhm = new TPZMHMeshControl(gmeshauto,coarseindices);
        MHM = mhm;
        TPZMHMeshControl &meshcontrol = *mhm;
        
   //     meshcontrol.SetLagrangeAveragePressure(Configuration.LagrangeMult);
        
        InsertMaterialObjects(meshcontrol);
        
        meshcontrol.SetInternalPOrder(Configuration.pOrderInternal);
        meshcontrol.SetSkeletonPOrder(Configuration.pOrderSkeleton);
        
		// Preenche o vetor de Interfaces do meshcontrol com as interfaces da malha descontinua (que é temporária), para formar o skeleton
		// para isso cria uma malha temporaria discontinua e
        meshcontrol.CreateSkeletonElements(skeleton);
        
        meshcontrol.DivideSkeletonElements(Configuration.numDivSkeleton);
        if(Configuration.Hybridize)
        {
            meshcontrol.Hybridize(secondskeleton, matpressure);
        }
        
        bool substructure = true;
        if (Configuration.Condensed == 0) {
            substructure = false;
        }
        meshcontrol.BuildComputationalMesh(substructure);
        
        ///JORGE
        /// Imprimindo malha computacional
        std::ofstream saida("cmeshh.txt");
        meshcontrol.GMesh()->Reference()->Print(saida);
        
        std::cout << "MHM Computational meshes created\n";

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

void RefineParticular(TPZGeoMesh *gmesh,int a,int b, int c) {
    TPZVec<TPZGeoEl *> subels;
    TPZGeoEl *gel;
    bool cont = true;
    while(cont) {
        gel = gmesh->ElementVec()[a];
        if(gel && !gel->HasSubElement() && gel->Dimension()==2) {
            gel->Divide(subels);
            gel = subels[2];
            gel->Divide(subels);
            cont = false;
        }
        a++;
    }
    cont = true;
    while(cont) {
        gel = gmesh->ElementVec()[b];
        if(gel && !gel->HasSubElement() && gel->Dimension()==2) {
            gel->Divide(subels);
            cont = false;
        }
        b++;
    }
    cont = true;
    while(cont) {
        gel = gmesh->ElementVec()[c];
        if(gel && !gel->HasSubElement() && gel->Dimension()==2) {
            gel->Divide(subels);
            cont = false;
        }
        c++;
    }
    
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
    
    if(example) material1->SetForcingFunction(example->ForcingFunction());
    
    if(example) material1->SetElasticityFunction(example->ConstitutiveLawFunction());
    TPZMaterial * mat1(material1);
    
    TPZMat1dLin *materialCoarse = new TPZMat1dLin(matCoarse);
    TPZFNMatrix<1,STATE> xk(2,2,0.),xb(2,2,0.),xc(2,2,0.),xf(2,1,0.);
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
    
    
    control.SwitchLagrangeMultiplierSign(true);
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
    TPZVec<long> coarseindices;
    
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
    long neq = cmeshauto->NEquations();
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

void RefinePightCorners(TPZGeoMesh *gmesh) {
    int nel = gmesh->NElements();
    TPZVec<TPZGeoEl *> subels;
    TPZVec<TPZGeoEl *> subsubels;
    long el;
    for (el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement() ||  gel->Dimension() != 2) {
            continue;
        }
        if(gel->Node(0).Coord(0)>0.69 && (gel->Node(0).Coord(1)>0.69 || gel->Node(0).Coord(1)<0.26))
            gel->Divide(subels);
    }

    nel = gmesh->NElements();
    for (el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement() ||  gel->Dimension() != 2) {
            continue;
        }
        if(gel->Node(0).Coord(0)>0.8 && (gel->Node(0).Coord(1)>0.69 || gel->Node(0).Coord(1)<0.26))
            gel->Divide(subels);
    }
}
