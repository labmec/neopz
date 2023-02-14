//
// Created by Philippe on 04/02/2023
//
//#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int64_t

#include "pzlog.h"

#ifdef MACOSX
#include <Accelerate/Accelerate.h>
#endif

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "TPZEigenSparseMatrix.h"
#include "TPZGeoMeshTools.h"
#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"
#include "Poisson/TPZMatPoisson.h"
#include "TPZLinearAnalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZVTKGenerator.h"
#include "TPZSpStructMatrix.h"

typedef Eigen::SparseMatrix<double,0,long long> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;


#ifdef PZ_LOG
static TPZLogger logger("pz.eigen");
#endif

// ----- Run tests with or without main -----
#define RUNWITHMAIN

void InvertUsingEigen();
void CreateSparse();
void AccelerateSparse();
void TestSparseClass();
void TestH1Problem();
TPZGeoMesh* CreateGeoMesh(const int dim, TPZVec<int> &nDivs, const int volId, const int bcId);
TPZCompMesh* CreateCMeshH1(TPZGeoMesh* gmesh, const int pOrder, const int volId, const int bcId);


#ifndef RUNWITHMAIN
#include <catch2/catch.hpp>

TEST_CASE("Eigen_inversion", "[eigen_test]") {
    InvertUsingEigen();
}
TEST_CASE("Sparse_mat", "[eigen_test]") {
    CreateSparse();
}
TEST_CASE("Sparse_macos", "[eigen_test]") {
    AccelerateSparse();
}
TEST_CASE("Sparse_class", "[eigen_test]") {
    TestSparseClass();
}

TEST_CASE("H1_simulation", "[eigen_test]") {
    TestH1Problem();
}

#else

int main(){
    TestH1Problem();
    return 0;
}

#endif


void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n)
{
#ifndef RUNWITHMAIN
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(25, 63); // define the range
    for(int i = 0; i<n; i++) {
        double sum = 0, sumabs = 0.;
        for(int j = 0; j<=i; j++) {
//            T a(i,j,distr(gen));
            T a(i,j,-1.);
//            sumabs += abs(a.value());
//            sum += a.value();
            coefficients.push_back(a);
        }
        T diag(i,i,n+1.);
        coefficients.push_back(diag);
        b[i] = 1;
    }
    std::cout << "coefficients\n";
    for(auto it=coefficients.begin(); it != coefficients.end(); it++)
        std::cout << "row " << it->row() << " col " << it->col() << " val " << it->value() << std::endl;
#endif
}

void InvertUsingEigen()
{
    std::cout << "We re all happy\n";
    int n = 3;  // size of the image
   
    // Assembly:
    std::vector<T> coefficients;            // list of non-zeros coefficients
    Eigen::VectorXd b(n);                   // the right hand side-vector resulting from the constraints
    buildProblem(coefficients, b, n);
    std::cout << "B = " << b << std::endl;
   
    SpMat A(n,n);
    A.setFromTriplets(coefficients.begin(), coefficients.end());
    std::cout << "mult " << A*b << std::endl;
    std::cout << " A " << A << std::endl;
    // Solving:
    Eigen::SimplicialCholesky<SpMat> chol(A);  // performs a Cholesky factorization of A
    Eigen::VectorXd x = chol.solve(b);         // use the factorization to solve for the given right hand side
   
    std::cout << "solution " << x << std::endl;
    // Export the result to a file:
//    saveAsBitmap(x, n, argv[1]);
   
    return 0;
}

void CreateSparse()
{
    long long ia[] = {0,1,3};
    long long ja[] = {0,0,1};
    double val[] = {1.,0.,2.};
    double bptr[] = {1.,2.};
    Eigen::Map<Eigen::VectorXd> b(bptr,2);
    b << 1., 2. ;
    int row = 2, col = 2, nnz = 3;
    Eigen::Map< SpMat> sparse(row,col,nnz,ia,ja,val);
    Eigen::SimplicialCholesky<SpMat> chol(sparse);  // performs a Cholesky factorization of A
    Eigen::VectorXd x = chol.solve(b);         // use the factorization to solve for the given right hand side
   
    std::cout << "solution " << x << std::endl;

}

void AccelerateSparse()
{
    long long ia[] = {0,1,3};
    long long ja[] = {0,0,1};
    double val[] = {1.,0.,2.};
    double bptr[] = {1.,2.};
    Eigen::Map<Eigen::VectorXd> b(bptr,2);
    b << 1., 2. ;
    int row = 2, col = 2, nnz = 3;
    Eigen::Map< SpMat > sparse(row,col,nnz,ia,ja,val);
//    Eigen::AccelerateCholeskyAtA<Eigen::Map< SpMat > > sparse(row,col,nnz,ia,ja,val);
    Eigen::SimplicialCholesky<SpMat> chol(sparse);  // performs a Cholesky factorization of A
    Eigen::VectorXd x = chol.solve(b);         // use the factorization to solve for the given right hand side
   
    std::cout << "solution " << x << std::endl;

}

void TestSparseClass() {
    const int64_t nrows = 2;
    int64_t ia[] = {0,1,3};
    int64_t ja[] = {0,0,1};
    STATE val[] = {1.,0.,2.};

    TPZFMatrix<STATE> F(nrows,1);
    F(0,0) = 1.;
    F(1,0) = 4.;
   
    
    TPZEigenSparseMatrix<STATE> spmat(nrows,nrows);
    spmat.SetData(ia, ja, val);
    const DecomposeType dt = ECholesky;
    spmat.Decompose(dt);
    spmat.SolveDirect(F, dt);
    std::cout << "solution " << F << std::endl;
    
#ifndef RUNWITHMAIN
    REQUIRE(F[0] == Approx(1.0));
    REQUIRE(F[1] == Approx(2.0));
#endif
}

void TestH1Problem() {
    
    const int volid = 1, bcid = -1;
    const int ndiv = 2;
    const int dim = 3;
    const int pOrder = 1;
    TPZVec<int> nDivs;
    if(dim == 2) nDivs = {ndiv,ndiv};
    else nDivs = {ndiv,ndiv,ndiv};
    TPZGeoMesh* gmesh = CreateGeoMesh(dim,nDivs,volid,bcid);
    TPZCompMesh* cmesh = CreateCMeshH1(gmesh,pOrder,volid,bcid);
    
    // ========> Solve H1
    TPZLinearAnalysis an(cmesh,false);
    constexpr int nThreads{0};
    TPZStructMatrixT<STATE>* matstruct;
    
    const bool useSparse = true;
    if(useSparse){
        matstruct = new TPZSpStructMatrix<STATE>(cmesh);
    }
    else{
        matstruct = new TPZSkylineStructMatrix<STATE>(cmesh);
    }
//    matstruct = new TPZSkylineStructMatrix<STATE> matstruct(cmesh);
//    TPZSSpStructMatrix<STATE> matstruct(cmesh);
    
    matstruct->SetNumThreads(nThreads);
    an.SetStructuralMatrix(*matstruct);
    TPZStepSolver<STATE> step;
    
//    step.SetDirect(ECholesky);//ELU //ECholesky // ELDLt
    step.SetDirect(ELU);
    an.SetSolver(step);
    
    an.Assemble();
    an.Solve();
    
    // ========> Print
    const std::string plotfile = "postprocess";//sem o .vtk no final
    constexpr int vtkRes{0};
    TPZVec<std::string> fields = {"Solution","Derivative"};
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    vtk.Do();

    // ========> Check if solution is constant 2.0
    TPZFMatrix<STATE>& sol = an.Solution();
    const int64_t nel = sol.Rows();
    STATE average = 0.;
    for(int64_t i = 0; i < nel; i++) average += sol[i];
    average /= nel;
    std::cout << "\n=> Sol average = " << average << std::endl;
#ifndef RUNWITHMAIN
    REQUIRE(average == Approx(2.));    
#endif
    
    delete matstruct;
    delete cmesh;
    delete gmesh;
}

TPZGeoMesh* CreateGeoMesh(const int dim, TPZVec<int> &nDivs, const int volId, const int bcId){
    
    TPZManVector<REAL,3> minX = {0,0,0};
    TPZManVector<REAL,3> maxX = {1,1,1};
    int nMats = 2*dim+1;

    //all bcs share the same id
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,bcId);
    matIds[0] = 1; // volume id
    
    MMeshType meshType;
    if(dim == 2){
        meshType = MMeshType::EQuadrilateral;
    }
    else if(dim == 3) {
        meshType = MMeshType::EHexahedral;
    }
    else{
        DebugStop();
    }
    TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                        matIds, nDivs, meshType,createBoundEls);
    
    std::ofstream out("geomesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    return gmesh;
}

auto ExactSolution = [](const TPZVec<REAL> &loc,TPZVec<STATE>&u,TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];
    const REAL aux = 1./sinh(M_SQRT2*M_PI);
    u[0] = sin(M_PI*x)*sin(M_PI*y)*sinh(M_SQRT2*M_PI*z)/sinh(M_SQRT2*M_PI);
    gradU(0,0) = -M_PI*cos(M_PI*x)*sin(M_PI*y)*sinh(M_SQRT2*M_PI*z)*aux;
    gradU(1,0) = -M_PI*cos(M_PI*y)*sin(M_PI*x)*sinh(M_SQRT2*M_PI*z)*aux;
    gradU(2,0) = -sqrt(2)*M_PI*cosh(M_SQRT2*M_PI*z)*sin(M_PI*x)*sin(M_PI*y)*aux;
};

TPZCompMesh* CreateCMeshH1(TPZGeoMesh* gmesh, const int pOrder, const int volId, const int bcId) {
    
    // ======> Create mesh and set polynomial order
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    const int dim = gmesh->Dimension();
    cmesh->SetDefaultOrder(pOrder);
    
    // ======> Set the volume materials
    TPZMatPoisson<>* poi = new TPZMatPoisson<>(volId,dim);
//    poi->SetExactSol(, ) // Does not need if div(grad(u))=0
    cmesh->InsertMaterialObject(poi);
    
    // ======> Set the boundary conditions
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,2.0);
    const int diri = 0;
    auto * BCond = poi->CreateBC(poi, bcId, diri, val1, val2);
//    BCond->SetForcingFunctionBC(ExactSolution, pOrder+1);
    cmesh->InsertMaterialObject(BCond);
        
    cmesh->SetDimModel(dim);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
    return cmesh;
}
