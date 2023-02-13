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

typedef Eigen::SparseMatrix<double,0,long long> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;


#ifdef PZ_LOG
static TPZLogger logger("pz.eigen");
#endif

// ----- Run tests with or without main -----
//#define RUNWITHMAIN

void InvertUsingEigen();
void CreateSparse();
void AccelerateSparse();
void TestSparseClass();


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
  AccelerateSparse();
}


#else

int main(){
    TestSparseClass();
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
}
