/**
 * @file
 * @brief Contains Unit Tests for methods of the matrices classes.
 */

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzblockdiag.h"
#include "pzbndmat.h"
#include "pzsbndmat.h"
#include "pzsfulmat.h"
#include "pzskylnsymmat.h"
#include "pzskylmat.h"
#include "pzysmp.h"
#include "pzsysmp.h"
#include "TPZKrylovEigenSolver.h"

#include <catch2/catch.hpp>

template<>
struct Catch::StringMaker<long double> {
    static std::string convert(long double ref);
    static int precision;
};

int Catch::StringMaker<long double>::precision = 10;

std::string Catch::StringMaker<long double>::convert(long double value) {
    std::ostringstream out;
    out.precision(precision);
    out << std::fixed << value;
    return out.str();
}

/**tests if the H matrix is Hessenberg, 
and if the vectors are an orthonormal basis*/
template<class matx, class TVar>
void TestArnoldiIteration(bool sym);


/**solves a dummy EVP using the arnoldi solver.
for hermitian matrices it is ensured that 
lambda - lambda_approx <= residual norm*/
template<class matx, class TVar>
void TestArnoldiSolver();

TEMPLATE_TEST_CASE("Arnoldi Iteration", "[eigen_tests]",
                   double)
{
  bool isSym{true};
  SECTION("Sym")
    {
      SECTION("TPZFMatrix")
        {TestArnoldiIteration<TPZFMatrix<TestType>,TestType>(isSym);}
      SECTION("TPZSkylNSymMatrix")
        {TestArnoldiIteration<TPZSkylNSymMatrix<TestType>,TestType>(isSym);}
      SECTION("TPZSkylMatrix")
        {TestArnoldiIteration<TPZSkylMatrix<TestType>,TestType>(isSym);}
      SECTION("TPZFYsmpMatrix")
        {TestArnoldiIteration<TPZFYsmpMatrix<TestType>,TestType>(isSym);}
      SECTION("TPZSYsmpMatrix")
        {TestArnoldiIteration<TPZSYsmpMatrix<TestType>,TestType>(isSym);}
    }
  isSym = false;
  
  SECTION("NSym")
    {
      SECTION("TPZFMatrix")
        {TestArnoldiIteration<TPZFMatrix<TestType>,TestType>(isSym);}
      SECTION("TPZSkylNSymMatrix")
        {TestArnoldiIteration<TPZSkylNSymMatrix<TestType>,TestType>(isSym);}
      SECTION("TPZFYsmpMatrix")
        {TestArnoldiIteration<TPZFYsmpMatrix<TestType>,TestType>(isSym);}
    }
}

TEMPLATE_TEST_CASE("Arnoldi Solver", "[eigen_tests]",
                   double)
{
  
  SECTION("TPZSkylNSymMatrix")
    {TestArnoldiSolver<TPZSkylNSymMatrix<TestType>,TestType>();}
  SECTION("TPZSkylMatrix")
    {TestArnoldiSolver<TPZSkylMatrix<TestType>,TestType>();}
  SECTION("TPZFYsmpMatrix")
    {TestArnoldiSolver<TPZFYsmpMatrix<TestType>,TestType>();}
  SECTION("TPZSYsmpMatrix")
    {TestArnoldiSolver<TPZSYsmpMatrix<TestType>,TestType>();}

}

template<class matx, class TVar>
void TestArnoldiIteration(bool sym)
{

  const auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  
  matx A;
  constexpr int dim{100};
  int dimKrylov{100};
  A.AutoFill(dim,dim,sym);
  TPZKrylovEigenSolver<TVar> arnoldi;
  arnoldi.SetKrylovDim(dimKrylov);
  TPZVec<TPZAutoPointer<TPZFMatrix<TVar>>> qVecs;
  TPZFMatrix<TVar> hMat;
  bool success = arnoldi.ArnoldiIteration(A,qVecs,hMat);
  
  REQUIRE(success);
  
  constexpr RTVar tol = 1e6*std::numeric_limits<RTVar>::epsilon();
  for(auto i = 0; i < dimKrylov; i++){
    const auto & qi = *qVecs[i];
    const auto norm = Norm(qi);
    CAPTURE(i,norm);
    REQUIRE(norm == Approx(1.0));
    for(auto j = 0; j < i; j++){
      const auto &qj = *qVecs[j];
      const auto qiqj = Dot(qi,qj);
      CAPTURE(j,qiqj);
      REQUIRE(qiqj == Approx(0.0).margin(tol));
    }
  }
  for(auto i = 0; i < dimKrylov; i++){
    for(auto j = 0; j < i-1; j++){
      const auto val = hMat.GetVal(i,j);
      CAPTURE(i,j,val);
      REQUIRE(val == Approx(0.0));
    }
  }

  dimKrylov = dim;
  arnoldi.SetKrylovDim(dimKrylov);
  success = arnoldi.ArnoldiIteration(A,qVecs,hMat);
  REQUIRE(success);
  TPZFMatrix<TVar> qMat(dim,dimKrylov);
  for(int i = 0; i < dimKrylov; i++){
    qMat.PutSub(0,i,*qVecs[i]);
  }
  TPZFMatrix<TVar> aMat = qMat * hMat;
  qMat.Transpose();
  aMat = aMat * qMat;
  std::cout<< "A = ("<<A.Rows()<<","<<A.Cols()<<")\n";
  std::cout<< "a = ("<<aMat.Rows()<<","<<aMat.Cols()<<")\n";
  for(auto i = 0; i < dim; i++)
    for(auto j = 0; j < dim; j++){
      CAPTURE(i,j,aMat(i,j),A.GetVal(i,j));
      REQUIRE(aMat(i,j)-A.GetVal(i,j)== Approx(0.0).margin(tol));
    }
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}

template<class matx, class TVar>
void TestArnoldiSolver()
{

  const auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  
  TPZAutoPointer<matx> A = new matx;
  constexpr int64_t dim{20};
  constexpr int dimKrylov{20};
  constexpr bool isSym{true};
  A->AutoFill(dim,dim,isSym);//generates adequate storage
  
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      if(!IsZero(A->GetVal(i,j))) A->PutVal(i,j,0);
    }
    A->PutVal(i,i,i+1);
  }
  TPZFMatrix<CTVar> cpma(dim,dim);
  for(int i = 0; i < dim; i++)
    cpma.PutVal(i,i,i+1);
  
  TPZKrylovEigenSolver<TVar> arnoldi;
  arnoldi.SetKrylovDim(dimKrylov);
  arnoldi.SetNEigenpairs(dimKrylov);
  arnoldi.SetMatrixA(A);
  TPZManVector <CTVar,10> w;
  TPZFMatrix <CTVar> eigenVectors;
  
  const bool success = arnoldi.SolveEigenProblem(w,eigenVectors) == 0;
  REQUIRE(success);

  const RTVar mult = sizeof(RTVar) == 4 ? 100 : 10;
  TPZFMatrix<CTVar> x(dim,1);
  for(auto i = 0; i < dimKrylov; i++) std::cout<<w[i]<<std::endl;

  
  for(auto i = 0; i < dimKrylov; i++){
    eigenVectors.GetSub(0, i, dim, 1, x);
    const CTVar analyticW = i+1;
    const auto res = Norm(cpma * x - w[i] * x);
    CAPTURE(res,analyticW,w[i]);
    REQUIRE(fabs(analyticW - w[i]) < res);
  }
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}
