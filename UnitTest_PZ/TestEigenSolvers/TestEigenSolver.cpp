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


/**solves an EVP and ensures that
lambda - lambda_approx <= residual norm*/
template<class matx, class TVar>
void TestArnoldiSolver(TPZAutoPointer<matx> A, const TPZFMatrix<CTVar> &sol);

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

#ifdef PZ_USING_LAPACK
TEMPLATE_TEST_CASE("Arnoldi Solver 1", "[eigen_tests]",
                   TPZFYsmpMatrix<double>,
                   TPZFMatrix<double>,
                   TPZFMatrix<std::complex<double>>
                   )
{

  /*
    The following matrices are solved:
    [ 1 1 ]
    [ 1 1 ]
    (real eigenpairs: 0, 2)
    [ 1 -1 ]
    [ 1  1 ]
    (complex eigenpairs: 1+j, 1-j)

    for complex types:
    [ 1 i]
    [-i 1]
    (real eigenpairs: 2, 0)
    [ 1 i ]
    [ i 1 ]
    (complex eigenpairs: 1+j, 1-j)
*/
  constexpr int dim{2};
  TPZAutoPointer<TestType> basis_mat = [dim]() -> TestType *{
    if constexpr (std::is_same_v<TPZFMatrix<double>,TestType> ||
                  std::is_same_v<TPZFMatrix<std::complex<double>>,TestType>){
      return new TestType(dim,dim,0);
    }else if constexpr (std::is_same_v<TPZFYsmpMatrix<double>,TestType> ||
                  std::is_same_v<TPZFYsmpMatrix<std::complex<double>>,TestType>){
      //full matrix
      auto mat = new TestType(dim,dim);
      TPZVec<int64_t> ia(dim+1,0), ja(dim*dim,0);
      for(int i = 0; i < dim; i++){
        ia[i] = dim*i;
        for(int j = 0; j < dim; j++){
          ja[dim*i+j] = j;
        }
      }
      ia[dim] = dim*dim;
      TPZVec<typename TestType::Type> a(dim*dim,0);
      mat->SetData(ia,ja,a);
      return mat;
    }else{
      return nullptr;
    }
  }();
  
  SECTION("real_mat_real_eigenpairs"){
    auto mat = basis_mat;
    mat->PutVal(0,0,1);
    mat->PutVal(1,1,1);
    mat->PutVal(0,1,1);
    mat->PutVal(1,0,1);
    TPZFMatrix<CType(typename TestType::Type)> sol(dim,1);
    sol(0,0) = 0;
    sol(1,0) = 2;
    if constexpr (std::is_same_v<typename TestType::Type,double>){
      TestArnoldiSolver<TestType, double>(mat, sol);
    }else{
      TestArnoldiSolver<TestType, std::complex<double>>(mat, sol);
    }
  }
  SECTION("real_mat_cplx_eigenpairs"){
    using namespace std::complex_literals;
    auto mat = basis_mat;
    mat->PutVal(0,0,1);
    mat->PutVal(1,1,1);
    mat->PutVal(0,1,-1);
    mat->PutVal(1,0,1);
    TPZFMatrix<CType(typename TestType::Type)> sol(dim,1);
    sol(0,0) = 1.-1i;
    sol(1,0) = 1.+1i;
    if constexpr (std::is_same_v<typename TestType::Type,double>){
      TestArnoldiSolver<TestType, double>(mat,sol);
    }else{
      TestArnoldiSolver<TestType, std::complex<double>>(mat,sol);
    }
  }
  if constexpr (std::is_same_v<typename TestType::Type,std::complex<double>>){
    using namespace std::complex_literals;
    SECTION("cplx_mat_real_eigenpairs"){
      auto mat = basis_mat;
      mat->PutVal(0,0,1);
      mat->PutVal(1,1,1);
      mat->PutVal(0,1,1i);
      mat->PutVal(1,0,-1i);
      TPZFMatrix<CType(typename TestType::Type)> sol(dim,1);
      sol(0,0) = 2.;
      sol(1,0) = 0.;
      TestArnoldiSolver<TestType, std::complex<double>>(mat,sol);
    }
    SECTION("cplx_mat_cplx_eigenpairs"){
      auto mat = basis_mat;
      mat->PutVal(0,0,1);
      mat->PutVal(1,1,1);
      mat->PutVal(0,1,1i);
      mat->PutVal(1,0,1i);
      TPZFMatrix<CType(typename TestType::Type)> sol(dim,1);
      sol(0,0) = 1.+1i;
      sol(1,0) = 1.-1i;
      TestArnoldiSolver<TestType, std::complex<double>>(mat,sol);
    }
  }
}

TEMPLATE_TEST_CASE("Arnoldi Solver 2", "[eigen_tests]",
                   TPZSkylNSymMatrix<double>,
                   TPZSkylMatrix<double>,
                   TPZFYsmpMatrix<double>,
                   TPZSYsmpMatrix<double>
                   )
{

  TPZAutoPointer<TestType> A = new TestType;
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

  TPZFMatrix<CType(typename TestType::Type)> sol(dim,1,0);
  for(int i = 0; i < dim; i++){sol(i,0) = i + 1;}
  
  TestArnoldiSolver<TestType, typename TestType::Type>(A,sol);
}


#endif

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
void TestArnoldiSolver(TPZAutoPointer<matx> A, const TPZFMatrix<CTVar> &sol)
{

  const auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  
  
  const int64_t dim = A->Rows();
  const int dimKrylov = dim;
  
  TPZFMatrix<CTVar> cpma(dim,dim);
  for(auto i = 0; i < dim; i++){
    for(auto j = 0; j < dim; j++){
      cpma(i,j) = A->GetVal(i,j);
    }
  }
  
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

  const RTVar tol = std::numeric_limits<RTVar>::epsilon()*10;
  for(auto i = 0; i < dimKrylov; i++){
    eigenVectors.GetSub(0, i, dim, 1, x);
    const auto currW = w[i];
    //let us find WHICH eigenvalue was approximated
    RTVar minres = 1e12;
    int minindex=-1;
    
    for(auto ix = 0; ix < dim; ix++){
      const CTVar analyticW = sol[ix];
      const auto res = std::abs(analyticW - currW);
      if(res < minres){
        minres = res;
        minindex = ix;
      }
    }
    const CTVar analyticW = sol[minindex];
    const auto res = Norm(cpma * x - currW * x);

    CAPTURE(res,analyticW,currW);
    REQUIRE((fabs(analyticW - currW) < res || res < tol));
  }
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}
