/**
 * @file
 * @brief Contains Unit Tests for algebraic operations using TPZMatrix<T> derived types
 */

#include "pzfmatrix.h"
#include "pzsfulmat.h"
#include "pzbndmat.h"
#include "pzsbndmat.h"
#include "pzskylnsymmat.h"
#include "pzskylmat.h"
#include "TPZYSMPMatrix.h"
#include "TPZSYSMPMatrix.h"
#ifdef PZ_USING_MKL
#include "TPZYSMPPardiso.h"
#include "TPZSYSMPPardiso.h"
#endif
#include "pzblockdiag.h"
#include "fad.h"


#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>

using namespace std::complex_literals;

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


template <class MAT, class TVar>
void EigenDecompositionAutoFill(int dim, SymProp sp);
template <class TVar>
void BasicEigenTests();
template <class MAT, class TVar>
void GeneralisedEigenvaluesWithAutoFill(int dim, SymProp sp);

#ifdef PZ_USING_LAPACK
/*There is no long double lapack interface in our code*/
TEMPLATE_TEST_CASE("Eigenvalues","[matrix_tests]",
                   float,
                   double,
                   std::complex<float>,
                   std::complex<double>
                   ) {
  using SCAL=TestType;
  using RSCAL=RType(TestType);
  SECTION("Basic tests (TPZFMatrix"){
    BasicEigenTests<SCAL>();
  }

  SymProp sp = GENERATE(SymProp::NonSym,SymProp::Sym,SymProp::Herm);
  SECTION(SymPropName(sp)){
    for(int dim = 3; dim < 10; dim++){
      SECTION("TPZFMatrix"){
        EigenDecompositionAutoFill<TPZFMatrix<SCAL>,SCAL>(dim,sp);
      }
      if(sp==SymProp::Herm){
        //lapack will only compute ev of hermitian band matrices
        SECTION("TPZSBMatrix"){
          EigenDecompositionAutoFill<TPZSBMatrix<SCAL>,SCAL>(dim,sp);
        }
      }
    }
  }
}

TEMPLATE_TEST_CASE("Generalised Eigenvalues","[matrix_tests]",
                   float,
                   double,
                   std::complex<float>,
                   std::complex<double>
                   ) {
  
  using SCAL=TestType;

  SymProp sp = GENERATE(SymProp::NonSym,SymProp::Sym,SymProp::Herm);
  SECTION(SymPropName(sp)){
    for(int dim = 3; dim < 10; dim++){
      SECTION("TPZFMatrix"){
        GeneralisedEigenvaluesWithAutoFill<TPZFMatrix<SCAL>,SCAL>(dim,sp);
      }
      if(sp==SymProp::Herm){
        //lapack will only compute ev of hermitian band matrices
        SECTION("TPZSBMatrix"){
          GeneralisedEigenvaluesWithAutoFill<TPZSBMatrix<SCAL>,SCAL>(dim,sp);
        }
      }
    }
  }
}

// TEMPLATE_TEST_CASE("Singular Value Decomposition (REAL)","[matrix_tests]",
//                 //    float,
//                    double
//                    ) {
//     SingularValueDecomposition_Real<TestType>();
// }
#endif


template <class TVar>
void BasicEigenTests() {
  const auto generalised = GENERATE(0,1);
  constexpr int dim{10};
  TPZFMatrix<TVar> ma(dim,dim);
  TPZFMatrix<TVar> mb(dim,dim);
  mb.Identity();
  TPZFMatrix<CTVar> cpma(dim,dim);

  TPZManVector < CTVar,10 > w;
  TPZFNMatrix < 100,CTVar > eigenVectors;  
  TPZFNMatrix< 10,CTVar > x(dim, 1, 0.);

  SECTION("DiagMatrixRealEigenvectors"){
    ma.Zero();
    cpma.Zero();
    cpma(0,0) = ma(0,0) = 1;
    cpma(dim-1,dim-1) = ma(dim-1,dim-1) = 1;
    
    
    const auto info = [&](){
      TPZFMatrix<TVar> tmpa(ma);
      TPZFMatrix<TVar> tmpb(mb);
      if(generalised){
        return tmpa.SolveGeneralisedEigenProblem(tmpb,w,eigenVectors);
      }
      else return tmpa.SolveEigenProblem(w, eigenVectors);
    }();
    REQUIRE(info == 0);
    TPZFNMatrix< 10,CTVar > actualx1(dim, 1, 0.);
    TPZFNMatrix< 10,CTVar > actualx2(dim, 1, 0.);
    actualx1(0,0) = 1;
    actualx2(dim-1,0) = 1;
    // std::cout<<"--------------------\n";
    // std::cout<<"---------t1---------\n";
    // std::cout<<"--------------------"<<std::endl;
    for (int i = 0; i < dim; i++) {
      eigenVectors.GetSub(0, i, dim, 1, x);
      // std::cout<<"w: "<<w[i]<<std::endl;
      // std::cout<<"v: ";
      // for(int j = 0; j < dim; j++) std::cout<<x(j,0)<<'\t';
      // std::cout<<std::endl;
      if(IsZero(w[i])){
        for(int j = 0; j < dim; j++){
          CAPTURE(i,j,x(j,0));
          REQUIRE(IsZero(x(j,0).imag()));
        }
      }else{
        const auto normres1 = Norm(x-actualx1);
        const bool res1 = IsZero(normres1);
      
        const auto normres2 = Norm(x-actualx2);
        const bool res2 = IsZero(normres2);
      
        CAPTURE(normres1);
        CAPTURE(normres2);
        const bool res = res1 || res2;
        REQUIRE(res);
      }
    }
  }
  
  SECTION("PerturbedIdentityMatrix"){
    const auto epsilon = std::numeric_limits<RTVar>::epsilon()/
      (10*std::numeric_limits<RTVar>::digits10);
    ma.Identity();
    cpma.Identity();
    cpma(dim-1,0) = ma(dim-1,0) = -epsilon;
    cpma(0,dim-1) = ma(0,dim-1) = epsilon;
    
    
    const auto info = [&](){
      TPZFMatrix<TVar> tmpa(ma);
      TPZFMatrix<TVar> tmpb(mb);
      if(generalised){
        return tmpa.SolveGeneralisedEigenProblem(tmpb,w,eigenVectors);
      }
      else return tmpa.SolveEigenProblem(w, eigenVectors);
    }();
    REQUIRE(info == 0);
    TPZFNMatrix< 10,CTVar > res(dim, 1, 0.);
    // std::cout<<"--------------------\n";
    // std::cout<<"---------t2---------\n";
    // std::cout<<"--------------------"<<std::endl;
    for (int i = 0; i < dim; i++) {
      eigenVectors.GetSub(0, i, dim, 1, x);
      // std::cout<<"w: "<<w[i]<<std::endl;
      // std::cout<<"v: ";
      // for(int j = 0; j < dim; j++) std::cout<<x(j,0)<<'\t';
      // std::cout<<std::endl;
      res = cpma * x - w[i] * x;
      const auto norm = Norm(res);
      CAPTURE(i,norm);
      REQUIRE(IsZero(norm));
    }
  }

}
  
template <class MAT, class TVar>
void EigenDecompositionAutoFill(int dim, SymProp sp) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  MAT ma;
  ma.AutoFill(dim, dim, sp);

  TPZFMatrix<CTVar> cpma(dim, dim);
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      cpma(i, j) = ma.Get(i, j);
    }
  }

  TPZVec<CTVar> w;
  TPZFMatrix<CTVar> eigenVectors;
  ma.SolveEigenProblem(w, eigenVectors);

  RTVar mult = 1.;
  if(std::is_same_v<RTVar,float>) {
    mult *= 12.; //This value is arbitrary
  }
  TPZFMatrix<CTVar> x(dim, 1, 0.);
  TPZFMatrix<CTVar> res(dim, 1, 0.);
  for (int i = 0; i < dim; i++) {
    eigenVectors.GetSub(0, i, dim, 1, x);
    res = cpma * x - w[i] * x;
    const auto norm = Norm(res);
    CAPTURE(i, norm);
    REQUIRE(IsZero(norm/mult));
  }
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}

template <class MAT, class TVar>
void GeneralisedEigenvaluesWithAutoFill(int dim, SymProp sp) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  MAT ma,mb;
  ma.AutoFill(dim, dim, sp);
  mb.AutoFill(dim, dim, sp);

  TPZFMatrix<CTVar> cpma(dim, dim), cpmb(dim,dim);
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      cpma(i, j) = ma.Get(i, j);
      cpmb(i, j) = mb.Get(i, j);
    }
  }

  TPZVec<CTVar> w;
  TPZFMatrix<CTVar> eigenVectors;
  ma.SolveGeneralisedEigenProblem(mb,w, eigenVectors);

  RTVar mult = 1;
  if(std::is_same_v<RTVar,float>) {
    mult *= 10.;
  }
  TPZFMatrix<CTVar> x(dim, 1, 0.);
  TPZFMatrix<CTVar> res(dim, 1, 0.);
  for (int i = 0; i < dim; i++) {
    eigenVectors.GetSub(0, i, dim, 1, x);
    res = cpma * x - w[i] * cpmb * x;
    const auto norm = Norm(res);
    CAPTURE(i, norm);
    REQUIRE(IsZero(norm/mult));
  }
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}
