/**
 * @file
 * @brief Contains Basic Unit Tests for TPZMatrix<T> derived types
 */

#include "pzfmatrix.h"
#include "pzsfulmat.h"
#include "pzbndmat.h"
#include "pzsbndmat.h"
#include "pzskylnsymmat.h"
#include "pzskylmat.h"
#include "TPZYSMPMatrix.h"
#include "TPZSYSMPMatrix.h"
#include "pzblockdiag.h"
#include "fad.h"


#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
/**
 * AUXILIARY FUNCTIONS
 **/
//returns true if two matrices are equal
template<class T>
bool CompareMatrices(const TPZMatrix<T> &m1, const TPZMatrix<T> &m2);


TEMPLATE_TEST_CASE("TPZFMatrix copy","[matrix_tests]", CSTATE, REAL, Fad<REAL>) {
  SECTION("Copy constructor"){
    TPZFMatrix<TestType> orig;
    orig.AutoFill(3, 3, SymProp::NonSym);
    TPZFMatrix<TestType> cp(orig);
    REQUIRE(CompareMatrices(cp,orig));
  }
  SECTION("Copy assignment"){
    TPZFMatrix<TestType> orig;
    orig.AutoFill(3,3,SymProp::NonSym);
    TPZFMatrix<TestType> cp;
    SECTION("From empty"){
      cp = orig;
      REQUIRE(CompareMatrices(cp,orig));
    }
      
    SECTION("From diff sizes"){
      const int sz = GENERATE(2,3,5);
      cp.AutoFill(sz,sz,SymProp::NonSym);
      cp = orig;
      REQUIRE(CompareMatrices(cp,orig));
    }
  }
}

TEMPLATE_TEST_CASE("TPZFMatrix move","[matrix_tests]", CSTATE, REAL, Fad<REAL>) {
  SECTION("Move constructor"){
    TPZFMatrix<TestType> orig;
    const int sz = GENERATE(2,3,5);
    orig.AutoFill(sz,sz,SymProp::NonSym);
    TPZFMatrix<TestType> mv(std::move(TPZFMatrix<TestType>(orig)));
    REQUIRE(CompareMatrices(mv,orig));
  }
  
  SECTION("Move assignment"){
    TPZFMatrix<TestType> orig;
    orig.AutoFill(3,3,SymProp::NonSym);
    TPZFMatrix<TestType> mv;
    SECTION("From empty"){
      mv = std::move(TPZFMatrix<TestType>(orig));
      REQUIRE(CompareMatrices(mv,orig));
    }
      
    SECTION("From diff sizes"){
      const int sz = GENERATE(2,3,5);
      mv.AutoFill(sz,sz,SymProp::NonSym);
      mv = std::move(TPZFMatrix<TestType>(orig));
      REQUIRE(CompareMatrices(mv,orig));
    }
  }
}

TEST_CASE("TPZFMatrix list ctor and assignment","[!shouldfail][matrix_tests]") {
  //we dont have any tests written for initializer lists
  SECTION("IMPLEMENT ME"){
    FAIL("IMPLEMENT ME");
  }
}

TEMPLATE_TEST_CASE("TPZFNMatrix copy","[matrix_tests]", REAL) {
  
  //we will try to copy smaller, same size and bigger matrices
  const int orig_sz = GENERATE(2,3,5);
  TPZFMatrix<TestType> orig;
  orig.AutoFill(orig_sz,orig_sz,SymProp::NonSym);
  SECTION("Copy constructor from TPZFMatrix"){
    TPZFNMatrix<9,TestType> cp(orig);
    REQUIRE(CompareMatrices(cp,orig));
  }
  SECTION("Copy assignment from TPZFMatrix"){
    TPZFNMatrix<9,TestType> cp;
    SECTION("From empty"){
      cp = orig;
      REQUIRE(CompareMatrices(cp,orig));
    } 
    SECTION("From diff sizes"){
      const int my_sz = GENERATE(2,3,5);
      cp.AutoFill(my_sz,my_sz,SymProp::NonSym);
      cp = orig;
      REQUIRE(CompareMatrices(cp,orig));
    }
  }
}

TEMPLATE_TEST_CASE("TPZFNMatrix move","[matrix_tests]", CSTATE, REAL, Fad<REAL>) {
  //we will try to move smaller, same size and bigger matrices
  const int orig_sz = GENERATE(2,3,5);
  TPZFMatrix<TestType> orig;
  orig.AutoFill(orig_sz,orig_sz,SymProp::NonSym);
  SECTION("Move constructor from TPZFMatrix"){
    TPZFNMatrix<9,TestType> cp(std::move(TPZFMatrix<TestType>(orig)));
    REQUIRE(CompareMatrices(cp,orig));
  }
  SECTION("Move assignment from TPZFMatrix"){
    TPZFNMatrix<9,TestType> cp;
    SECTION("From empty"){
      cp = std::move(TPZFMatrix<TestType>(orig));
      REQUIRE(CompareMatrices(cp,orig));
    } 
    SECTION("From diff sizes"){
      const int my_sz = GENERATE(2,3,5);
      cp.AutoFill(my_sz,my_sz,SymProp::NonSym);
      cp = std::move(TPZFMatrix<TestType>(orig));
      REQUIRE(CompareMatrices(cp,orig));
    }
  }
  SECTION("Move constructor from TPZFNMatrix using static mem(bigger)"){
    TPZFNMatrix<9,TestType> cp(std::move(TPZFNMatrix<16,TestType>(orig)));
    REQUIRE(CompareMatrices(cp,orig));
  }
  SECTION("Move constructor from TPZFNMatrix using static mem(exact)"){
    TPZFNMatrix<9,TestType> cp(std::move(TPZFNMatrix<9,TestType>(orig)));
    REQUIRE(CompareMatrices(cp,orig));
  }
  SECTION("Move constructor from TPZFNMatrix using dynamic mem"){
    TPZFNMatrix<9,TestType> cp(std::move(TPZFNMatrix<4,TestType>(orig)));
    REQUIRE(CompareMatrices(cp,orig));
  }
  SECTION("Move assignment from TPZFNMatrix using static mem(bigger)"){
    TPZFNMatrix<9,TestType> cp;
    SECTION("From empty"){
      cp = std::move(TPZFNMatrix<16,TestType>(orig));
      REQUIRE(CompareMatrices(cp,orig));
    } 
    SECTION("From diff sizes"){
      const int my_sz = GENERATE(2,3,5);
      cp.AutoFill(my_sz,my_sz,SymProp::NonSym);
      cp = std::move(TPZFNMatrix<16,TestType>(orig));
      REQUIRE(CompareMatrices(cp,orig));
    }
  }
  SECTION("Move assignment from TPZFNMatrix using static mem(exact)"){
    TPZFNMatrix<9,TestType> cp;
    SECTION("From empty"){
      cp = std::move(TPZFNMatrix<9,TestType>(orig));
      REQUIRE(CompareMatrices(cp,orig));
    } 
    SECTION("From diff sizes"){
      const int my_sz = GENERATE(2,3,5);
      cp.AutoFill(my_sz,my_sz,SymProp::NonSym);
      cp = std::move(TPZFNMatrix<9,TestType>(orig));
      REQUIRE(CompareMatrices(cp,orig));
    }
  }
  SECTION("Move assignment from TPZFNMatrix using dynamic mem"){
    TPZFNMatrix<9,TestType> cp;
    SECTION("From empty"){
      cp = std::move(TPZFNMatrix<4,TestType>(orig));
      REQUIRE(CompareMatrices(cp,orig));
    } 
    SECTION("From diff sizes"){
      const int my_sz = GENERATE(2,3,5);
      cp.AutoFill(my_sz,my_sz,SymProp::NonSym);
      cp = std::move(TPZFNMatrix<4,TestType>(orig));
      REQUIRE(CompareMatrices(cp,orig));
    }
  }
}

TEST_CASE("Resize tests","[matrix_tests][!shouldfail]"){
  SECTION("IMPLEMENT ME"){
    FAIL("IMPLEMENT ME");
  }
}

TEMPLATE_PRODUCT_TEST_CASE("Generate hermitian matrix","[matrix_tests]",
                           (TPZFMatrix,TPZSFMatrix,
                            TPZFBMatrix,TPZSBMatrix,
                            TPZSkylNSymMatrix,TPZSkylMatrix,
                            TPZFYsmpMatrix,TPZSYsmpMatrix),(STATE,CSTATE))
{
  using SCAL = typename TestType::Type;
  const auto nrows = 10;
  const auto ncols = 10;
  TestType mat;
  mat.AutoFill(nrows,ncols,SymProp::Herm);
  for (auto i = 0; i < nrows; i++) {
    auto j = i;
    if constexpr (is_complex<SCAL>::value){
      CAPTURE(i,j,mat.Get(i,j));
      REQUIRE(mat.Get(i,j).imag() == (RType(SCAL))0);
    }
    j++;
    for(; j < ncols; j++){
      CAPTURE(i,j,mat.Get(i,j));
      REQUIRE(mat.Get(i,j)-std::conj(mat.Get(j,i)) == (SCAL)0);
    }
  }
  REQUIRE(mat.VerifySymmetry()==SymProp::Herm);
}

TEMPLATE_PRODUCT_TEST_CASE("Generate symmetric matrix","[matrix_tests]",
                           (TPZFMatrix,TPZSFMatrix,
                            TPZFBMatrix,TPZSBMatrix,
                            TPZSkylNSymMatrix,TPZSkylMatrix,
                            TPZFYsmpMatrix,TPZSYsmpMatrix),(STATE,CSTATE))
{
  using SCAL = typename TestType::Type;
  const auto nrows = 10;
  const auto ncols = 10;
  TestType mat;
  mat.AutoFill(nrows,ncols,SymProp::Sym);
  for (auto i = 0; i < nrows; i++) {
    for(auto j=i; j < ncols; j++){
      CAPTURE(i,j,mat.Get(i,j),mat.Get(j,i));
      REQUIRE(mat.Get(i,j)-mat.Get(j,i) == (SCAL)0);
    }
  }
  const auto sp = mat.VerifySymmetry();
  if constexpr (is_complex<SCAL>::value){
    REQUIRE(sp==SymProp::Sym);
  }else{
    //it makes no difference for real matrices
    REQUIRE((sp==SymProp::Sym || sp==SymProp::Herm));
  }
}

TEMPLATE_PRODUCT_TEST_CASE("Generate diag dominant matrix","[matrix_tests]",
                           (TPZFMatrix,TPZSFMatrix,
                            TPZFBMatrix,TPZSBMatrix,
                            TPZSkylNSymMatrix,TPZSkylMatrix,
                            TPZFYsmpMatrix,TPZSYsmpMatrix),(STATE,CSTATE))
{
  const auto nrows = 10;
  const auto ncols = 10;
  TestType mat;
  SymProp sp = GENERATE(SymProp::Sym, SymProp::Herm, SymProp::NonSym);
  if constexpr(std::is_same_v<TestType,TPZSFMatrix<STATE>> || std::is_same_v<TestType,TPZSFMatrix<CSTATE>> ||
               std::is_same_v<TestType,TPZSBMatrix<STATE>> || std::is_same_v<TestType,TPZSBMatrix<CSTATE>> ||
               std::is_same_v<TestType,TPZSkylMatrix<STATE>> || std::is_same_v<TestType,TPZSkylMatrix<CSTATE>> ||
               std::is_same_v<TestType,TPZSYsmpMatrix<STATE>> || std::is_same_v<TestType,TPZSYsmpMatrix<CSTATE>>
               ){
    if(sp==SymProp::NonSym){
      return;
    }
  }
  mat.AutoFill(nrows,ncols,sp);

  REAL sum{0};
  for (int i = 0; i < mat.Rows(); i++) {
    sum = 0.0;
    for (int j = 0; j < mat.Cols(); j++) {
      if (i != j)
        sum += fabs(mat.Get(i, j));
    }
    CAPTURE(i);
    CAPTURE(SymPropName(sp));
    CAPTURE(fabs(mat.Get(i,i)));
    CAPTURE(sum);
    REQUIRE(fabs(mat.Get(i, i)) > sum);
  }
}
TEMPLATE_TEST_CASE("Symmetry test","[matrix_tests]",
                   TPZSFMatrix<CSTATE>,
                   TPZSBMatrix<CSTATE>,
                   TPZSkylMatrix<CSTATE>,
                   TPZSYsmpMatrix<CSTATE>){
  constexpr int nrows = 3;
  TPZAutoPointer<TPZMatrix<CSTATE>> mat{nullptr};
  if constexpr(std::is_same_v<TestType,TPZSFMatrix<CSTATE>>){
    auto mymat = new TPZSFMatrix<CSTATE>(nrows);
    mat = mymat;
  } else if constexpr(std::is_same_v<TestType,TPZSBMatrix<CSTATE>>){
    constexpr int band = 2;
    auto mymat = new TPZSBMatrix<CSTATE>(nrows,band);
    mat = mymat;
  } else if constexpr(std::is_same_v<TestType,TPZSkylMatrix<CSTATE>>){
    auto mymat = new TPZSkylMatrix<CSTATE>(nrows);
    TPZVec<int64_t> skylvec(nrows,0);
    mymat->SetSkyline(skylvec);
    mat = mymat;
  } else if constexpr(std::is_same_v<TestType,TPZSYsmpMatrix<CSTATE>>){
    auto mymat = new TPZSYsmpMatrix<CSTATE>(nrows,nrows);
    TPZVec<int64_t> ia = {0,3,5,6}, ja = {0,1,2,1,2,2};
    TPZVec<CSTATE> aa = {0,0,0,0,0,0};
    mymat->SetData(ia,ja,aa);
    mat = mymat;
  }else {
    FAIL("Expected different matrix type");
  }
  SECTION("Symmetric"){
    mat->SetSymmetry(SymProp::Sym);
    REQUIRE(mat->GetSymmetry() == SymProp::Sym);
    const auto val = CSTATE(0,1);
    mat->PutVal(2,0,val);
    auto cval = mat->GetVal(0,2);
    REQUIRE(val == cval);
  }
  SECTION("Hermitian"){
    mat->SetSymmetry(SymProp::Herm);
    REQUIRE(mat->GetSymmetry() == SymProp::Herm);
    const auto val = CSTATE(0,1);
    mat->PutVal(2,0,val);
    auto cval = mat->GetVal(0,2);
    REQUIRE(val == std::conj(cval));
  }
}

TEMPLATE_PRODUCT_TEST_CASE("Transpose tests","[matrix_tests]",
                           (TPZFMatrix,TPZBlockDiagonal,TPZFBMatrix),(STATE,CSTATE)){
  using MAT = TestType;
  using SCAL = typename TestType::Type;

  constexpr int nr{10};
  const int nc = GENERATE(5,10,15);
  const SymProp sp = GENERATE(SymProp::Sym, SymProp::Herm, SymProp::NonSym);
  MAT ma;
  //only square matrices can be symmetric/hermitian
  if(sp !=SymProp::NonSym && nc!=10){return;}
  
  ma.AutoFill(nr,nc,sp);

  TPZFMatrix<SCAL> matrans(nc,nr,0),matranstrans(nr,nc,0);

  bool conj = GENERATE(false,true);

  ma.Transpose(&matrans,conj);


  for(int dummy = 0; dummy < 1; dummy++){//stupid for just so we can use break
    if constexpr (std::is_same_v<SCAL,CSTATE>){
      if(conj){
        for(int ir = 0; ir < nr; ir++){
          for(int ic = 0; ic < nc; ic++){
            CAPTURE(ir,ic,ma.Get(ir,ic), matrans.Get(ic,ir));
            REQUIRE(ma.Get(ir,ic)-std::conj(matrans.Get(ic,ir)) == (SCAL)0);
          }
        }
        break;
      }
    }
    for(int ir = 0; ir < nr; ir++){
      for(int ic = 0; ic < nc; ic++){
        CAPTURE(ir,ic,ma.Get(ir,ic), matrans.Get(ic,ir));
        REQUIRE(ma.Get(ir,ic)-matrans.Get(ic,ir) == (SCAL)0);
      }
    }
  }
  
  
  matrans.Transpose(&matranstrans,conj);
  for(int ir = 0; ir < nr; ir++){
    for(int ic = 0; ic < nc; ic++){
      CAPTURE(ir,ic,ma.Get(ir,ic), matranstrans.Get(ir,ic));
      REQUIRE(ma.Get(ir,ic)-matranstrans.Get(ir,ic) == (SCAL)0);
    }
  }  
}



template<class T>
bool CompareMatrices(const TPZMatrix<T> &m1, const TPZMatrix<T> &m2){
  if(m1.Rows()!=m2.Rows() || m1.Cols() != m2.Cols()){
    return false;
  }
  const int nr = m1.Rows();
  const int nc = m1.Cols();
  for(int ic = 0; ic < nc; ic++){
    for(int ir = 0; ir < nr; ir++){
      if(m1.Get(ic,ir)!=m2.Get(ic,ir)){
        return false;
      }
    }
  }
  return true;
}