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


template <class MAT>
void TestInverse(const bool sym_storage, const SymProp sp, const DecomposeType dec);
template<class MAT>
void TestAdd(const MAT &mat, const bool use_operator);
template<class MAT>
void TestSubtract(const MAT &mat, const bool use_operator);

template<class MAT>
struct SymmetricStorage : std::false_type {};
template<class TVar>
struct SymmetricStorage<TPZSFMatrix<TVar>> : std::true_type {};
template<class TVar>
struct SymmetricStorage<TPZSBMatrix<TVar>> : std::true_type {};
template<class TVar>
struct SymmetricStorage<TPZSkylMatrix<TVar>> : std::true_type {};
template<class TVar>
struct SymmetricStorage<TPZSYsmpMatrix<TVar>> : std::true_type {};
#ifdef PZ_USING_MKL
template<class TVar>
struct SymmetricStorage<TPZSYsmpMatrixPardiso<TVar>> : std::true_type {};

template<class MAT>
inline constexpr bool IsSymmetricStorage = SymmetricStorage<MAT>::value;
#endif

TEMPLATE_PRODUCT_TEST_CASE("Inverse (nsym,REAL)","[matrix_tests]",
                           (
                             TPZFMatrix,
                             TPZFBMatrix,
                             TPZSkylNSymMatrix,
#ifdef PZ_USING_MKL
                             TPZFYsmpMatrixPardiso,
#endif
                             TPZBlockDiagonal
                            ),
                           (float,
                            double,
                            long double)
                           ) {

  using MAT = TestType;
  using SCAL = typename MAT::Type;
#ifdef PZ_USING_MKL
  //pardiso is only for doubles
  if(
    std::is_same_v<MAT,TPZFYsmpMatrixPardiso<float>> ||
    std::is_same_v<MAT,TPZFYsmpMatrixPardiso<long double>>||
    std::is_same_v<MAT,TPZFYsmpMatrixPardiso<std::complex<float>>>||
    std::is_same_v<MAT,TPZFYsmpMatrixPardiso<std::complex<long double>>>){return;}
#endif
  DecomposeType dec = GENERATE(ELU,ECholesky,ELDLt);
  SECTION(DecomposeTypeName(dec)){
    SymProp sp = GENERATE(SymProp::NonSym,SymProp::Sym,SymProp::Herm);
    SECTION(SymPropName(sp)){
      if((dec==ECholesky && sp==SymProp::NonSym) ||
         (dec==ELDLt && sp==SymProp::NonSym)){
        return;
      }
      if(
        (std::is_same_v<MAT,TPZFBMatrix<SCAL>> ||
         std::is_same_v<MAT,TPZSkylNSymMatrix<SCAL>> ||
         std::is_same_v<MAT,TPZBlockDiagonal<SCAL>>) &&
        dec != ELU){
        return;
      }
  
      TestInverse<MAT>(false,sp,dec);
    }
  }
}

template <class MAT>
void TestInverse();
TEMPLATE_PRODUCT_TEST_CASE("Inverse (nsym,CPLX)","[matrix_tests]",
                           (
                             TPZFMatrix,
                             TPZFBMatrix,
                             TPZSkylNSymMatrix,
#ifdef PZ_USING_MKL
                             TPZFYsmpMatrixPardiso,
#endif
                             TPZBlockDiagonal
                            ),
                           (std::complex<float>,
                            std::complex<double>,
                            std::complex<long double>)
                           ) {
  using MAT = TestType;
  using SCAL = typename MAT::Type;
  
#ifdef PZ_USING_MKL
  //pardiso is only for doubles
  if(
    std::is_same_v<TestType,TPZFYsmpMatrixPardiso<float>> ||
    std::is_same_v<TestType,TPZFYsmpMatrixPardiso<long double>>||
    std::is_same_v<TestType,TPZFYsmpMatrixPardiso<std::complex<float>>>||
    std::is_same_v<TestType,TPZFYsmpMatrixPardiso<std::complex<long double>>>){return;}
#endif
  DecomposeType dec = GENERATE(ELU,ECholesky,ELDLt);
  SECTION(DecomposeTypeName(dec)){
    SymProp sp = GENERATE(SymProp::NonSym,SymProp::Sym,SymProp::Herm);
    SECTION(SymPropName(sp)){
      if((dec==ECholesky && sp==SymProp::NonSym) ||
         (dec==ELDLt && sp==SymProp::NonSym) ||
         (is_complex<SCAL>::value && dec!=ELU && sp==SymProp::Sym)){
        return;
      }
      //some more
      if(
        (std::is_same_v<MAT,TPZFBMatrix<SCAL>> ||
         std::is_same_v<MAT,TPZSkylNSymMatrix<SCAL>> ||
         std::is_same_v<MAT,TPZBlockDiagonal<SCAL>>) &&
        dec != ELU){
        return;
      }
      TestInverse<MAT>(false,sp,dec);
    }
  }
}

template <class MAT>
void TestInverse();
TEMPLATE_PRODUCT_TEST_CASE("Inverse (sym,REAL)","[matrix_tests]",
                           (
                             TPZSFMatrix,
                             TPZSBMatrix,
                             TPZSkylMatrix,
#ifdef PZ_USING_MKL
                             TPZSYsmpMatrixPardiso
#endif
                            ),
                           (float,
                            double,
                            long double)
                           ) {
  using MAT = TestType;
  using SCAL = typename MAT::Type;
#ifdef PZ_USING_MKL
  //pardiso is only for doubles
  if(
    std::is_same_v<TestType,TPZSYsmpMatrixPardiso<float>> ||
    std::is_same_v<TestType,TPZSYsmpMatrixPardiso<long double>>||
    std::is_same_v<TestType,TPZSYsmpMatrixPardiso<std::complex<float>>>||
    std::is_same_v<TestType,TPZSYsmpMatrixPardiso<std::complex<long double>>>){return;}
#endif
  DecomposeType dec = GENERATE(ECholesky,ELDLt);
  SECTION(DecomposeTypeName(dec)){
    const SymProp sp = SymProp::Herm;
    TestInverse<MAT>(true,sp,dec);
  }
}

template <class MAT>
void TestInverse();


//for debugging: once this is fixed, uncomment next test
TEST_CASE("Inverse TPZSBMatrix cplx","[matrix_tests][!shouldfail]"){
  SECTION("Find out what is happening"){
    using MAT=TPZSBMatrix<std::complex<double>>;
    const SymProp sp=SymProp::Herm;
    MAT ma1;
    ma1.AutoFill(10,10,sp);
    TestInverse<MAT>(true,sp,ELDLt);
  }
}
TEMPLATE_PRODUCT_TEST_CASE("Inverse (sym,CPLX)","[matrix_tests]",
                           (
                             TPZSFMatrix,
                             // TPZSBMatrix,
                             TPZSkylMatrix,
#ifdef PZ_USING_MKL
                             TPZSYsmpMatrixPardiso
#endif
                            ),
                           (std::complex<float>,
                            std::complex<double>,
                            std::complex<long double>)
                           ) {
  using MAT = TestType;
  using SCAL = typename MAT::Type;
#ifdef PZ_USING_MKL
  //pardiso is only for doubles
  if(
    std::is_same_v<TestType,TPZSYsmpMatrixPardiso<float>> ||
    std::is_same_v<TestType,TPZSYsmpMatrixPardiso<long double>>||
    std::is_same_v<TestType,TPZSYsmpMatrixPardiso<std::complex<float>>>||
    std::is_same_v<TestType,TPZSYsmpMatrixPardiso<std::complex<long double>>>){return;}
#endif
  DecomposeType dec = GENERATE(ECholesky,ELDLt);
  SECTION(DecomposeTypeName(dec)){
    SymProp sp = GENERATE(SymProp::Sym, SymProp::Herm);
    SECTION(SymPropName(sp)){
#ifdef PZ_USING_MKL
      //we can only solve symmetric complex matrices with sym storage using pardiso
      if(sp==SymProp::Sym && !std::is_same_v<MAT,TPZSYsmpMatrixPardiso<std::complex<double>>>){
        return;
      }
#endif
      TestInverse<MAT>(true,sp,dec);
    }
  }
}

TEMPLATE_TEST_CASE("DotNorm","[matrix_tests]",float,double,long double,
                   std::complex<float>,std::complex<double>,std::complex<long double>){
  using TVar=TestType;
  constexpr int dim{10};
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  TPZFMatrix<TVar> ma1(dim,1,0);
  TPZFMatrix<TVar> ma2(dim,1,0);
  
  for(int i = 0; i < dim; i++) ma1.PutVal(i,0,i+1);

  SECTION("dot product with zero vec"){
    TVar res = Dot(ma1,ma2);
    REQUIRE(res == (TVar)0.);
  }
  
  for(int i = 0; i < dim; i++) ma2.PutVal(i,0,i+1);

  const auto res_check = [dim]{
    TVar res = 0;
    for(int i = 0; i < dim; i++) res+= (i+1)*(i+1);
    return res;
  }();
  
  SECTION("dot vs norm comparison"){
    TVar res = Dot(ma1,ma2);
    REQUIRE(res == res_check);
  }
  
  //complex tests
  if constexpr (std::is_same_v<TVar,CTVar>){
    using namespace std::complex_literals;
    CAPTURE(dim);

    constexpr RTVar tol = std::numeric_limits<RTVar>::epsilon()/
      (10*std::numeric_limits<RTVar>::digits10);

    for(int i = 0; i < dim; i++){
      const TVar val = (TVar)(i+1);
      ma1.PutVal(i,0,val*(TVar)1.0i);
      ma2.PutVal(i,0,val);
    }
    const auto res = Dot(ma1,ma2);
    CAPTURE(res);
    CAPTURE(res_check);
    REQUIRE(res == res_check * (TVar)1i);
    
    for(int i = 0; i < dim; i++){
      const TVar val = (TVar)(i+1) * (TVar) 1i;
      ma2.PutVal(i,0, val);
      ma2.PutVal(i,0,-val);
    }
    //ma1 is the conjugate of ma2
    const auto dot1 = Dot(ma1,ma2);
    const auto dot2 = Dot(ma2,ma1);
    const auto norm1 = Norm(ma1);
    const auto norm2 = Norm(ma2);
    CAPTURE(dot1);
    CAPTURE(dot2);
    REQUIRE(std::abs(dot1) == Catch::Approx(norm1*norm1).margin(tol));
    REQUIRE(std::abs(dot1) == Catch::Approx(norm2*norm2).margin(tol));
    REQUIRE(std::abs(dot1-dot2) == Catch::Approx(0).margin(tol));
  }
  
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}


//for debugging: once this is fixed, uncomment next test
TEST_CASE("Add TPZSkylNSymMatrix","[matrix_tests][!shouldfail]"){
  SECTION("With operator it executes, without it doesnt"){
    TPZSkylNSymMatrix<double> ma1;
    ma1.AutoFill(10,10,SymProp::NonSym);
    TestAdd(ma1,false);
    TestAdd(ma1,true);
  }
}
TEMPLATE_PRODUCT_TEST_CASE("Add(nsym)","[matrix_tests]",
                           (TPZFMatrix,TPZFBMatrix,TPZFYsmpMatrix
                            // ,TPZSkylNSymMatrix
                            ),(float,double,long double,
                               std::complex<float>,std::complex<double>,
                               std::complex<long double>)){
  using MAT=TestType;
  using SCAL=typename MAT::Type;
  const SymProp sp=SymProp::NonSym;
  constexpr int row{10};
  const int col = GENERATE(5,10,15);
  if(std::is_same_v<MAT,TPZFBMatrix<SCAL>> && col!=10){return;}
  MAT ma1;
  ma1.AutoFill(row, col, sp);
  SECTION("Add method"){
    TestAdd(ma1,false);
  }
  SECTION("Add operator"){
    TestAdd(ma1,true);
  }
}

//for debugging: once this is fixed, uncomment next test
TEST_CASE("Add TPZSFMatrix","[matrix_tests][!shouldfail]"){
  SECTION("Find out why it is failing"){
    TPZSFMatrix<double> ma1;
    ma1.AutoFill(10,10,SymProp::Sym);
    TestAdd(ma1,false);
    TestAdd(ma1,true);
  }
}
TEMPLATE_PRODUCT_TEST_CASE("Add(sym)","[matrix_tests]",
                           (
                             // TPZSFMatrix,
                             TPZSBMatrix,
                             TPZSYsmpMatrix,TPZSkylMatrix),(float,double,long double,
                                                            std::complex<float>,std::complex<double>,
                                                            std::complex<long double>)){
  using MAT=TestType;
  const SymProp sp=GENERATE(SymProp::Sym,SymProp::Herm);
  constexpr int row{10};
  constexpr int col{10};
  SECTION(SymPropName(sp)){
    MAT ma1;
    ma1.AutoFill(row, col, sp);
    SECTION("Add method"){
      TestAdd(ma1,false);
    }
    SECTION("Add operator"){
      TestAdd(ma1,true);
    }
  }
}

//for debugging: once this is fixed, uncomment next test
TEST_CASE("Subtract TPZSkylNSymMatrix","[matrix_tests][!shouldfail]"){
  SECTION("With operator it executes, without it doesnt"){
    TPZSkylNSymMatrix<double> ma1;
    ma1.AutoFill(10,10,SymProp::NonSym);
    TestSubtract(ma1,false);
    TestSubtract(ma1,true);
  }
}

TEMPLATE_PRODUCT_TEST_CASE("Subtract(nsym)","[matrix_tests]",
                           (TPZFMatrix,
                            TPZFBMatrix,
                            TPZFYsmpMatrix
                            // ,TPZSkylNSymMatrix
                            ),(float,double,long double,
                               std::complex<float>,std::complex<double>,
                               std::complex<long double>)){
  using MAT=TestType;
  using SCAL=typename MAT::Type;
  const SymProp sp=SymProp::NonSym;
  constexpr int row{10};
  const int col = GENERATE(5,10,15);
  if(std::is_same_v<MAT,TPZFBMatrix<SCAL>> && col!=10){return;}
  MAT ma1;
  ma1.AutoFill(row, col, sp);
  SECTION("Subtract method"){
    TestSubtract(ma1,false);
  }
  SECTION("Subtract operator"){
    TestSubtract(ma1,true);
  }
}

//for debugging: once this is fixed, uncomment next test
TEST_CASE("Subtract TPZSFMatrix","[matrix_tests][!shouldfail]"){
  SECTION("Find out why it is failing"){
    TPZSFMatrix<double> ma1;
    ma1.AutoFill(10,10,SymProp::Sym);
    TestSubtract(ma1,false);
    TestSubtract(ma1,true);
  }
}

TEMPLATE_PRODUCT_TEST_CASE("Subtract(sym)","[matrix_tests]",
                           (
                             // TPZSFMatrix,
                             TPZSBMatrix,
                             TPZSYsmpMatrix,TPZSkylMatrix),(float,double,long double,
                                                            std::complex<float>,std::complex<double>,
                                                            std::complex<long double>)){
  using MAT=TestType;
  const SymProp sp=GENERATE(SymProp::Sym,SymProp::Herm);
  constexpr int row{10};
  constexpr int col{10};
  SECTION(SymPropName(sp)){
    MAT ma1;
    ma1.AutoFill(row, col, sp);
    SECTION("Subtract method"){
      TestSubtract(ma1,false);
    }
    SECTION("Subtract operator"){
      TestSubtract(ma1,true);
    }
  }
}

TEMPLATE_TEST_CASE("MultiplyByScalar","[matrix_tests][!shouldfail]",
                   // float,
                   // double,
                   // long double,
                   // std::complex<float>,
                   // std::complex<double>,
                   std::complex<long double>
                   ) {
  SECTION("IMPLEMENT ME"){
    FAIL("IMPLEMENT ME");
  }
}

TEMPLATE_TEST_CASE("Multiply","[matrix_tests][!shouldfail]",
                   // float,
                   // double,
                   // long double,
                   // std::complex<float>,
                   // std::complex<double>,
                   std::complex<long double>
                   ) {
  SECTION("IMPLEMENT ME"){
    FAIL("IMPLEMENT ME");
  }
}


//for debugging: once this is fixed, uncomment next test
TEST_CASE("MultAdd TPZFBMatrix","[matrix_tests][!shouldfail]"){
  SECTION("Find out why it is failing"){
    FAIL("no clue");
  }
}
TEST_CASE("MultAdd TPZSBMatrix","[matrix_tests][!shouldfail]"){
  SECTION("Find out why it is failing"){
    FAIL("no clue");
  }
}
TEST_CASE("MultAdd TPZSkylNSymMatrix","[matrix_tests][!shouldfail]"){
  SECTION("Find out why it is failing"){
    FAIL("no clue");
  }
}
TEST_CASE("MultAdd TPZSkylMatrix","[matrix_tests][!shouldfail]"){
  SECTION("Find out why it is failing"){
    FAIL("no clue");
  }
}


TEMPLATE_PRODUCT_TEST_CASE("MultAdd","[matrix_tests]",
                           (
                             TPZFMatrix,
                             // TPZFBMatrix,
                             // TPZSBMatrix,
                             // TPZSkylNSymMatrix,
                             // TPZSkylMatrix,
                             TPZFYsmpMatrix,
                             TPZSYsmpMatrix
#ifdef PZ_USING_MKL
                             ,TPZFYsmpMatrixPardiso,
                             TPZSYsmpMatrixPardiso
#endif
                            ),
                           (
                             float,
                             double,
                             long double,
                             std::complex<float>,
                             std::complex<double>,
                             std::complex<long double>)
                   ) {


  using MAT=TestType;
  using SCAL=typename MAT::Type;
  using RSCAL = RType(SCAL);

#ifdef PZ_USING_MKL
  //pardiso is only for doubles
  if(
    std::is_same_v<MAT,TPZFYsmpMatrixPardiso<float>> ||
    std::is_same_v<MAT,TPZFYsmpMatrixPardiso<long double>>||
    std::is_same_v<MAT,TPZFYsmpMatrixPardiso<std::complex<float>>>||
    std::is_same_v<MAT,TPZFYsmpMatrixPardiso<std::complex<long double>>>||
    std::is_same_v<MAT,TPZSYsmpMatrixPardiso<float>> ||
    std::is_same_v<MAT,TPZSYsmpMatrixPardiso<long double>>||
    std::is_same_v<MAT,TPZSYsmpMatrixPardiso<std::complex<float>>>||
    std::is_same_v<MAT,TPZSYsmpMatrixPardiso<std::complex<long double>>>
     ){return;}
#endif
  SCAL alpha{0.0};
  SCAL beta{0.0};
  TPZFMatrix<SCAL> x,y,z;
  MAT mat;
  const int nr{10};
  const int opt=GENERATE(0,1,2);
  SECTION("Zero sized A"){
    y.AutoFill(nr,1,SymProp::NonSym);
    alpha=0.0;
    beta=-1.0;
    x.Resize(0,1);
    if(opt){
      mat.Resize(0,nr);
    }else{
      mat.Resize(nr,0);
    }
    SECTION("z=beta*y"){
      mat.MultAdd(x,y,z,alpha,beta,opt);
      for(int i = 0; i < nr; i++){
        const auto yv = y.GetVal(i,0);
        const auto zv = z.GetVal(i,0);
        CAPTURE(i,yv,zv);
        REQUIRE(zv==-yv);
      }
    }
    SECTION("invalid x"){
      x.AutoFill(nr,1,SymProp::NonSym);
      REQUIRE_THROWS(mat.MultAdd(x,y,z,alpha,beta,opt));
    }
    SECTION("invalid y"){
      x.AutoFill(nr,1,SymProp::NonSym);
      mat.AutoFill(nr,nr,SymProp::Herm);
      y.AutoFill(2*nr,1,SymProp::NonSym);
      REQUIRE_THROWS(mat.MultAdd(x,y,z,alpha,beta,opt));
    }
  }
  SECTION("square case"){
    std::map<int,std::string> optname;
    optname[0]="NoTransp";
    optname[1]="Transp";
    optname[2]="TranspConj";
    SECTION(optname[opt]){
      const auto sp = GENERATE(SymProp::NonSym,SymProp::Herm);
      SECTION(SymPropName(sp)){
        if (std::is_same_v<MAT,TPZFMatrix<std::complex<long double>>>
            && opt==2){
          int var{2};
          var++;
        }
        if constexpr(IsSymmetricStorage<MAT>){
          if(sp==SymProp::NonSym){return;}
        }
        mat.AutoFill(nr,nr,sp);
        if constexpr (is_complex<SCAL>::value){
          alpha = (SCAL)1.0i;
        }else{
          alpha = 1.0;
        }
        beta = -alpha;
        TPZFMatrix<SCAL> aux = mat;
        if(opt==1){
          aux.Transpose();
        }else if(opt==2){
          aux.ConjTranspose();
        }

        DecomposeType dec{ENoDecompose};
        if constexpr(IsSymmetricStorage<MAT>){
          dec=ELDLt;
        }else{
          dec=ELU;
        }
      
        y.Resize(nr,nr);
        y.Identity();
        x=y;
        aux.SolveDirect(x,dec);
        mat.MultAdd(x,y,z,alpha,beta,opt);
        constexpr RSCAL tol = 1000*std::numeric_limits<RSCAL>::epsilon()/
          (std::numeric_limits<RSCAL>::digits10);
        const auto normz = Norm(z);
        CAPTURE(opt,normz,tol);
        REQUIRE(normz == Catch::Approx(0).margin(tol));
      }
    }
  }
}

template<class MAT>
void TestInverse(const bool sym_storage, const SymProp sp, const DecomposeType dec)
{
  using SCAL = typename MAT::Type;
  using RSCAL = RType(SCAL);
  
  MAT ma;
  const int dim = 5;
  ma.AutoFill(dim, dim, sp);

  CAPTURE(SymPropName(sp),DecomposeTypeName(dec),sym_storage);
  
  // Making ma copy because ma is modified by Inverse method (it's decomposed)
  MAT cpma(ma);
  TPZFMatrix<SCAL> inv(dim, dim), invkeep;
  TPZFMatrix<SCAL> res(inv);

  //it must fail in a few cases
  if(dec == ELU && sym_storage){
    REQUIRE_THROWS(ma.Inverse(inv, dec));
  }

  // getting inverse twice
  ma.Inverse(inv,dec);
  invkeep = inv;
  inv.Inverse(res, dec);
  bool check = true;
  /// Checking whether the res matrix is identical to m1 matrix
  auto oldPrecision = Catch::StringMaker<RSCAL>::precision;
  Catch::StringMaker<RSCAL>::precision = std::numeric_limits<RSCAL>::max_digits10;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      RSCAL diff =
        fabs(cpma.Get(i, j) - res.Get(i, j));
        
      bool loccheck = IsZero(diff / (RSCAL) 100.);
      if (loccheck == false) {
        CAPTURE(i, j, diff);
        std::cout << " i " << i << " j " << j << "diff " << diff << std::endl;
      }
      check &= loccheck;
    }
  }
  if(!check) {
    cpma.Print(" mat ",std::cout);
    inv.Print(" inv mat ", std::cout);
    res.Print(" inv inv mat ",std::cout);
  }
  REQUIRE(check);
  Catch::StringMaker<RSCAL>::precision = oldPrecision;
  
}


template<class MAT>
void TestAdd(const MAT& ma1, const bool use_operator){
  using TVar=typename MAT::Type;
  
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  
  //this is to ensure that they have the same sparsity pattern
  MAT ma2(ma1);
  MAT res;
  if(use_operator){
    res = ma1+ma2;
  }else{
    ma1.Add(ma2,res);
  }
    
  bool check = true;
  constexpr RTVar tol = [](){
    if constexpr (std::is_same_v<RTVar,float>) return (RTVar)10;
    else return (RTVar)1;
  }();
  const int row = ma1.Rows();
  const int col = ma1.Cols();
  
  for (int i = 0; i < row; i++) {
    for (int j = 0; (j < col) && check; j++) {
      CAPTURE(i,j,res.Get(i,j),(TVar)2*ma1.Get(i,j));
      if (!IsZero((res.Get(i, j)-(ma1.Get(i,j)+ma2.Get(i,j)))/tol)) {
        check = false;
      }
      REQUIRE(check);
    }
  }
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}

template<class MAT>
void TestSubtract(const MAT& ma1, const bool use_operator){
  using TVar=typename MAT::Type;

  
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  
  //this is to ensure that they have the same sparsity pattern
  MAT ma2(ma1);
  MAT res;
  if(use_operator){
    res = ma1-ma2;
  }else{
    ma1.Subtract(ma2,res);
  }
    
  bool check = true;
  constexpr RTVar tol = [](){
    if constexpr (std::is_same_v<RTVar,float>) return (RTVar)10;
    else return (RTVar)1;
  }();
  const int row = ma1.Rows();
  const int col = ma1.Cols();
  
  for (int i = 0; i < row; i++) {
    for (int j = 0; (j < col) && check; j++) {
      CAPTURE(i,j,res.Get(i,j));
      if (!IsZero((res.Get(i, j))/tol)) {
        check = false;
      }
      REQUIRE(check);
    }
  }
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}