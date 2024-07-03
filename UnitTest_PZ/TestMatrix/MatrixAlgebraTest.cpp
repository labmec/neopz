/**
 * @file
 * @brief Contains Unit Tests for algebraic operations using TPZMatrix<T> derived types
 */

#include "TestMatrixHeaders.h"

template <class MAT>
void TestInverse(MAT &mat, const DecomposeType dec);
template<class MAT>
void TestAdd(const MAT &mat, const bool use_operator);
template<class MAT>
void TestSubtract(const MAT &mat, const bool use_operator);
template<class MAT>
void TestMultiplyByScalar(const MAT &mat, const bool use_operator);
template<class TVar>
void BlockDiagLUPivot();
template<class TVar>
void BlockDiagZeroSizedBlock();
template<class TVar>
void SparseBlockDiagInverse();
template<class TVar>
void SparseBlockColorInverse();

TEMPLATE_PRODUCT_TEST_CASE("Inverse","[matrix_tests]",
                           (
                             TPZFMatrix,
                             TPZSFMatrix,
                             TPZFBMatrix,
                             TPZSBMatrix,
                             TPZSkylNSymMatrix,
                             TPZSkylMatrix,
#ifdef PZ_USING_MKL
                             TPZFYsmpMatrixPardiso,
                             TPZSYsmpMatrixPardiso,
#endif
                             TPZBlockDiagonal
                            ),
                           (float,
                            double,
                            long double,
                            std::complex<float>,
                            std::complex<double>,
                            std::complex<long double>)
                           ) {
  using MAT = TestType;
  using SCAL = typename MAT::Type;

  constexpr bool sym_storage=IsSymmetricStorage<MAT>;
#ifdef PZ_USING_MKL
  //pardiso is only for doubles
  if((std::is_same_v<MAT,TPZFYsmpMatrixPardiso<SCAL>> &&
      !std::is_same_v<double,RType(SCAL)>)||
     (std::is_same_v<MAT,TPZSYsmpMatrixPardiso<SCAL>> &&
      !std::is_same_v<double,RType(SCAL)>)){return;}
#endif
  DecomposeType dec = GENERATE(ELU,ECholesky,ELDLt);
  SECTION(DecomposeTypeName(dec)){
    SymProp sp = GENERATE(SymProp::NonSym,SymProp::Sym,SymProp::Herm);
    SECTION(SymPropName(sp)){
      if((dec==ECholesky && sp!=SymProp::Herm) ||
         (dec==ELDLt && sp!=SymProp::Herm)){
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
      if(sym_storage && dec==ELU){return;}
      MAT ma;
      const int dim = GENERATE(5,10);
      ma.AutoFill(dim, dim, sp);
      if (std::is_same_v<MAT, TPZSBMatrix<SCAL>> && is_complex<SCAL>::value && dec==ELDLt)
      {
        // ELDLt is not implemented for complex numbers in TPZSBMatrix. Please Implement it!
        REQUIRE_THROWS(TestInverse(ma, dec));
      }
      else
      {
        TestInverse(ma, dec);
      }
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
  SECTION("It fails with and without operator"){
    TPZSkylNSymMatrix<double> ma1;
    ma1.AutoFill(10,10,SymProp::NonSym);
    TestAdd(ma1,true);
    TestAdd(ma1,false);
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

TEMPLATE_PRODUCT_TEST_CASE("Add(sym)","[matrix_tests]",
                           (
                             TPZSFMatrix,
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
  SECTION("It fails with and without operator"){
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

TEMPLATE_PRODUCT_TEST_CASE("Subtract(sym)","[matrix_tests]",
                           (
                             TPZSFMatrix,
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

//for debugging: once this is fixed, uncomment next test
TEST_CASE("MultiplyByScalar TPZSkylNSymMatrix","[matrix_tests][!shouldfail]"){
  SECTION("It fails with and without operator"){
    TPZSkylNSymMatrix<double> ma1;
    ma1.AutoFill(10,10,SymProp::NonSym);
    TestMultiplyByScalar(ma1,false);
    TestMultiplyByScalar(ma1,true);
  }
}

TEMPLATE_PRODUCT_TEST_CASE("MultiplyByScalar", "[matrix_tests][!shouldfail]",
                           (TPZFMatrix,
                            TPZFBMatrix,
                            TPZFYsmpMatrix,
                            // TPZSkylNSymMatrix,
                            TPZSFMatrix,
                            TPZSBMatrix,
                            TPZSYsmpMatrix,
#ifdef PZ_USING_MKL
                            TPZFYsmpMatrixPardiso,
                            TPZSYsmpMatrixPardiso,
#endif
                            TPZSkylMatrix),
                           (float, double, long double,
                            std::complex<float>, std::complex<double>,
                            std::complex<long double>))
{
  using MAT = TestType;
  using SCAL = typename MAT::Type;

  constexpr bool sym_storage = IsSymmetricStorage<MAT>;
#ifdef PZ_USING_MKL
  // pardiso is only for doubles
  if ((std::is_same_v<MAT, TPZFYsmpMatrixPardiso<SCAL>> &&
       !std::is_same_v<double, RType(SCAL)>) ||
      (std::is_same_v<MAT, TPZSYsmpMatrixPardiso<SCAL>> &&
       !std::is_same_v<double, RType(SCAL)>))
  {
    return;
  }
#endif
  SymProp sp = GENERATE(SymProp::NonSym, SymProp::Sym, SymProp::Herm);
  SECTION(SymPropName(sp))
  {
    if (sym_storage && sp == SymProp::NonSym)
    {
      return;
    }

    MAT ma;
    const int dim = GENERATE(5, 10);
    ma.AutoFill(dim, dim, sp);

    SECTION("MultiplyByScalar method")
    {
      TestMultiplyByScalar(ma, false);
    }
    SECTION("MultiplyByScalar operator")
    {
      TestMultiplyByScalar(ma, true);
    }
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
  SECTION("Multiply method IMPLEMENT ME"){
    FAIL("IMPLEMENT ME");
  }
  SECTION("Multiply operator IMPLEMENT ME"){
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
  if((std::is_same_v<MAT,TPZFYsmpMatrixPardiso<SCAL>> &&
      !std::is_same_v<double,RType(SCAL)>)||
     (std::is_same_v<MAT,TPZSYsmpMatrixPardiso<SCAL>> &&
      !std::is_same_v<double,RType(SCAL)>)){return;}
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
      const auto sp = GENERATE(SymProp::NonSym,SymProp::Sym,SymProp::Herm);
      SECTION(SymPropName(sp)){
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
        aux.SolveDirect(x,ELU);
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

TEMPLATE_TEST_CASE("Additional Block Diagonal tests","[matrix_tests]",
                   float,
                   double,
                   // long double,
                   std::complex<float>,
                   std::complex<double>
                   // std::complex<long double>
                   ){
  SECTION("LU Pivot"){
    BlockDiagLUPivot<TestType>();
  }
  SECTION("BlockDiag block with zero size"){
    BlockDiagZeroSizedBlock<TestType>();
  }
  SECTION("Sparse Block Inverse"){
    SparseBlockDiagInverse<TestType>();
    SparseBlockColorInverse<TestType>();
  }
}


#include "pzseqsolver.h"
#include "pzstepsolver.h"

template<class MAT>
void TestInverse(MAT &ma, const DecomposeType dec)
{
  using SCAL = typename MAT::Type;
  using RSCAL = RType(SCAL);
  const auto dim = ma.Rows();
  REQUIRE(dim==ma.Cols());
  // Making ma copy because ma is modified by Inverse method (it's decomposed)
  MAT cpma(ma);
  TPZFMatrix<SCAL> inv(dim, dim), invkeep;
  TPZFMatrix<SCAL> res(inv);
  // getting inverse twice
  ma.Inverse(inv,dec);
  invkeep = inv;
  //since it is a full matrix, we use LU to be sure
  inv.Inverse(res, ELU);
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
    invkeep.Print(" inv mat ", std::cout);
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

template<class MAT>
void TestMultiplyByScalar(const MAT &ma, const bool use_operator)
{
  using TVar=typename MAT::Type;
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;

  bool check = true;
  constexpr RTVar tol = []()
  {
    if constexpr (std::is_same_v<RTVar, float>)
      return (RTVar)10;
    else
      return (RTVar)1;
  }();
  const int row = ma.Rows();
  const int col = ma.Cols();

  RTVar alpha = GENERATE(0.1, 2.0, 10.0);
  MAT res;
  if (use_operator)
  {
    res = ma * alpha;
  }
  else
  {
    ma.MultiplyByScalar(alpha,res);
  }

  for (int i = 0; i < row; i++)
  {
    for (int j = 0; (j < col) && check; j++)
    {
      CAPTURE(i, j, res.GetVal(i, j), ma.GetVal(i, j), alpha);
      if (!IsZero((res.Get(i, j) - ma.Get(i, j) * alpha) / tol))
      {
        check = false;
      }
      REQUIRE(check);
    }
  }
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}

template<class TVar>
void BlockDiagLUPivot(){
  //we will create a matrix with 2x2 blocks that will require pivoting
  const int nblocks = 5;
  const int dim = 2*nblocks;
  TPZFMatrix<TVar> fmat(dim,dim,0);
  for(int i = 0; i < nblocks; i++){
    fmat.PutVal(2*i, 2*i, 0);
    fmat.PutVal(2*i, 2*i+1, 1);
    fmat.PutVal(2*i+1, 2*i, 1);
    fmat.PutVal(2*i+1, 2*i+1, 0);
  }
  TPZVec<int> blocksizes(nblocks,2);
  TPZBlockDiagonal<TVar> blckmat(blocksizes,fmat);
  TestInverse(blckmat, ELU);
}

template<class TVar>
void BlockDiagZeroSizedBlock(){
  //we will create a matrix with 2x2 blocks that will require pivoting
  const int nblocks = 5;
  const int dim = 2*nblocks;
  TPZFMatrix<TVar> fmat(dim,dim,0);
  for(int i = 0; i < nblocks; i++){
    fmat.PutVal(2*i, 2*i, 0);
    fmat.PutVal(2*i, 2*i+1, 1);
    fmat.PutVal(2*i+1, 2*i, 1);
    fmat.PutVal(2*i+1, 2*i+1, 0);
  }
  //create additional block with size 0 at the end of block list
  TPZVec<int> blocksizes(nblocks+1,2);
  blocksizes[nblocks] = 0;
  REQUIRE_NOTHROW(TPZBlockDiagonal<TVar>(blocksizes,fmat));
}

template<class TVar>
void SparseBlockDiagInverse(){

  constexpr int neq = 10;
  constexpr int nblocks = 3;
  constexpr int bsize = 2;

  TPZFMatrix<TVar> fmat(neq,neq);
  fmat.AutoFill(neq, neq, SymProp::NonSym);
  //now we will extract a few blocks from the full matrix
  //non-null equations
  TPZVec<int64_t> blockgraph =
    {
      0,1,//first block
      4,5,6,//second block
      8,9//third block
    };

  //index in blockgraph of first equation of each block
  TPZVec<int64_t> blockgraphindex =
    {
      0,//first block
      2,//second block
      5,//third block
      7//size of blockgraph
    };

  TPZSparseBlockDiagonal<TVar> blmat(blockgraph,blockgraphindex,neq),
    blcp(blockgraph,blockgraphindex, neq);
  //now we remove the coupling from the full matrix
  for (auto bleq : blockgraph){
    int64_t block{0}, blockind{0};
    blmat.FindBlockIndex(bleq, block, blockind);
    for(int ieq : blockgraph){
      int64_t blockother{0};
      blmat.FindBlockIndex(ieq, blockother, blockind);
      if(blockother != -1 && blockother != block){
        fmat.PutVal(bleq, ieq, 0);
        fmat.PutVal(ieq, bleq, 0);
      }
    }
  }
  blmat.BuildFromMatrix(fmat);
  blcp.BuildFromMatrix(fmat);

  //rhs, initial solution, residual and sol update
  TPZFMatrix<TVar> rhs(neq,1,0), u0(neq,1,0), res(neq,1,0),
    du(neq,1,0), resblck(neq,1,0);
  rhs.AutoFill(neq,1,SymProp::NonSym);

  fmat.Residual(u0, rhs, res);

  // fmat.Print("full",std::cout,EMathematicaInput);
  // blmat.Print("blck",std::cout,EMathematicaInput);
  // blmat.Print("blck",std::cout,EFormatted);
  
  du = res;
  blmat.Solve_LU(&du);

  
  constexpr RTVar tol =std::numeric_limits<RTVar>::epsilon()*2000;
  //now we check if the residual is null for every eq in a block
  SECTION("Residual (full mat)"){
    fmat.Residual(du,rhs,res);
    for(auto ieq : blockgraph){
      CAPTURE(ieq);
      CAPTURE(res.GetVal(ieq,0));
      REQUIRE((std::abs(res.GetVal(ieq,0)) == Catch::Approx(0).margin(tol)));
    }
  }

  SECTION("Residual (block mat)"){
    blcp.Residual(du,rhs,resblck);
    for(auto ieq : blockgraph){
      CAPTURE(ieq);
      CAPTURE(resblck.GetVal(ieq,0));
      REQUIRE((std::abs(resblck.GetVal(ieq,0)) == Catch::Approx(0).margin(tol)));
    }
  }
}

template<class TVar>
void SparseBlockColorInverse(){

  constexpr int neq = 10;
  constexpr int nblocks = 4;
  //expected solution is {1,1,1,1,1,1,1,1,1,1}^T
  TPZAutoPointer<TPZFMatrix<TVar>> fmat =
    new TPZFMatrix<TVar>({
        {1, 2, 0, 1, 0, 0, 0, 0, 0, 0},
        {3, 1, 0, 2, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 4, 1, 0, 0, 0, 0},
        {3, 2, 0, 3, 0, 0, 0, 0, 0, 0},
        {0, 0, 2, 0, 1, 1, 0, 0, 0, 0},
        {0, 0, 3, 0, 5, 4, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 2, 0},
        {0, 0, 0, 0, 0, 0, 0, 6, 0, 2},
        {0, 0, 0, 0, 0, 0, 3, 0, 4, 0},
        {0, 0, 0, 0, 0, 0, 0, 2, 0, 5}
    });

  // TPZAutoPointer<TPZFYsmpMatrix<TVar>> smat = new TPZFYsmpMatrix<TVar>(neq,neq);
  // TPZVec<int64_t> ia = {0,3,6,9,12,15,18,20,22,24,26};
  // TPZVec<int64_t> ja = {0,1,3,0,1,3,2,4,5,0,1,3,2,4,5,2,4,5,6,8,7,9,6,8,7,9};
  // TPZVec<TVar> aa =    {1,2,1,3,1,2,1,4,1,3,2,3,2,1,1,3,5,4,1,2,6,2,3,4,2,5};
  // smat->SetData(ia,ja,aa);
  
  TPZFMatrix<TVar> rhs = {4,6,6,8,4,12,3,8,7,7};
  //now we will extract a few blocks from the full matrix
  //non-null equations
  TPZVec<int64_t> blockgraph =
    {
      0,1,3,//first block
      2,4,5,//second block
      6,8,//third block
      7,9//fourth block
    };

  //index in blockgraph of first equation of each block
  TPZVec<int64_t> blockgraphindex =
    {
      0,//first block
      3,//second block
      6,//third block
      8,//fourth block
      10//size of blockgraph
    };
  constexpr int numcolors{2};
  //block colors
  TPZVec<int> colors =
    {
      0,//first block
      1,//second block
      0,//third block
      1//fourth bock
    };

  TPZSequenceSolver<TVar> seqsolv;
  for(int c = 0; c < numcolors; c++){
    
    TPZAutoPointer<TPZSparseBlockDiagonal<TVar>> blmat
      = new TPZSparseBlockDiagonal<TVar>(blockgraph,blockgraphindex,neq,c,colors);
    blmat->BuildFromMatrix(*fmat);
    
    TPZStepSolver<TVar> step(blmat);
    step.SetDirect(ELU);
    seqsolv.AppendSolver(step);
  }
  seqsolv.SetMatrix(fmat);
  TPZFMatrix<TVar> res = rhs;



  constexpr RTVar tol =std::numeric_limits<RTVar>::epsilon()*2000;
  /*
  //just to ensure that the solution is correct
  TPZAutoPointer<TPZMatrix<TVar>> fmatcp = fmat->Clone();
  fmatcp->SolveDirect(res,ELU);
  for(int i = 0; i < neq; i++){
    REQUIRE(std::abs(res.GetVal(i,0)-(TVar)1.)==Catch::Approx(0).margin(tol));
  }
  */
  seqsolv.Solve(rhs, res);
  for(int i = 0; i < neq; i++){
    REQUIRE(std::abs(res.GetVal(i,0)-(TVar)1.)==Catch::Approx(0).margin(tol));
  }
}
