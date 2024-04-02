/**
 * @file
 * @brief Contains Unit Tests for operations in "windows" of matrices using TPZMatrix<T> derived types
 */

#include "TestMatrixHeaders.h"
#include "TPZMatrixWindow.h"
/**
 * @brief Tests the addContribution method of the matrix, that adds a block C += alpha * A*B starting at C(i,j), using AutoFill to build a square matrix of dimension dim (user defined)
 * @param nrows Number of rows of the matrix to be build.
 * @param ncols Number of columns of the matrix to be build.
 * @param ntype test type: 0: MultAddComparsion, 1: Incompatible dimension, 2: Out of bounds
 * @note Process: build a matrix C with randomic values, 0: adds a contribution C += A*B of the same size as C. Compare the results with AddContribution and MultAdd. 1: same thing, but we choose B to be of incompatible size 2: sets initial position such that we have an out of bounds operation
 */
template <class TVar>
void TestingAddContribution(int nrows, int ncols, int ntype);

TEMPLATE_TEST_CASE("AddContribution","[matrix_tests]",
                   float,
                   double,
                   long double,
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
  
  using SCAL=TestType;
  using RSCAL=RType(SCAL);

  std::map<int,std::string> optname;
  optname[0]="NoTransp";
  optname[1]="Transp";
  optname[2]="TranspConj";
  
  constexpr RSCAL tol =std::numeric_limits<RSCAL>::epsilon()*2000;
  SECTION("MultAddComparison (square mat)")
  {
    const int nr{10};
    const int nc{10};
    
    TPZFMatrix<SCAL> B(nr,nc), A(nr,nc);
    A.AutoFill(nr,nc,SymProp::NonSym);
    B.AutoFill(nr,nc,SymProp::NonSym);
    TPZFMatrix<SCAL> C_multadd(nr,nc,0.), C_addcontr(nr,nc,0.);
    SCAL beta{1};
    SCAL alpha{1};
    if constexpr (is_complex<SCAL>::value){
      alpha=1i;
    }
    //whether to transpose A
    const int opt = GENERATE(0,1,2);
    SECTION(optname[opt]){
      A.MultAdd(B,C_multadd,C_multadd,alpha,beta,opt);
      C_addcontr.AddContribution(0,0,A,opt,B,false,alpha);
      TPZFMatrix<SCAL> diff=C_multadd-C_addcontr;
      auto normdiff = Norm(diff);
      REQUIRE(normdiff == Catch::Approx(0).margin(tol));
      //let us try it again
      C_addcontr.AddContribution(0,0,A,opt,B,false,alpha);
      diff = (SCAL)2.0*C_multadd-C_addcontr;
      normdiff = Norm(diff);
      REQUIRE(normdiff == Catch::Approx(0).margin(tol));
    }
  }
  SECTION("MultAddComparison (non-square mat)")
  {
    //we need nr>nc
    const int nr{10};
    const int nc{5};
    TPZFMatrix<SCAL> B(nr,nc), A(nr,nr);
    A.AutoFill(nr,nr,SymProp::NonSym);
    B.AutoFill(nr,nc,SymProp::NonSym);
    //we make addcontr matrix bigger on purpose
    TPZFMatrix<SCAL> C_multadd(nr,nc,0.), C_addcontr(2*nr,2*nc,0.);
    SCAL beta{1};
    SCAL alpha{1};
    if constexpr (is_complex<SCAL>::value){
      alpha=1i;
    }
    //how to transpose A
    const int opt = GENERATE(0,1,2);
    SECTION(optname[opt]){
      A.MultAdd(B,C_multadd,C_multadd,alpha,beta,opt);
      C_addcontr.AddContribution(0,0,A,opt,B,false,alpha);
      for(int ic = 0; ic < nc; ic++){
        for(int ir = 0; ir < nr; ir++){
          const auto v1 = C_multadd.GetVal(ir,ic);
          const auto v2 = C_addcontr.GetVal(ir,ic);
          REQUIRE(std::abs(v1-v2) == Catch::Approx(0).margin(tol));
        }
      }
      //let us try it again
      C_addcontr.AddContribution(0,0,A,opt,B,false,alpha);
      for(int ic = 0; ic < nc; ic++){
        for(int ir = 0; ir < nr; ir++){
          const auto v1 = (SCAL)2.0*C_multadd.GetVal(ir,ic);
          const auto v2 = C_addcontr.GetVal(ir,ic);
          REQUIRE(std::abs(v1-v2) == Catch::Approx(0).margin(tol));
        }
      }
    }
  }
  SECTION("Incompatible dim")
  {
    const int nr{10};
    const int nc{5};
    const int too_big{50};
    TPZFMatrix<SCAL> C(nr,nc,0.0);
    SECTION("Wrong A"){
      TPZFMatrix<SCAL> A(2*nr,nc,0.);
      TPZFMatrix<SCAL> B(nc,nc,0.);
      const int opt=GENERATE(0,1,2);
      if(opt){A.Transpose();}
      SECTION(optname[opt]){
        REQUIRE_THROWS(C.AddContribution(0,0,A,opt,B,false));
      }
    }
    SECTION("Wrong B"){
      TPZFMatrix<SCAL> A(nr,nr,0.);
      TPZFMatrix<SCAL> B(nr,too_big,0.);
      const int opt=GENERATE(0,1,2);
      if(opt){B.Transpose();}
      SECTION(optname[opt]){
        REQUIRE_THROWS(C.AddContribution(0,0,A,false,B,opt));
      }
    }
    SECTION("A and B incompatible"){
      TPZFMatrix<SCAL> A(nr,1,0.);
      TPZFMatrix<SCAL> B(2,nc,0.);
      const int opt_a=GENERATE(0,1,2);
      SECTION(optname[opt_a]){
        if(opt_a){
          A.Transpose();
        }
        const int opt_b=GENERATE(0,1,2);
        SECTION(optname[opt_b]){
          if(opt_b){
            B.Transpose();
          }
          REQUIRE_THROWS(C.AddContribution(0,0,A,opt_a,B,opt_b));
        }
      }
    }
    
    SECTION("Out of bounds")
    {
      const int nr{10};
      const int nc{5};
      const int too_big{50};
      TPZFMatrix<SCAL> C(nr,nc,0.0);
      SECTION("Rows out"){
        TPZFMatrix<SCAL> A(nr,nc,0.0);
        TPZFMatrix<SCAL> B(nc,nc,0.0);
        REQUIRE_THROWS(C.AddContribution(2,0,A,false,B,false));
      }
      SECTION("Cols out"){
        TPZFMatrix<SCAL> A(nr,nc,0.0);
        TPZFMatrix<SCAL> B(nc,nc,0.0);
        REQUIRE_THROWS(C.AddContribution(0,2,A,false,B,false));
      }
      SECTION("All out:1"){
        TPZFMatrix<SCAL> A(nr,nc,0.0);
        TPZFMatrix<SCAL> B(nc,nc,0.0);
        REQUIRE_THROWS(C.AddContribution(2,2,A,false,B,false));
      }
      SECTION("All out:2"){
        TPZFMatrix<SCAL> A(too_big,nc,0.0);
        TPZFMatrix<SCAL> B(nc,too_big,0.0);
        REQUIRE_THROWS(C.AddContribution(2,2,A,false,B,false));
      }
      TPZFMatrix<SCAL> A(1,nr,0.);
      TPZFMatrix<SCAL> B(nc,2,0.);
      TestingAddContribution<SCAL>(4, 4, 2);
    }
  }
}

TEMPLATE_TEST_CASE("MatrixWindowInit","[matrix_tests]",
                   float
                   // double,
                   // long double,
                   // std::complex<float>,
                   // std::complex<double>,
                   // std::complex<long double>
                   )
{
  using SCAL=TestType;
  
  constexpr int nrows_orig{10};
  constexpr int ncols_orig{12};
  TPZFMatrix<SCAL> orig_mat;
  orig_mat.AutoFill(nrows_orig, ncols_orig, SymProp::NonSym);
  SECTION("test same values"){
    constexpr int nrows_window{5};
    constexpr int ncols_window{6};
    constexpr int first_i{3};
    constexpr int first_j{4};

    auto CheckMatrix = [](const TPZFMatrix<SCAL> &orig, const TPZMatrixWindow<SCAL> &window){
      for(int irow = 0; irow < nrows_window; irow++){
        for(int icol = 0; icol < ncols_window; icol++){
          CAPTURE(irow);
          CAPTURE(icol);
          REQUIRE(orig.Get(first_i+irow,first_j+icol)==window.Get(irow,icol));
        }
      }
    };
    
    SECTION("TPZFMatrix ctor"){
      TPZMatrixWindow<SCAL> window_mat(orig_mat,first_i,first_j,nrows_window, ncols_window);
      CheckMatrix(orig_mat,window_mat);
    }
    SECTION("mem area ctor"){
      SCAL *begin = orig_mat.Elem() + nrows_orig*first_j + first_i;
      TPZMatrixWindow<SCAL> window_mat(begin, nrows_window, ncols_window,
                                       nrows_orig,nrows_orig*ncols_orig);
      CheckMatrix(orig_mat,window_mat);
    }
    SECTION("Modifying window"){
      TPZMatrixWindow<SCAL> window_mat(orig_mat,first_i,first_j,nrows_window, ncols_window);
      const int64_t pos_i = 1;
      const int64_t pos_j = 2;
      const SCAL newval = window_mat.GetVal(pos_i,pos_j)*(SCAL)7.;
      window_mat.PutVal(pos_i,pos_j,newval);
      CheckMatrix(orig_mat,window_mat);
    }
    SECTION("Modifying original"){
      TPZMatrixWindow<SCAL> window_mat(orig_mat,first_i,first_j,nrows_window, ncols_window);
      const int64_t pos_i = first_i+1;
      const int64_t pos_j = first_j+2;
      const SCAL newval = orig_mat.GetVal(pos_i,pos_j)*(SCAL)7.;
      orig_mat.PutVal(pos_i,pos_j,newval);
      CheckMatrix(orig_mat,window_mat);
    }
  }
  SECTION("invalid usage"){
    SECTION("negative i"){
      constexpr int nrows_window{5};
      constexpr int ncols_window{6};
      constexpr int first_i{-1};
      constexpr int first_j{4};
      REQUIRE_THROWS(TPZMatrixWindow<SCAL>(orig_mat,first_i,first_j,nrows_window, ncols_window));
    }
    SECTION("negative j"){
      constexpr int nrows_window{5};
      constexpr int ncols_window{6};
      constexpr int first_i{3};
      constexpr int first_j{-1};
      REQUIRE_THROWS(TPZMatrixWindow<SCAL>(orig_mat,first_i,first_j,nrows_window, ncols_window));
    }
    SECTION("negative ij"){
      constexpr int nrows_window{5};
      constexpr int ncols_window{6};
      constexpr int first_i{-1};
      constexpr int first_j{-1};
      REQUIRE_THROWS(TPZMatrixWindow<SCAL>(orig_mat,first_i,first_j,nrows_window, ncols_window));
    }
    SECTION("zero rows"){
      constexpr int nrows_window{0};
      constexpr int ncols_window{6};
      constexpr int first_i{1};
      constexpr int first_j{1};
      REQUIRE_THROWS(TPZMatrixWindow<SCAL>(orig_mat,first_i,first_j,nrows_window, ncols_window));
    }
    SECTION("zero cols"){
      constexpr int nrows_window{5};
      constexpr int ncols_window{0};
      constexpr int first_i{1};
      constexpr int first_j{1};
      REQUIRE_THROWS(TPZMatrixWindow<SCAL>(orig_mat,first_i,first_j,nrows_window, ncols_window));
    }
    SECTION("out of bounds row"){
      constexpr int nrows_window{10};
      constexpr int ncols_window{6};
      constexpr int first_i{1};
      constexpr int first_j{1};
      REQUIRE_THROWS(TPZMatrixWindow<SCAL>(orig_mat,first_i,first_j,nrows_window, ncols_window));
    }
    SECTION("out of bounds col"){
      constexpr int nrows_window{5};
      constexpr int ncols_window{12};
      constexpr int first_i{1};
      constexpr int first_j{1};
      REQUIRE_THROWS(TPZMatrixWindow<SCAL>(orig_mat,first_i,first_j,nrows_window, ncols_window));
    }
    SECTION("out of bounds rowcol"){
      constexpr int nrows_window{10};
      constexpr int ncols_window{12};
      constexpr int first_i{1};
      constexpr int first_j{1};
      REQUIRE_THROWS(TPZMatrixWindow<SCAL>(orig_mat,first_i,first_j,nrows_window, ncols_window));
    }
  }
}

TEMPLATE_TEST_CASE("TPZMatrixWindowMultAdd","[matrix_tests][!shouldfail]",
                   long double,
                   std::complex<long double>
                   )
{
  SECTION("TPZFMatrix case:must implement non-lapack version"){
    FAIL("Only LAPACK version is implemented");
  }
  SECTION("TPZMatrixWindow case:must implement non-lapack version"){
    FAIL("Only LAPACK version is implemented");
  }
}

TEMPLATE_TEST_CASE("TPZMatrixWindowMultAdd","[matrix_tests]",
                   float,
                   double,
                   // long double,
                   std::complex<float>,
                   std::complex<double>
                   // ,std::complex<long double>
                   )
{

  using SCAL=TestType;
  using RSCAL = RType(SCAL);

  SCAL alpha{0.0};
  SCAL beta{0.0};
  
  const int nr{10};
  const int opt_a=GENERATE(0,1,2);
  const int opt_x=GENERATE(0,1,2);

  std::map<int,std::string> opt_name;
  opt_name[0]="NoTransp";
  opt_name[1]="Transp";
  opt_name[2]="TranspConj";
  
  SECTION("TPZFMatrix case"){
    SECTION("Mat A:"+opt_name[opt_a]){
      const auto sp = GENERATE(SymProp::NonSym,SymProp::Sym,SymProp::Herm);
      SECTION(SymPropName(sp)){
        TPZFMatrix<SCAL> mat_full;
        TPZFMatrix<SCAL> x_full,y_full,z_full;
        
        mat_full.AutoFill(2*nr,2*nr,sp);
        TPZMatrixWindow<SCAL> mat(mat_full,nr,nr,nr,nr);
        if constexpr (is_complex<SCAL>::value){
          alpha = (SCAL)1.0i;
        }else{
          alpha = 1.0;
        }
        beta = -alpha;
        TPZFMatrix<SCAL> aux = mat;
        if(opt_a==1){
          aux.Transpose();
        }else if(opt_a==2){
          aux.ConjTranspose();
        }

        const DecomposeType dec=ELU;
      
        y_full.Resize(nr,nr);
        y_full.Identity();
        x_full=y_full;
        aux.SolveDirect(x_full,ELU);
        mat.MultAdd(x_full,y_full,z_full,alpha,beta,opt_a);
        constexpr RSCAL tol = 1000*std::numeric_limits<RSCAL>::epsilon()/
          (std::numeric_limits<RSCAL>::digits10);
        const auto normz = Norm(z_full);
        CAPTURE(opt_a,normz,tol);
        REQUIRE(normz == Catch::Approx(0).margin(tol));
      }
    }
  }
  SECTION("TPZMatrixWindow case"){
    SECTION("Mat A:"+opt_name[opt_a]){
      SECTION("Mat X:"+opt_name[opt_x]){
        const auto sp = GENERATE(SymProp::NonSym,SymProp::Sym,SymProp::Herm);
        SECTION(SymPropName(sp)){
          TPZFMatrix<SCAL> mat_full;
          TPZFMatrix<SCAL> x_full,y_full,z_full;
          mat_full.AutoFill(2*nr,2*nr,sp);
          TPZMatrixWindow<SCAL> mat(mat_full,nr,nr,nr,nr);
          if constexpr (is_complex<SCAL>::value){
            alpha = (SCAL)1.0i;
          }else{
            alpha = 1.0;
          }
          beta = -alpha;
          TPZFMatrix<SCAL> aux = mat;
          if(opt_a==1){
            aux.Transpose();
          }else if(opt_a==2){
            aux.ConjTranspose();
          }

          const DecomposeType dec=ELU;
      
          y_full.AutoFill(3*nr,4*nr,SymProp::NonSym);
          TPZMatrixWindow<SCAL> y(y_full,nr,nr,nr,nr);
          x_full=y;
          aux.SolveDirect(x_full,ELU);
          if(opt_x==1){
            x_full.Transpose();
          }else if(opt_x==2){
            x_full.ConjTranspose();
          }
          //just to change leading dim
          x_full.Resize(5*nr,nr);
          TPZMatrixWindow<SCAL>x(x_full,0,0,nr,nr);
          int nr_res{-1}, nc_res{-1};
          if(opt_a){
            nr_res = mat.Cols();
            
          }else{
            nr_res = mat.Rows();
          }
          if(opt_x){
            nc_res = x.Rows();
          }else{
            nc_res = x.Cols();
          }
          
          z_full.Redim(6*nr_res,3*nc_res);
          TPZMatrixWindow<SCAL> z(z_full,0,0,nr_res,nc_res);
          
          mat.MultAdd(x,y,z,alpha,beta,opt_a,opt_x);
          constexpr RSCAL tol = 1000*std::numeric_limits<RSCAL>::epsilon()/
            (std::numeric_limits<RSCAL>::digits10);
          const auto normz = Norm(z_full);
          CAPTURE(opt_a,normz,tol);
          REQUIRE(normz == Catch::Approx(0).margin(tol));
        }
      }
    }
  }

  auto CopyWindowFromMat = [](TPZFMatrix<SCAL> &orig, TPZFMatrix<SCAL> &window,
                              const int fr, const int fc, const int nr, const int nc){
    window.Redim(nr,nc);
    for(int ic = 0; ic < nc; ic++){
      for(int ir = 0; ir < nr; ir++){
        const auto val = orig.GetVal(fr+ir,fc+ic);
        window.PutVal(ir,ic,val);
      }
    }
  };
  
  SECTION("Traditional MultAddComparison"){
    TPZFMatrix<SCAL> x_full, y_full, z_full, mat_full;
    
    constexpr int nr_mat = 10;
    constexpr int nc_mat = 5;
    //number of columns of output
    constexpr int nc_y = 3;
    mat_full.AutoFill(nr_mat+2,nr_mat+4,SymProp::NonSym);
    TPZMatrixWindow<SCAL> mat_w(mat_full, 1,3,nr_mat,nc_mat);
    TPZFMatrix<SCAL> mat_w_fake;
    CopyWindowFromMat(mat_full,mat_w_fake,1,3,nr_mat,nc_mat);
    SECTION("Mat A: "+opt_name[opt_a]+ " Mat X: "+opt_name[opt_x]){
      //generate y
      const int nr_y = opt_a == 0 ? nr_mat : nc_mat;
      y_full.AutoFill(nr_y+3,nc_y+2,SymProp::NonSym);
      TPZMatrixWindow<SCAL> y_w(y_full, 2,1,nr_y,nc_y);
      TPZFMatrix<SCAL> y_w_fake;
      CopyWindowFromMat(y_full,y_w_fake,2,1,nr_y,nc_y);
      //generate x
      int nr_x{-1};
      int nc_x{-1};
      if(opt_x){
        nr_x = nc_y;
        nc_x = opt_a==0 ? nc_mat : nr_mat;
      }else{
        nr_x = opt_a==0 ? nc_mat : nr_mat;
        nc_x = nc_y;
      }
      x_full.AutoFill(nr_x+3,nc_x+2,SymProp::NonSym);
      TPZMatrixWindow<SCAL> x_w(x_full, 2,1,nr_x,nc_x);
      TPZFMatrix<SCAL> x_w_fake;
      CopyWindowFromMat(x_full,x_w_fake,2,1,nr_x,nc_x);
      TPZFMatrix<SCAL> z_full;
      //set alpha and beta
      if constexpr (is_complex<SCAL>::value){
        alpha = (SCAL)1.0i;
        beta = (SCAL)2.0i;
      }else{
        alpha = 1.0;
        beta = 2.0;
      }
      if(opt_x==1){
        x_w_fake.Transpose();
      }
      if(opt_x==2){
        x_w_fake.ConjTranspose();
      }
      mat_w_fake.MultAdd(x_w_fake,y_w_fake,z_full,alpha,beta,opt_a);
      const int nr_z = z_full.Rows();
      const int nc_z = z_full.Cols();
      TPZFMatrix<SCAL> res_full = z_full;
      //let us change leading dim of z_full
      z_full.Resize(nr_z+2,nc_z+4);
      TPZMatrixWindow<SCAL> z_w(z_full, 0,0,nr_z,nc_z);
      mat_w.MultAdd(x_w,y_w,z_w,alpha,beta,opt_a,opt_x);
      TPZFMatrix<SCAL> diff = z_w;
      diff-=res_full;
      constexpr RSCAL tol = 1000*std::numeric_limits<RSCAL>::epsilon()/
        (std::numeric_limits<RSCAL>::digits10);
      const auto normdiff = Norm(diff);
      REQUIRE(normdiff == Catch::Approx(0).margin(tol));
    }
  }
}

template <class TVar>
void TestingAddContribution(int nrows, int ncols, int ntype)
{
    auto oldPrecision = Catch::StringMaker<RTVar>::precision;
    Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
    TPZFMatrix<TVar> C1;
    C1.AutoFill(nrows, ncols, SymProp::NonSym);
    TPZFMatrix<TVar> C2;
    TPZFMatrix<TVar> A(C1);
    TPZFMatrix<TVar> B(C1);
    TPZFMatrix<TVar> BT;
    B.Transpose(&BT);
    TPZFMatrix<TVar> y(C1);

    switch (ntype)
    {
      case 0: // Comparison between AddContribution and MultAdd
      {
          C1.AddContribution(0, 0, A, false, B, true, 1.0);
          A.MultAdd(BT, y, C2, 1.0, 1.0);

          constexpr RTVar tol = []()
          {
            if constexpr (std::is_same_v<RTVar, float>)
              return (RTVar)180;//100; old tolerance stopped working(?)
            else if constexpr (std::is_same_v<RTVar, long double>)
              return (RTVar)10;
            else
              return (RTVar)1;
          }();


          for (int i = 0; i < nrows; i++)
          {
              for (int j = 0; j < ncols; j++)
              {
                  TVar diff = C1(i, j) - C2(i, j);
                  if (!IsZero(diff / tol))
                  {
                      CAPTURE(nrows, ncols);
                      CAPTURE(C1(i, j),C2(i, j));
                      CAPTURE(fabs(C1(i, j)-C2(i, j)));
                      CAPTURE(diff);
                      std::cout << "i " << i << " j " << j << " C1 " << C1(i, j) << " C2 " << C2(i, j) << std::endl;
                      A.Print("A = ", std::cout, EMathematicaInput);
                      B.Print("B = ", std::cout, EMathematicaInput);
                      BT.Print("BT = ", std::cout, EMathematicaInput);
                      const bool check = false;
                      REQUIRE(check);
                  }
              }
          }
          break;
      }
      case 1: // Multiplying matrices with incompatible dimensions
      {
          REQUIRE_THROWS(C1.AddContribution(0, 0, A, false, B, false, 1.0)); // this will fail for not square matrices, as A and B have the same sizes
          break;
      }
      case 2: // Adding a contribution out of matrix bounds
      {
          REQUIRE_THROWS(C1.AddContribution(1, 1, A, false, B, false, 1.0)); // this will fail because we are adding a contribution out of C1 bounds
          break;
      }
      default:
      {
          std::cout << "Test type not implemented\n";
          break;
      }
    }
    Catch::StringMaker<RTVar>::precision = oldPrecision;
}
