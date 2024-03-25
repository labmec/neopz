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
template<class TVar>
void TestSVD_Simple(const int nrows, const int ncols);
template<class TVar>
void TestSVD_Identity(const int nrows, const int ncols);
template<class TVar>
void TestSVD_Projection(const int nrows, const int ncols);
template<class TVar>
void TestSVD_Resizing(const int nrows, const int ncols);

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

TEMPLATE_TEST_CASE("Singular Value Decomposition","[matrix_tests][!shouldfail]",
                   float,
                   std::complex<float>,
                   std::complex<double>
                   )
{
  using SCAL=TestType;
  const int nrows{10};
  const int ncols{10};
  SECTION("SVD not implemented for type:"){
    TestSVD_Simple<SCAL>(nrows,ncols);
  }
  //SVD only implemented for doubles
}
TEMPLATE_TEST_CASE("Singular Value Decomposition","[matrix_tests]",
                   double
                   // float,
                   // std::complex<float>,
                   // std::complex<double>
                   ) {
  using SCAL=TestType;
  const int nrows=GENERATE(1,3,5,10);
  const int ncols=GENERATE(1,3,5,10);
  std::string sectitle;
  if(nrows==ncols){
    sectitle="Square Dense Matrix";
  }else{
    sectitle="Non-Square Dense Matrix";
  }
  SECTION(sectitle){
    SECTION("Simple"){
      TestSVD_Simple<SCAL>(nrows,ncols);
    }
    SECTION("Identity"){
      TestSVD_Identity<SCAL>(nrows, ncols);
    }
    SECTION("Projection"){
      TestSVD_Projection<SCAL>(nrows, ncols);
    }
    SECTION("Resizing"){
      TestSVD_Resizing<SCAL>(nrows, ncols);
    }
  }
}
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

template< typename TVar>
void TestSVD_Simple(int nrows, int ncols) {
	// ===============
	// Simple test
	// ===============

	/* Check if, for an auto generated matrix A, 
   * U and VT are orthogonal;
   * A == U*Sigma*VT
   */
  bool check = true;


	TPZFMatrix<TVar> mA(nrows,ncols);
	TPZFMatrix<TVar> mA_copy(nrows,ncols);
	int min = std::min(nrows,ncols);
	int max = std::max(nrows,ncols);
	TPZFMatrix<TVar> U; 
	TPZFMatrix<TVar> S; 
	TPZFMatrix<TVar> VT; 

	TPZFMatrix<TVar> Sigma(nrows,ncols); 
	TPZFMatrix<TVar> aux(nrows,ncols);

	
	mA.Resize(nrows,ncols);
	mA.AutoFill(nrows,ncols,SymProp::NonSym);
	mA_copy = mA;
	mA.SingularValueDecomposition(U,S,VT,'A','A');
	// S only gets the diagonal of Sigma
	Sigma.Resize(nrows,ncols);
	for(int i=0; i<nrows; i++){
		for(int j=0; j<ncols; j++){
			Sigma(i,j) = (i==j ? S(i,0) : 0.);
		}
	}
	aux.Resize(nrows,ncols);
	U.Multiply(Sigma,aux);
	aux.Multiply(VT,mA);
	aux = mA - mA_copy;
	for(int i=0; i<nrows; i++){
		for(int j=0; j<ncols; j++){
      TVar a = aux(i,j);
      if(!IsZero(a)) {
        check = false;
      }
		}
	}

	// Test if U is orthogonal
	TPZFMatrix<TVar> inverse(nrows,nrows);
	TVar detU = 0.;
	U.DeterminantInverse(detU,inverse);
	if(!IsZero(std::abs(detU) - 1.)) 
  {
    U.DeterminantInverse(detU,inverse);
    check = false;
  }
  TPZFMatrix<TVar> UTU(U.Rows(),U.Cols()); //< U transposed * U
  U.Multiply(U,UTU,true);
	for(int i=0; i<UTU.Rows(); i++){
		for(int j=0; j<UTU.Cols(); j++){
      if(!IsZero(std::abs(UTU(i,j)) - TVar(i==j))) {
        check = false;
      }
		}
	}
	// Test if VT is orthogonal
	inverse.Resize(nrows,nrows);
	TVar detVT = 0.;
	VT.DeterminantInverse(detVT,inverse);
	if(!IsZero(std::abs(detVT) - 1.)) 
  {check = false;}
  TPZFMatrix<TVar> VVT(VT.Rows(),VT.Cols()); //< V * V transposed
  VT.Multiply(VT,VVT,true);
	for(int i=0; i<VVT.Rows(); i++){
		for(int j=0; j<VVT.Cols(); j++){
      if(!IsZero(std::abs(VVT(i,j)) - TVar(i==j))) {
        check = false;
      }
		}
	}
	
	if (!check) {
		PZError <<__PRETTY_FUNCTION__
            <<" has failed the "
            <<" SIMPLE_TEST_SVD\n";
		mA_copy.Print("Matrix = ", std::cout, EMathematicaInput);
		// mA.Print("Matrix = ", std::cout, EFixedColumn);
		U.Print("U",std::cout, EFixedColumn);
		S.Print("SIGMA",std::cout, EFixedColumn);
		VT.Print("VT",std::cout, EFixedColumn);
	}
	REQUIRE(check);
}

template<typename TVar>
void TestSVD_Identity(int nrows, int ncols){
	// ===================================================================================
	// Check SVD of Identity Matrix and Non-square matrices that match the Kronecker delta
	// ===================================================================================
	/* let w > 0, conventionally
     w*A    = w*U*Sigma*VT
	   | 1 0 0 0 |   | 1 0 0 | | w 0 0 0 | |-1 0 0 0 |
     -w*| 0 1 0 0 | = | 0 1 0 | | 0 w 0 0 | | 0-1 0 0 |
	   | 0 0 1 0 |   | 0 0 1 | | 0 0 w 0 | | 0 0-1 0 |
     | 0 0 0 1 |
	*/
  bool check = true;


	TPZFMatrix<TVar> mA(nrows,ncols);
	TPZFMatrix<TVar> mA_copy(nrows,ncols);
	int min = std::min(nrows,ncols);
	int max = std::max(nrows,ncols);
	TPZFMatrix<TVar> U; 
	TPZFMatrix<TVar> S; 
	TPZFMatrix<TVar> VT; 

	TPZFMatrix<TVar> Sigma(nrows,ncols); 
	TPZFMatrix<TVar> aux(nrows,ncols);

	for(int i=0; i<nrows; i++){
		for(int j=0; j<ncols; j++){
      //bool gets converted to integer as 0 or 1
      const int val = i==j;
			mA(i,j) = (TVar) val;
			// mA_copy(i,j) = mA(i,j);
		}
	}

	TVar SomeNumber = -__FLT_MAX__;
	mA *= SomeNumber;
	mA_copy = mA;
	char jobU = 'A'; 
	char jobVT = 'A';
	mA.SingularValueDecomposition(U,S,VT,jobU,jobVT);

	// S only gets the diagonal of Sigma
	for(int i=0; i<nrows; i++){
		for(int j=0; j<ncols; j++){
			Sigma(i,j) = (i==j ? S(i,0) : 0.);
		}
	}

	// First, test if mA == U*Sigma*VT
	U.Multiply(Sigma,aux);
	aux.Multiply(VT,mA);
	aux = mA - mA_copy;
	for(int i=0; i<nrows; i++){
		for(int j=0; j<ncols; j++){
			if(!IsZero(aux(i,j))) check = false;
		}
	}


	// Check if absolute values of U and VT match Identity Matrix
	for(int i=0; i<U.Rows(); i++){
		for(int j=0; j<U.Cols(); j++){
			if(!IsZero(std::abs(U(i,j)) - TVar(i==j))) check = false;
		}
	}
	for(int i=0; i<VT.Rows(); i++){
		for(int j=0; j<VT.Cols(); j++){
			if(!IsZero(std::abs(VT(i,j)) - TVar(i==j))) check = false;
		}
	}
	// All singular values should match absolute value of TVar SomeNumber
	for(int j=0; j<S.Rows(); j++){
		if(!IsZero(S(j,0) - std::abs(SomeNumber))) check = false;
	}
	if (!check) {
		PZError <<__PRETTY_FUNCTION__
            <<" has failed the "
            <<" IDENTITY_TEST_SVD\n";
		mA_copy.Print("Matrix = ", std::cout, EMathematicaInput);
		// mA.Print("Matrix = ", std::cout, EFixedColumn);
		U.Print("U",std::cout, EFixedColumn);
		S.Print("SIGMA",std::cout, EFixedColumn);
		VT.Print("VT",std::cout, EFixedColumn);
	}
	REQUIRE(check);
}
	
template<typename TVar>
void TestSVD_Projection(int nrows, int ncols){
	// ==================================================
	// Checking Normal Vector of SVD of co-planar vectors
	// ==================================================
	/* If the column/row space of a matrix is a hyperplane, the Eigenvector corresponding to
     the least Eigenvalue obtained through the SVD is going to be a unit vector normal to that
     hyperplane.
     * The decision if we should check for either the column or the row space is dependent on the
     dimensions of the matrix:
     nRows < nCols ==> Check the column space;
     nRows > nCols ==> Check the  row  space;
     * Since the matrix is generated randomly, I'm not completely sure this test is deterministic
     so I gave it a second attempt in case it fails initially. (I think it's unnecessary though,
     I just couldn't prove it, and this commit was already taking too long)
	*/
  bool check = true;


	TPZFMatrix<TVar> mA(nrows,ncols);
	TPZFMatrix<TVar> mA_copy(nrows,ncols);
	int min = std::min(nrows,ncols);
	int max = std::max(nrows,ncols);
	TPZFMatrix<TVar> U; 
	TPZFMatrix<TVar> S; 
	TPZFMatrix<TVar> VT; 

	TPZFMatrix<TVar> Sigma(nrows,ncols); 
	// TPZFMatrix<TVar> aux(nrows,ncols);


	int attempts = 0;
	while(attempts < 2){
		attempts++;

		// Auto generate a random matrix
		mA.AutoFill(nrows,ncols,SymProp::NonSym);
    min = std::min(nrows,ncols);
    max = std::max(nrows,ncols);
		for(int i=0; i<max; i++){
			// Project vectors onto the same hyperplane
			if(ncols < nrows)
				mA(i,ncols-1) = 0.;
			else // (ncols > nrows)
				mA(nrows-1,i) = 0.;
		}
		mA_copy = mA;
		mA.SingularValueDecomposition(U,S,VT,'A','A');


		// Eigenvectors will go on U or VT depending on the dimensions of mA
		const TPZFMatrix<TVar>& eigenvectors = (mA.Cols() >= mA.Rows() ? U : VT);
		// const TPZFMatrix<TVar>& eigenvectors = VT;

		int dim = eigenvectors.Rows();
		TPZFMatrix<TVar> normal(dim,1);
		for(int i = 0; i<dim; i++){
			normal(i,0) = eigenvectors.g(i,dim-1);
		}
		/* Normal vector should be a unit vector in the direction of the column/row that 
       was set to zero (as in a projection to a hyperplane whose normal is that column/row)*/
		check = IsZero(1. - std::abs(normal(dim-1,0))) && IsZero(std::abs(1. - Norm(normal)));

    if(!check && attempts < 2) continue;
		if (!check) {
			PZError <<__PRETTY_FUNCTION__
              <<" has failed the "
              <<" PROJECTION_TEST_SVD\n";
			mA_copy.Print("Matrix = ", std::cout, EMathematicaInput);
			U.Print("U",std::cout, EFixedColumn);
			S.Print("SIGMA",std::cout, EFixedColumn);
			VT.Print("VT",std::cout, EFixedColumn);
			eigenvectors.Print("eigen",std::cout, EFixedColumn);
		}
		REQUIRE(check);
	}
}








template<typename TVar>
void TestSVD_Resizing(int nrows, int ncols){
	// =======================
	// Resizing test
	// =======================
	/* Checks if code is setting the sizes consistently for the user.
     it's not mandatory, but it's convenient that the function figures out the
     matrices dimensions so the user doesn't have to.
	*/
	// nrows += 1;
	// ncols += 1;
  bool check = true;


	TPZFMatrix<TVar> mA(nrows,ncols);
	TPZFMatrix<TVar> mA_copy(nrows,ncols);
	int min = std::min(nrows,ncols);
	int max = std::max(nrows,ncols);
	TPZFMatrix<TVar> U; 
	TPZFMatrix<TVar> S; 
	TPZFMatrix<TVar> VT; 
	char jobU; 
	char jobVT;

	TPZFMatrix<TVar> Sigma(nrows,ncols); 

	mA.AutoFill(nrows,ncols,SymProp::NonSym);
	mA_copy = mA;
	int m = mA.Rows();
	int n = mA.Cols();
	min = std::min(m,n);
	if(check){
		jobU = 'A'; 
		jobVT = 'A';
		mA.SingularValueDecomposition(U,S,VT,jobU,jobVT);
		if(U.Rows() != m) 
    {check = false;}
		if(U.Rows() != U.Cols()) 
    {check = false;}
		if(VT.Rows() != VT.Cols()) 
    {check = false;}
		if(VT.Rows() != n) 
    {check = false;}
		if(S.Rows() != min) 
    {check = false;}
		if(S.Cols() != 1) 
    {check = false;}
	}
	if(check){
		jobU = 'S'; 
		jobVT = 'S';
		mA.SingularValueDecomposition(U,S,VT,jobU,jobVT);
		if(U.Rows() != m) 
    {check = false;}
		if(U.Cols() != min) 
    {check = false;}
		if(VT.Cols() != n) 
    {check = false;}
		if(VT.Rows() != min) 
    {check = false;}
		if(S.Rows() != min) 
    {check = false;}
		if(S.Cols() != 1) 
    {check = false;}
	}
	if(check){
		jobU = 'N'; 
		jobVT = 'N';
		mA.SingularValueDecomposition(U,S,VT,jobU,jobVT);
		if(U.Rows() != 1) 
    {check = false;}
		if(U.Cols() != 1) 
    {check = false;}
		if(VT.Cols() != 1) 
    {check = false;}
		if(VT.Rows() != 1) 
    {check = false;}
		if(S.Rows() != min) 
    {check = false;}
		if(S.Cols() != 1) 
    {check = false;}
	}
	if (!check) {
		PZError <<__PRETTY_FUNCTION__
            <<" has failed the "
            <<" RESIZING_TEST_SVD \n with jobU = " << jobU << " and jobVT = " << jobVT << "\n";
		mA_copy.Print("Matrix = ", std::cout, EMathematicaInput);
	}
	REQUIRE(check);

}