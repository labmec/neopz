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


namespace testmatrix{
/**
 * @brief Tests whether a matrix matr is diagonally dominant. For fixed i, checks the condition |Aii| > sum_j(|Aij|) on j!=i.
 * @param matr Matrix to check, it can to be non square matrix.
 * @note matx is a class of the type matrix matr.
 */
template <class matx>
void CheckDiagonalDominantMatrix(matx &matr);

template <class matx>
void TestGeneratingDiagonalDominantMatrix(bool isSymmetric);

template <class matx, class TVar>
void TestGeneratingHermitianMatrix();
    
/**
 * @brief Tests the Inverse method of the matrix to any matrix types. It uses the AutoFill method to create a square matrix with
 * @param dim Dimension of the square matrix to be build.
 * @param symmetric Whether to build a symmetric matrix
 * @param dec Decomposition method to be used (See enum DecomposeType at pzmatrix.h)
 * @note Process: build square matrix with randomic values, compute its inverse and the inverse of the inverse.
 * Then, checks whether the first and last matrices are identical.
 */

/** Thanks it:
 * NOTE1: I knowed that Inverse function change data of the current matrix. It stores some decomposition (LU or Cholesky ...)
 * NOTE2: If a matrix is decomposed, all times when you call Inverse, it isn't calculated because the matrix has decomposed flag as TRUE
 * NOTE3: The first function is for square matrices and the second function is for rectangular matrices.
 */
template <class matx, class TVar>
void TestingInverseWithAutoFill(int dim, int symmetric, DecomposeType dec);
    
/**
 * @brief Tests the operator * of the matrix, using AutoFill to build a square matrix of dimension dim (user defined)
 * @param dim Dimension of the square matrix to be build.
 * @param symmetric Whether to build a symmetric matrix
 * @note Process: build square matrix with randomic values, compute its square twice using * operator with two different copies of itself and checks whether the result is the same
 */
template <class matx, class TVar>
void TestingMultiplyOperatorWithAutoFill(int dim, int symmetric);
/**
 * @brief Tests the Multiply method of the matrix, using AutoFill to build a square matrix of dimension dim (user defined)
 * @param dim Dimension of the square matrix to be build.
 * @param symmetric Whether to build a symmetric matrix
 * @note Process: build square matrix with randomic values, compute its square twice using Multiply function with two different copies of itself and checks whether the result is the same
 */

/**
 * Check if the result is a identity matrix.
 */
template <class matx, class TVar>
void TestingMultiplyWithAutoFill(int dim, int symmetric);
/**
 * @brief Tests the Multiply method of the matrix, using AutoFill to build a square matrix of dimension dim (user defined)
 * @param dim Dimension of the square matrix to be build.
 * @param symmetric Whether to build a symmetric matrix
 * @note Process: build square matrix with randomic values, compute its square twice using Multiply function with two different copies of itself and checks whether the result is the same
 */

/**
 * Check if the result is a identity matrix.
 */
template <class matx, class TVar>
void TestingTransposeMultiply(int row, int col, int symmetric);
/**
 * @brief Tests the Transpose method of the matrix, using AutoFill to build a matrix of dimension row x cols (user defined)
 * @param rows Number of rows
 * @param cols Number of columns
 * @param symmetric Whether to build a symmetric matrix
 * @note Process: build matrix with randomic values, compute its transpose and transpose again then compare the first and last matrices.
 */
template <class matx, class TVar>
void TestingTransposeWithAutoFill(int rows, int cols, int symmetric);
/**
 * @brief Tests the MultAdd method of the matrix, that calculates z = beta * y + alpha * A * x , using AutoFill to build a square matrix of dimension dim (user defined)
 * @param dim Dimension of the square matrix to be build.
 * @param symmetric Whether to build a symmetric matrix
 * @param dec Decomposition method to be used (See enum DecomposeType at pzmatrix.h)
 * @note Process: build square matrix with randomic values, compute its inverse and uses MultAdd for calculating I - A*A^-1 and check whether you have any non-zero entries on the result
 */
template <class matx, class TVar>
void TestingMultAdd(int dim, int symmetric, DecomposeType dec);
#ifdef PZ_USING_LAPACK

/**
 * @brief Tests the Eigenvalues/eigenvectors of the generalised eigenproblem Av=wBv to any matrix types. It uses the AutoFill method to create a square matrix with
 * @param dim Dimension of the square matrix to be build.
 * @param symmetric Whether to build a symmetric matrix
 * @note Process: uses LAPACK's routine */
template <class matx, class TVar>
void TestingGeneralisedEigenValuesWithAutoFill(int dim, int symmetric);
/// Testing Eigenvalues of a matrix

/**
 * @brief Tests the Eigenvalues/eigenvectors of the generalised eigenproblem Av=wBv to any matrix types. It uses the AutoFill method to create a square matrix with
 * @param dim Dimension of the square matrix to be build.
 * @param symmetric Whether to build a symmetric matrix
 * @note Process: uses LAPACK's routine */
template <class matx, class TVar>
void TestingEigenDecompositionAutoFill(int dim, int symmetric);

/** @brief Check if, for an auto generated matrix A, 
* U and VT are orthogonal;
* A == U*Sigma*VT **/
template <class TMatrix, typename TVar>
void TestSVD_Simple(int nrows, int ncols);

/** @brief Checks if SVD conserves Identity Matrix
* U and VT are orthogonal;
* A == U*Sigma*VT **/
template<class TMatrix, typename TVar>
void TestSVD_Identity(int nrows, int ncols);

/** @brief If the column/row space of a matrix is a hyperplane, the Eigenvector corresponding 
 * to the least Eigenvalue obtained through the SVD should be a unit vector normal to that
 * hyperplane. **/
template<class TMatrix, typename TVar>
void TestSVD_Projection(int nrows, int ncols);

/** @brief Checks if implementation is setting the sizes consistently for the user.
 * it's not mandatory, but it's convenient that the function figures out the
 * matrices dimensions so the user doesn't have to. **/
template<class TMatrix, typename TVar>
void TestSVD_Resizing(int nrows, int ncols);

/** @brief Runs all SVD tests**/
template<class TMatrix, typename TVar>
void TestSVD(int nrows, int ncols);

#endif
    
    template<typename TVar>
    void Inverse(){
        constexpr int dim{10};

        SECTION("Cholesky"){
#ifdef PZ_USING_MKL
            if constexpr (std::is_same<RTVar,double>::value){
                SECTION("TPZSYsmpMatrix"){
                    TestingInverseWithAutoFill<TPZSYsmpMatrix<TVar>,TVar>(dim, 1, ECholesky);
                }
            }
#endif
            SECTION("TPZFMatrix"){
                TestingInverseWithAutoFill<TPZFMatrix<TVar>,TVar>(dim, 1, ECholesky);
            }
            SECTION("TPZSFMatrix"){
                TestingInverseWithAutoFill<TPZSFMatrix<TVar>,TVar>(dim, 1, ECholesky);
            }
            SECTION("TPZFBMatrix"){
                TestingInverseWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, 1, ECholesky);
            }
            SECTION("TPZSBMatrix"){
                TestingInverseWithAutoFill<TPZSBMatrix<TVar>,TVar>(dim, 1,ECholesky);
            }
            SECTION("TPZSkylMatrix"){
                TestingInverseWithAutoFill<TPZSkylMatrix<TVar>,TVar>(dim, 1,ECholesky);
            }
            SECTION("TPZNSymSkylMatrix"){
                TestingInverseWithAutoFill<TPZSkylMatrix<TVar>,TVar>(dim, 1,ECholesky);
            }
            
        }
        SECTION("LDLt"){
#ifdef PZ_USING_MKL
            if constexpr (std::is_same<RTVar,double>::value){
                SECTION("TPZSYsmpMatrix"){
                    TestingInverseWithAutoFill<TPZSYsmpMatrix<TVar>,TVar>(dim, 1, ELDLt);
                }
            }
#endif
            SECTION("TPZFMatrix"){
                TestingInverseWithAutoFill<TPZFMatrix<TVar>,TVar>(dim, 1, ELDLt);
            }
            SECTION("TPZSFMatrix"){
                TestingInverseWithAutoFill<TPZSFMatrix<TVar>,TVar>(dim, 1, ELDLt);
            }
            SECTION("TPZFBMatrix"){
                TestingInverseWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, 1, ECholesky);
            }
            SECTION("TPZSBMatrix"){
                TestingInverseWithAutoFill<TPZSBMatrix<TVar>,TVar>(dim, 1, ELDLt);
            }
            SECTION("TPZSkylMatrix"){
                TestingInverseWithAutoFill<TPZSkylMatrix<TVar>,TVar>(dim, 1,ELDLt);
            }
            SECTION("TPZNSymSkylMatrix"){
                TestingInverseWithAutoFill<TPZSkylMatrix<TVar>,TVar>(dim, 1,ELDLt);
            }
        }
        SECTION("LU"){
            SECTION("TPZFMatrix"){
                TestingInverseWithAutoFill<TPZFMatrix<TVar>,TVar>(dim, 0, ELU);
            }
            SECTION("TPZFBMatrix"){
                TestingInverseWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, 0, ELU);
            }
            SECTION("TPZBlockDiagonal"){
                TestingInverseWithAutoFill<TPZBlockDiagonal<TVar>,TVar>(dim, 0,ELU);
            }
        
            SECTION("TPZFNMatrix"){
                TestingInverseWithAutoFill<TPZFNMatrix<9,TVar>,TVar>(dim, 0, ELU);
            }
        
            SECTION("TPZSkylNSymMatrix"){
                TestingInverseWithAutoFill<TPZSkylNSymMatrix<TVar>,TVar>(dim, 0,ELU);
            }
        }
    }

    template<class TVar>
    void MultiplyTranspose(){
      for (auto dim = 5; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingTransposeMultiply<TPZFMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSFMatrix"){
              TestingTransposeMultiply<TPZSFMatrix<TVar>, TVar>(dim, dim, 1);
          }
          SECTION("TPZFBMatrix"){
              TestingTransposeMultiply<TPZFBMatrix<TVar>,TVar>(dim, dim, 0);
          }
          SECTION("TPZSBMatrix"){
              TestingTransposeMultiply<TPZSBMatrix<TVar>,TVar>(dim, dim,1);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingTransposeMultiply<TPZFYsmpMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingTransposeMultiply<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, 1);
          }
          SECTION("TPZSkylNSymMatrix"){
              TestingTransposeMultiply<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, 0);
          }
          SECTION("TPZSkylMatrix"){
              TestingTransposeMultiply<TPZSkylMatrix<TVar>, TVar>(dim, dim, 1);
          }
      }
    }

    template<class TVar>
    void Multiply(){
      for (auto dim = 20; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingMultiplyWithAutoFill<TPZFMatrix<TVar>, TVar>(dim, 0);
          }
          SECTION("TPZSFMatrix"){
              TestingMultiplyWithAutoFill<TPZSFMatrix<TVar>, TVar>(dim, 1);
          }
          SECTION("TPZFBMatrix"){
              TestingMultiplyWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, 0);
          }
          SECTION("TPZSBMatrix"){
              TestingMultiplyWithAutoFill<TPZSBMatrix<TVar>,TVar>(dim,1);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingMultiplyWithAutoFill<TPZFYsmpMatrix<TVar>, TVar>(dim, 0);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingMultiplyWithAutoFill<TPZSYsmpMatrix<TVar>, TVar>(dim, 1);
          }
          SECTION("TPZSkylNSymMatrix"){
              TestingMultiplyWithAutoFill<TPZSkylNSymMatrix<TVar>, TVar>(dim, 0);
          }
          SECTION("TPZSkylMatrix"){
              TestingMultiplyWithAutoFill<TPZSkylMatrix<TVar>, TVar>(dim, 1);
          }
      }
    }


    template<class TVar>
    void Hermitian(){
        SECTION("TPZFMatrix"){
            TestGeneratingHermitianMatrix<TPZFMatrix<TVar>,TVar>();
        }        
        SECTION("TPZFBMatrix"){
            TestGeneratingHermitianMatrix<TPZFBMatrix<TVar>,TVar>();
        }
        SECTION("TPZSkylNSymMatrix"){
            TestGeneratingHermitianMatrix<TPZSkylNSymMatrix<TVar>,TVar>();
        }
        SECTION("TPZFYsmpMatrix"){
            TestGeneratingHermitianMatrix<TPZFYsmpMatrix<TVar>,TVar>();
        }

        SECTION("TPZSFMatrix"){
            TestGeneratingHermitianMatrix<TPZSFMatrix<TVar>,TVar>();
        }
        SECTION("TPZSBMatrix"){
            TestGeneratingHermitianMatrix<TPZSBMatrix<TVar>,TVar>();
        }
        SECTION("TPZSkylMatrix"){
            TestGeneratingHermitianMatrix<TPZSkylMatrix<TVar>,TVar>();
        }
        SECTION("TPZSYsmpMatrix"){
            TestGeneratingHermitianMatrix<TPZSYsmpMatrix<TVar>,TVar>();
        }
    }
    
    template<class TVar>
    void DiagonalDominant(){
        bool isSymmetric{false};
        SECTION("TPZFMatrix"){            
            TestGeneratingDiagonalDominantMatrix<TPZFMatrix<TVar>>(isSymmetric);
        }
        SECTION("TPZFBMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZFBMatrix<TVar>>(isSymmetric);
        }
        SECTION("TPZSkylNSymMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZSkylNSymMatrix<TVar>>(isSymmetric);
        }
        SECTION("TPZFYsmpMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZFYsmpMatrix<TVar>>(isSymmetric);
        }
        isSymmetric = true;
        SECTION("TPZSFMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZSFMatrix<TVar>>(isSymmetric);
        }
        SECTION("TPZSBMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZSBMatrix<TVar>>(isSymmetric);
        }
        SECTION("TPZSkylMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZSkylMatrix<TVar>>(isSymmetric);
        }
        SECTION("TPZSYsmpMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZSYsmpMatrix<TVar>>(isSymmetric);
        }
        SECTION("TPZBlockDiagonal"){
            // Unit Test for block diagonal matrix
            TPZVec<int> blocks(13);
            for (int i = 0; i < 13; i++)
                blocks[i] = 15 + (i % 4);
            TPZBlockDiagonal<TVar> mabd(blocks);
            mabd.AutoFill(50, 50, 0);
            CheckDiagonalDominantMatrix(mabd);
        }
    }
    
    template <class TVar>
    void MultiplyOperatorWithAutoFill() {
      for (int dim = 3; dim < 100; dim += 5) {
        TestingMultiplyOperatorWithAutoFill<TPZFMatrix<TVar>, TVar>(dim, 0);
      }
    }
    
    template <class TVar>
    void TransposeWithAutoFill() {
      for (int rows = 3; rows < 4; rows += 5) {
        for (int cols = 3; cols < 100; cols += 5) {
          TestingTransposeWithAutoFill<TPZFMatrix<TVar>, TVar>(rows, cols, 0);
        }
      }
    }

    template<class TVar>
    void TestMultAdd(){
      for (int dim = 5; dim < 6; dim += 10) {
          SECTION("TPZFMatrix"){
              TestingMultAdd<TPZFMatrix<TVar>, TVar>(dim, 0, ELU);
          }
          SECTION("TPZSFMatrix"){
              TestingMultAdd<TPZSFMatrix<TVar>, TVar>(dim, 1, ECholesky);
          }
          SECTION("TPZFBMatrix"){
              TestingMultAdd<TPZFBMatrix<TVar>, TVar>(dim, 0, ELU);
          }
          SECTION("TPZSBMatrix"){
              TestingMultAdd<TPZSBMatrix<TVar>, TVar>(dim, 1, ECholesky);
          }
          SECTION("TPZSkylMatrix"){
              TestingMultAdd<TPZSkylMatrix<TVar>, TVar>(dim, 1, ECholesky);
          }
          SECTION("TPZSkylNSymMatrix"){
              TestingMultAdd<TPZSkylNSymMatrix<TVar>, TVar>(dim, 0, ELU);
          }
      }
    }
#ifdef PZ_USING_LAPACK
    template <class TVar> void GeneralisedEigenvaluesAutoFill() {
        for (int dim = 5; dim < 6; dim += 10) {
            SECTION("TPZFMatrix"){
                TestingGeneralisedEigenValuesWithAutoFill<TPZFMatrix<TVar>, TVar>(dim, 1);
            }
            SECTION("TPZSBMatrix"){
                TestingGeneralisedEigenValuesWithAutoFill<TPZSBMatrix<TVar>, TVar>(dim, 1);
            }
        }
    }
    
    template<class TVar>
    void EigenDecompositionAutoFill(){
        for (int dim = 5; dim < 6; dim += 10) {
            SECTION("TPZSBMatrix sym"){
                TestingEigenDecompositionAutoFill<TPZSBMatrix<TVar>, TVar>(dim, 1);
            }
            SECTION("TPZFMatrix sym"){
                TestingEigenDecompositionAutoFill<TPZFMatrix<TVar>, TVar>(dim, 1);
            }
            SECTION("TPZFMatrix nsym"){
                TestingEigenDecompositionAutoFill<TPZFMatrix<TVar>, TVar>(dim, 0);
            }
        }
    }
    
    template<class TVar>
    void SingularValueDecomposition_Real(){
        SECTION("Square Dense Matrix"){
            TestSVD<TPZFMatrix<TVar>, TVar>(1,1);
            TestSVD<TPZFMatrix<TVar>, TVar>(5,5);
        }
        SECTION("Non-Square Dense Matrix"){
            TestSVD<TPZFMatrix<TVar>, TVar>(6,4);
            TestSVD<TPZFMatrix<TVar>, TVar>(3,7);
        }
    }
#endif
};


TEMPLATE_TEST_CASE("Inverse (REAL)","[matrix_tests]",
                   float,
                   double
                   ,long double
                   ) {
    testmatrix::Inverse<TestType>();
}

TEMPLATE_TEST_CASE("Inverse (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
          ) {
    testmatrix::Inverse<TestType>();
}

TEMPLATE_TEST_CASE("Multiply transpose (REAL)","[matrix_tests]",
                   float,
                   double,
                   long double
                   ) {
    testmatrix::MultiplyTranspose<TestType>();
}
TEMPLATE_TEST_CASE("Multiply transpose (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::MultiplyTranspose<TestType>();
}

TEMPLATE_TEST_CASE("Multiply (REAL)","[matrix_tests]",
                   float,
                   double,
                   long double
                   ) {
    testmatrix::Multiply<TestType>();
}

TEMPLATE_TEST_CASE("Multiply (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::Multiply<TestType>();
}

TEMPLATE_TEST_CASE("Diagonal dominant (REAL)","[matrix_tests]",
                   float,
                   double,
                   long double
                   ) {
    for(int i = 0; i < 5; i++)
        testmatrix::DiagonalDominant<TestType>();
}

TEMPLATE_TEST_CASE("Diagonal dominant (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    for(int i = 0; i < 5; i++)
        testmatrix::DiagonalDominant<TestType>();
}

TEMPLATE_TEST_CASE("Hermitian (REAL)","[matrix_tests]",
                   float,
                   double,
                   long double
                   ) {
    for(int i = 0; i < 5; i++)
        testmatrix::Hermitian<TestType>();
}

TEMPLATE_TEST_CASE("Hermitian (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    for(int i = 0; i < 5; i++)
        testmatrix::Hermitian<TestType>();
}

TEMPLATE_TEST_CASE("Multiply operator (REAL)","[matrix_tests]",
                   float,
                   double,
                   long double
                   ) {
    testmatrix::MultiplyOperatorWithAutoFill<TestType>();
}

TEMPLATE_TEST_CASE("Multiply operator (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::MultiplyOperatorWithAutoFill<TestType>();
}

TEMPLATE_TEST_CASE("Transpose (REAL)","[matrix_tests]",
                   float,
                   double,
                   long double
                   ) {
    testmatrix::TransposeWithAutoFill<TestType>();
}

TEMPLATE_TEST_CASE("Transpose (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::TransposeWithAutoFill<TestType>();
}


TEMPLATE_TEST_CASE("MultAdd (REAL)","[matrix_tests]",
                   float,
                   double,
                   long double
                   ) {
    testmatrix::TestMultAdd<TestType>();
}

TEMPLATE_TEST_CASE("MultAdd (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::TestMultAdd<TestType>();
}

#ifdef PZ_USING_LAPACK
/*There is no long double lapack interface in our code*/
TEMPLATE_TEST_CASE("Eigenvalues (REAL)","[matrix_tests]",
                   float,
                   double
                   ) {
    
    testmatrix::EigenDecompositionAutoFill<TestType>();
}

TEMPLATE_TEST_CASE("Eigenvalues (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>
                   ) {
    
    testmatrix::EigenDecompositionAutoFill<TestType>();
}

TEMPLATE_TEST_CASE("Generalised Eigenvalues (REAL)","[matrix_tests]",
                   float,
                   double
                   ) {
    testmatrix::GeneralisedEigenvaluesAutoFill<TestType>();
}

TEMPLATE_TEST_CASE("Generalised Eigenvalues (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>
                   ) {
    testmatrix::GeneralisedEigenvaluesAutoFill<TestType>();
}

TEMPLATE_TEST_CASE("Singular Value Decomposition (REAL)","[matrix_tests]",
                //    float,
                   double
                   ) {
    testmatrix::SingularValueDecomposition_Real<TestType>();
}
#endif



namespace testmatrix{

    
template <class matx>
void CheckDiagonalDominantMatrix(matx &matr) {
    REAL sum;
    for (int i = 0; i < matr.Rows(); i++) {
        sum = 0.0;
        for (int j = 0; j < matr.Cols(); j++) {
            if (i != j)
                sum += fabs(matr.GetVal(i, j));
        }
        CAPTURE(fabs(matr.GetVal(i,i)));
        CAPTURE(sum);
        REQUIRE(fabs(matr.GetVal(i, i)) > sum);
        if (!(fabs(matr.GetVal(i, i)) > sum)) {
            std::cout << "line i " << i << " failed\n";
            // matr.Print("matrix = ", std::cout, EMathematicaInput);
            return;
        }
    }
}

template<class matx>
void TestGeneratingDiagonalDominantMatrix(bool symmetric){
    for(int dim = 3; dim<100;dim+=7){
        const auto nrows = dim;
        const auto ncols = dim;
        matx mat;
        mat.AutoFill(nrows,ncols,symmetric);
        CheckDiagonalDominantMatrix(mat);
    }
}
    
template <class matx, class TVar>
void TestGeneratingHermitianMatrix() {
    const auto nrows = 10;
    const auto ncols = 10;
    constexpr bool symmetric{true};
    matx mat;
    mat.AutoFill(nrows,ncols,symmetric);
    for (auto i = 0; i < nrows; i++) {
        auto j = i;
        if constexpr (is_complex<TVar>::value){
            CAPTURE(i,j,mat.GetVal(i,j));
            REQUIRE(mat.GetVal(i,j).imag() == (RTVar)0);
        }
        j++;
        for(; j < ncols; j++){
            CAPTURE(i,j,mat.GetVal(i,j));
            REQUIRE(mat.GetVal(i,j)-std::conj(mat.GetVal(j,i)) == (TVar)0);
        }
    }
}

template <class matx, class TVar>
void TestingInverseWithAutoFill(int dim, int symmetric, DecomposeType dec) {
  int i, j;
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  matx ma;
  ma.AutoFill(dim, dim, symmetric);  
  // Making ma copy because ma is modified by Inverse method (it's decomposed)
  matx cpma(ma);
  TPZFMatrix<TVar> inv(dim, dim), invkeep;
  TPZFMatrix<TVar> res(inv);
  // getting inverse twice
  ma.Inverse(inv, dec);
  invkeep = inv;
  inv.Inverse(res, dec);
  //    ma.Print("skyl2 =",std::cout,EMathematicaInput);
  bool check = true;
  /// Checking whether the res matrix is identical to m1 matrix
  for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
          RTVar diff =
              fabs(cpma.GetVal(i, j) - res.GetVal(i, j));
        
          bool loccheck = IsZero(diff / (RTVar) 10.);
          if (loccheck == false) {
              CAPTURE(i, j, diff);
              std::cout << "diff " << diff << std::endl;
          }
          check &= loccheck;
      }
  }
    
  REQUIRE(check);
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}

template <class matx, class TVar>
void TestingMultiplyOperatorWithAutoFill(int dim, int symmetric) {
    // ma times inv must to be a identity matrix
    matx ma;
    ma.AutoFill(dim, dim, symmetric);

    TPZFMatrix<TVar> duplicate(ma);
    TPZFMatrix<TVar> square, square2;

    square2 = duplicate*duplicate;
    square = ma*duplicate;
    auto oldPrecision = Catch::StringMaker<RTVar>::precision;
    // Checking whether both matrices are equal
    bool check = true;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            RTVar diff = fabs(square.GetVal(i, j) - square2.GetVal(i, j));
            if (!IsZero(diff)) {
                CAPTURE(i,j,diff);
                check = false;
            }
        }
    }
    REQUIRE(check);
    Catch::StringMaker<RTVar>::precision = oldPrecision;
}

template <class matx, class TVar>
void TestingMultiplyWithAutoFill(int dim, int symmetric) {
    // ma times inv must to be a identity matrix
    matx ma;
    ma.AutoFill(dim, dim, symmetric);

    TPZFMatrix<TVar> duplicate(ma), square, square2;

    //    ma.Print("SkylineNS");
    //    duplicate.Print("FullMat");
    ma.Multiply(duplicate, square);
    duplicate.Multiply(duplicate, square2);
    auto oldPrecision = Catch::StringMaker<RTVar>::precision;
    // Checking whether result matrix is the identity matrix
    bool check = true;
    constexpr RTVar tol = [](){
        if constexpr (std::is_same_v<RTVar,float>) return (RTVar)100;
        else return (RTVar)1;
    }();
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            if (!IsZero(fabs(square(i, j) - square2(i, j))/tol)) {
                CAPTURE(i,j,square(i,j));
                check = false;
            }
        }
    }
    if(!check){
        std::cout<<typeid(matx).name()<<std::endl;
        std::cout << __PRETTY_FUNCTION__ << "failed \n";
    }

    REQUIRE(check);
    Catch::StringMaker<RTVar>::precision = oldPrecision;
}


template <class matx, class TVar>
void TestingTransposeMultiply(int row, int col, int symmetric) {
    // ma times inv must to be a identity matrix
    matx ma;
    ma.AutoFill(row, col, symmetric);

    TPZFMatrix<TVar> duplicate(ma);
    TPZFMatrix<TVar> dup2(ma), square(col, col), square2(col, col);
    ma.MultAdd(dup2, dup2, square, 1., 0., 1);
    duplicate.MultAdd(dup2, dup2, square2, 1., 0., 1);

    //    ma.Print("ma",std::cout,EMathematicaInput);
    //    duplicate.Print("duplicate",std::cout,EMathematicaInput);
    //    square.Print("square",std::cout,EMathematicaInput);
    //    square2.Print("square2",std::cout,EMathematicaInput);
    // Checking whether result matrix is the identity matrix
    bool check = true;
    constexpr RTVar tol = [](){
        if constexpr (std::is_same_v<RTVar,float>) return (RTVar)10;
        else return (RTVar)1;
    }();
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < col; j++) {
            if (!IsZero(fabs(square(i, j) - square2(i, j))/tol)) {
                check = false;
            }
        }
    }
    if (!check) {
        std::cout<<typeid(matx).name()<<std::endl;
        std::cout << __PRETTY_FUNCTION__ << "failed \n";
        // duplicate.Print("Full");
        // ma.Print("Matrix");
        // duplicate -=ma;
        // duplicate.Print("Diff");
        // square2.Print("Full");
        // square.Print("Matrix");
        // square2 -= square;
        // square2.Print("Diff");
    }
    REQUIRE(check);

}

template <class matx, class TVar>
void TestingTransposeWithAutoFill(int rows, int cols, int symmetric) {
    int i, j;

    matx ma;
    ma.AutoFill(rows, cols, symmetric);

    matx matransp(cols, rows);
    matx matransptransp(ma);

    // getting inverse twice
    matransptransp.Transpose(&matransp);
    matransp.Transpose(&matransptransp);

    /// Checking whether the res matrix is identical to m1 matrix
    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            REQUIRE(IsZero(ma.GetVal(i, j) - matransptransp.GetVal(i, j)));
}

template <class matx, class TVar>
void TestingMultAdd(int dim, int symmetric, DecomposeType dec) {
    int i, j;

    matx ma;
    ma.AutoFill(dim, dim, symmetric);

    TPZFMatrix<TVar> cpy(ma);
    TPZFMatrix<TVar> inv(dim, dim);
    TPZFMatrix<TVar> y(ma), z;
    y.Identity();
    // getting inverse twice
    cpy.Inverse(inv, dec);

    //    virtual void MultAdd(const TPZFMatrix<TVar> & x,const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z,
    //                         const TVar alpha=1., const TVar beta = 0., const int opt = 0) const;


    ma.MultAdd(inv, y, z, 1., -1.);
    constexpr RTVar tol = [](){
        if constexpr (std::is_same_v<RTVar,float>) return (RTVar)100;
        else return (RTVar)1;
    }();
    /// Checking whether the res matrix is identical to m1 matrix
    bool check = true;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            TVar zval = z(i, j);
            if (!IsZero(zval/tol)) {
                check = false;
            }
        }
    }
    REQUIRE(check);
}

#ifdef PZ_USING_LAPACK

template <class matx, class TVar>
void TestingGeneralisedEigenValuesWithAutoFill(int dim, int symmetric) {

    matx ma, mb;
    ma.AutoFill(dim, dim, symmetric);
    mb.AutoFill(dim, dim, symmetric);

    // Making ma and mb copies because ma and mb are modified by the eigenvalues method
    bool check = true;
    matx cpmaOriginal(ma);
    matx cpma(ma), cpmb(mb);
    TPZFMatrix<std::complex<double> > cpfma(dim, dim), cpfmb(dim, dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            cpfma(i, j) = ma(i, j);
            cpfmb(i, j) = mb(i, j);
        }
    }

    TPZVec < std::complex<double> > w;
    TPZFMatrix < std::complex<double> > eigenVectors;
    cpma.SolveGeneralisedEigenProblem(cpmb, w, eigenVectors);
    //    ma.Print("a",std::cout , EMathematicaInput);
    //    mb.Print("b",std::cout , EMathematicaInput);
    //    std::cout<<"w = ";
    //    std::cout<<w<<std::endl;
    //    eigenVectors.Print("eV",std::cout , EMathematicaInput);
    double mult = 10;
    if (sizeof (RTVar) == 4) {
        mult *= 10.;
    }
    for (int i = 0; i < dim; i++) {
        TPZFMatrix< std::complex<double> > res(dim, dim, 0.);
        TPZFMatrix< std::complex<double> > x(dim, 1, 0.);
        eigenVectors.GetSub(0, i, dim, 1, x);

        res = cpfma * x - w[i] * cpfmb * x;
        for (int j = 0; j < dim; j++) {
            bool loccheck = IsZero(RTVar(res(j, 0).real() / mult)) && IsZero(RTVar(res(j, 0).imag() / mult));
            if (loccheck == false) {
                std::cout << "diff " << res(j, 0) << std::endl;
            }
            check &= loccheck;
        }
    }


    if (!check) {
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" has failed\n";
        
        cpmaOriginal.Print("Matrix = ", std::cout, EMathematicaInput);
    }
    REQUIRE(check);
}

template <class matx, class TVar>
void TestingEigenDecompositionAutoFill(int dim, int symmetric) {
matx ma;
    ma.AutoFill(dim, dim, symmetric);
    //    ma.Print("a",std::cout , EMathematicaInput);

    // Making ma and mb copies because ma and mb are modified by the eigenvalues method
    bool check = true;
    matx cpmaOriginal(ma);
    //    if(symmetric){
    TPZFMatrix< std::complex<double> > cpma(dim, dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            cpma(i, j) = ma.GetVal(i, j);
        }
    }

    TPZVec < std::complex<double> > w;
    TPZFMatrix < std::complex<double> > eigenVectors;
    ma.SolveEigenProblem(w, eigenVectors);

    //    cpma.Print("a = ",std::cout , EMathematicaInput);
    //    std::cout<<w<<std::endl;
    //    eigenVectors.Print("eV = ",std::cout,EMathematicaInput);
    double mult = 10.;
    if (sizeof (RTVar) == 4) {
        mult *= 10.;
    }
    for (int i = 0; i < dim; i++) {
        TPZFMatrix< std::complex<double> > res(dim, 1, 0.);
        TPZFMatrix< std::complex<double> > x(dim, 1, 0.);
        for (int j = 0; j < dim; j++) {
            x(j, 0) = eigenVectors(j, i);
        }
        for (int k = 0; k < dim; k++) {
            res(k, 0) = -(w[i]) * x(k, 0);
            for (int l = 0; l < dim; l++) {
                res(k, 0) += cpma(k, l) * x(l, 0);
            }
        }
        for (int j = 0; j < dim; j++) {
            bool loccheck = IsZero(RTVar(res(j, 0).real() / mult)) && IsZero(RTVar(res(j, 0).imag() / mult));
            if (loccheck == false) {
                std::cout << "diff " << res(j, 0) << "i = " << i << " j = " << j << std::endl;
            }
            check &= loccheck;
        }
    }
    //    }

    if (!check) {
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" has failed\n";
        cpmaOriginal.Print("Matrix = ", std::cout, EMathematicaInput);
    }
    REQUIRE(check);
}


template<class TMatrix, typename TVar>
void TestSVD_Simple(int nrows, int ncols) {
	// ===============
	// Simple test
	// ===============

	/* Check if, for an auto generated matrix A, 
	* U and VT are orthogonal;
	* A == U*Sigma*VT
	*/
	using namespace std;
	bool check = true;


	TMatrix mA(nrows,ncols);
	TMatrix mA_copy(nrows,ncols);
	int min = std::min(nrows,ncols);
	int max = std::max(nrows,ncols);
	TMatrix U; 
	TMatrix S; 
	TMatrix VT; 

	TMatrix Sigma(nrows,ncols); 
	TMatrix aux(nrows,ncols);

	
	mA.Resize(nrows,ncols);
	mA.AutoFill(nrows,ncols,false);
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
			if(!IsZero(aux(i,j))) check = false;
		}
	}

	// Test if U is orthogonal
	TMatrix inverse(nrows,nrows);
	TVar detU = 0.;
	U.DeterminantInverse(detU,inverse);
	if(!IsZero(std::abs(detU) - 1.)) 
		{check = false;}
	// Test if VT is orthogonal
	inverse.Resize(nrows,nrows);
	TVar detVT = 0.;
	VT.DeterminantInverse(detVT,inverse);
	if(!IsZero(std::abs(detVT) - 1.)) 
		{check = false;}
	
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










template<class TMatrix, typename TVar>
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
	using namespace std;
	bool check = true;


	TMatrix mA(nrows,ncols);
	TMatrix mA_copy(nrows,ncols);
	int min = std::min(nrows,ncols);
	int max = std::max(nrows,ncols);
	TMatrix U; 
	TMatrix S; 
	TMatrix VT; 

	TMatrix Sigma(nrows,ncols); 
	TMatrix aux(nrows,ncols);

	for(int i=0; i<nrows; i++){
		for(int j=0; j<ncols; j++){
			mA(i,j) = (TVar) i==j;
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
	












template<class TMatrix, typename TVar>
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
	using namespace std;
	bool check = true;


	TMatrix mA(nrows,ncols);
	TMatrix mA_copy(nrows,ncols);
	int min = std::min(nrows,ncols);
	int max = std::max(nrows,ncols);
	TMatrix U; 
	TMatrix S; 
	TMatrix VT; 

	TMatrix Sigma(nrows,ncols); 
	// TMatrix aux(nrows,ncols);


	int attempts = 0;
	while(attempts < 2){
		attempts++;

		// Auto generate a random matrix
		mA.AutoFill(nrows,ncols,false);
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
		const TMatrix& eigenvectors = (mA.Cols() >= mA.Rows() ? U : VT);
		// const TMatrix& eigenvectors = VT;

		int dim = eigenvectors.Rows();
		TMatrix normal(dim,1);
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








template<class TMatrix, typename TVar>
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
	using namespace std;
	bool check = true;


	TMatrix mA(nrows,ncols);
	TMatrix mA_copy(nrows,ncols);
	int min = std::min(nrows,ncols);
	int max = std::max(nrows,ncols);
	TMatrix U; 
	TMatrix S; 
	TMatrix VT; 
	char jobU; 
	char jobVT;

	TMatrix Sigma(nrows,ncols); 

	mA.AutoFill(nrows,ncols,false);
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
		PZError << "\nYou can fix your implementation with the following snippet:\n\n\n"
				<< "// Setup matrix sizes. This section is not mandatory, but it's good to do the work for the user\n"
				<< "int m = this->Rows();\n"
				<< "int n = this->Cols();\n"
				<< "int min = std::min(m,n);\n"
				<< "switch(jobU){\n"
				<< "    case 'A': U.Resize(m,m);    break;\n"
				<< "    case 'S': U.Resize(m,min);  break;\n"
				<< "    case 'N': U.Resize(1,1);    break;\n"
				<< "    default: PZError << \"\nInvalid input jobU = \'\" << jobU << \"\' unrecognized;\n\"; DebugStop();\n"
				<< "}\n"
				<< "switch(jobVT){\n"
				<< "    case 'A': VT.Resize(n,n);   break;\n"
				<< "    case 'S': VT.Resize(min,n); break;\n"
				<< "    case 'N': VT.Resize(1,1);   break;\n"
				<< "    default: PZError << \"\nInvalid input jobVT = \'\" << jobVT << \"\' unrecognized;\n\"; DebugStop();\n"
				<< "}\n"
				<< "S.Resize(min,1);\n\n\n";
		mA_copy.Print("Matrix = ", std::cout, EMathematicaInput);
	}
	REQUIRE(check);

}

template<class TMatrix, typename TVar>
void TestSVD(int nrows, int ncols){
	TestSVD_Simple<TMatrix, TVar >(nrows,ncols);
	TestSVD_Identity<TMatrix, TVar >(nrows,ncols);
	TestSVD_Projection<TMatrix, TVar >(nrows,ncols);
	TestSVD_Resizing<TMatrix, TVar >(nrows,ncols);
}

#endif // PZ_USING_LAPACK

};
