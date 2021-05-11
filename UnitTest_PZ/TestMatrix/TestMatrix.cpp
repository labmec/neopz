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
            if constexpr(!is_complex<TVar>::value){//FAILING
                SECTION("TPZSFMatrix"){
                    TestingInverseWithAutoFill<TPZSFMatrix<TVar>,TVar>(dim, 1, ECholesky);
                }
            }
            SECTION("TPZSBMatrix"){
                TestingInverseWithAutoFill<TPZSBMatrix<TVar>,TVar>(dim, 1,ECholesky);
            }
            if constexpr(!is_complex<TVar>::value){//FAILING
                SECTION("TPZSkylMatrix"){
                    TestingInverseWithAutoFill<TPZSkylMatrix<TVar>,TVar>(dim, 1,ECholesky);
                }
            }
            
        }
        //FAILING
        if constexpr(!is_complex<TVar>::value){
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
            SECTION("TPZSBMatrix"){
                TestingInverseWithAutoFill<TPZSBMatrix<TVar>,TVar>(dim, 1, ELDLt);
            }
            SECTION("TPZSkylMatrix"){
                TestingInverseWithAutoFill<TPZSkylMatrix<TVar>,TVar>(dim, 1,ELDLt);
            }
        }
        }
        SECTION("LU"){
            SECTION("TPZFBMatrix"){
                TestingInverseWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, 0, ELU);
            }
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
#endif
};


TEMPLATE_TEST_CASE("Inverse (REAL)","[matrix_tests]",
                   float,
                   double,
                   long double
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
    testmatrix::DiagonalDominant<TestType>();
}

TEMPLATE_TEST_CASE("Diagonal dominant (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::DiagonalDominant<TestType>();
}

TEMPLATE_TEST_CASE("Hermitian (REAL)","[matrix_tests]",
                   float,
                   double,
                   long double
                   ) {
    testmatrix::Hermitian<TestType>();
}

TEMPLATE_TEST_CASE("Hermitian (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
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

// TEMPLATE_TEST_CASE("MultAdd (CPLX)","[matrix_tests]",
//                    std::complex<float>,
//                    std::complex<double>,
//                    std::complex<long double>
//                    ) {
//     testmatrix::TestMultAdd<TestType>();
// }

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
#endif

//TEST_CASE("nonsingular_test","[matrix_tests]")
//{
//	TestingDiagonalDominant<TPZFBMatrix>();
//}



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
            REQUIRE(mat.GetVal(i,j).imag() == (RTVar)0);
        }
        j++;
        for(; j < ncols; j++){
            REQUIRE(mat.GetVal(i,j)-std::conj(mat.GetVal(j,i)) == (TVar)0);
        }
    }
}

template <class matx, class TVar>
void TestingInverseWithAutoFill(int dim, int symmetric, DecomposeType dec) {
  int i, j;

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
      // if (loccheck == false) {
      //     std::cout << "diff " << diff << std::endl;
      // }
      check &= loccheck;
    }
    }
    if (!check) {
        std::cout << __PRETTY_FUNCTION__;
        std::cout << "failed with dec type: ";
        switch (dec) {
        case ELU: std::cout << "LU\n"; break;
        case ELDLt: std::cout << "LDLt\n"; break;
        case ECholesky: std::cout << "Cholesky\n"; break;
        case ELUPivot: std::cout << "LU Pivot\n"; break;
        case ENoDecompose: DebugStop();
        }
        std::cout << std::flush;
        // cpma.Print("Matrix = ", std::cout, EMathematicaInput);
        // invkeep.Print("Inv = ", std::cout, EMathematicaInput);
    }
    REQUIRE(check);
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
    // Checking whether both matrices are equal
    bool check = true;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            TVar diff = fabs(square.GetVal(i, j) - square2.GetVal(i, j));
            if (!IsZero(diff)) {
                check = false;
            }
        }
    }
    REQUIRE(check);
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
    // Checking whether result matrix is the identity matrix
    bool check = true;
    constexpr RTVar tol = [](){
        if constexpr (std::is_same_v<RTVar,float>) return (RTVar)100;
        else return (RTVar)1;
    }();
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            if (!IsZero(fabs(square(i, j) - square2(i, j))/tol)) {
                check = false;
            }
        }
    }
    if(!check){
        std::cout<<typeid(matx).name()<<std::endl;
        std::cout << __PRETTY_FUNCTION__ << "failed \n";
    }

    REQUIRE(check);

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
    /// Checking whether the res matrix is identical to m1 matrix
    bool check = true;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            TVar zval = z(i, j);
            if (!IsZero(zval)) {
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
#endif

};
