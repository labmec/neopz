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
  
/*
Test Dot Product and Norm
*/
template <class TVar>
void TestingDotNorm(int dim);
/*
Test Add operation
*/
template <class matx, class TVar>
void TestingAdd(int row, int col, int symmetric);
/*
Test Subtract operation
*/
template <class matx, class TVar>
void TestingSubtract(int row, int col, int symmetric);
/*
Test MultiplyByScalar operation
*/
template <class matx, class TVar>
void TestingMultiplyByScalar(int row, int col, int symmetric);

template <class matx, class TVar>
void TestingAddOperator(int row, int col, int symmetric);
/*
Test Subtract operation
*/
template <class matx, class TVar>
void TestingSubtractOperator(int row, int col, int symmetric);
/*
Test MultiplyByScalar operation
*/
template <class matx, class TVar>
void TestingMultiplyByScalarOperator(int row, int col, int symmetric);
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

/**
 * @brief Tests the Eigenvalues/eigenvectors of a couple of known matrices
 * @note Process: uses LAPACK's routine */
template <class TVar>
void BasicEigenTests();


/**
 * @brief Tests the Eigenvalues/eigenvectors of the standard eigenproblem Av=wv to any matrix types. It uses the AutoFill method to create a square matrix with
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
                TestingInverseWithAutoFill<TPZFMatrix<TVar>,TVar>(dim, 1, ELU);
            }
            SECTION("TPZFBMatrix"){
                TestingInverseWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, 0, ELU);
                TestingInverseWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, 1, ELU);
            }
            SECTION("TPZBlockDiagonal"){
                TestingInverseWithAutoFill<TPZBlockDiagonal<TVar>,TVar>(dim, 0,ELU);
                TestingInverseWithAutoFill<TPZBlockDiagonal<TVar>,TVar>(dim, 1,ELU);
            }
        
            SECTION("TPZFNMatrix"){
                TestingInverseWithAutoFill<TPZFNMatrix<9,TVar>,TVar>(dim, 0, ELU);
                TestingInverseWithAutoFill<TPZFNMatrix<9,TVar>,TVar>(dim, 1, ELU);
            }
        
            SECTION("TPZSkylNSymMatrix"){
                TestingInverseWithAutoFill<TPZSkylNSymMatrix<TVar>,TVar>(dim, 0,ELU);
                TestingInverseWithAutoFill<TPZSkylNSymMatrix<TVar>,TVar>(dim, 1,ELU);
            }
#ifdef PZ_USING_MKL
            if constexpr (std::is_same<RTVar,double>::value){
                SECTION("TPZFYsmpMatrix"){
                    TestingInverseWithAutoFill<TPZFYsmpMatrix<TVar>,TVar>(dim, 0, ELU);
                }
                SECTION("TPZFYsmpMatrix"){
                    TestingInverseWithAutoFill<TPZFYsmpMatrix<TVar>,TVar>(dim, 1, ELU);
                }
            }
            
#endif            
        }
    }

  template<class TVar>
  void DotNorm(){
    for (auto dim = 5; dim < 100; dim += 500) {
      SECTION("TPZFMatrix"){
        TestingDotNorm<TVar>(dim);
      }
    }
  }

  
  template<class TVar>
    void Add(){
      for (auto dim = 5; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingAdd<TPZFMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSFMatrix"){
              TestingAdd<TPZSFMatrix<TVar>, TVar>(dim, dim, 1);
          }
          SECTION("TPZFBMatrix"){
              TestingAdd<TPZFBMatrix<TVar>,TVar>(dim, dim, 0);
          }
          SECTION("TPZSBMatrix"){
              TestingAdd<TPZSBMatrix<TVar>,TVar>(dim, dim,1);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingAdd<TPZFYsmpMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingAdd<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, 1);
          }
          // SECTION("TPZSkylNSymMatrix"){
          //     TestingAdd<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, 0);
          // }
          SECTION("TPZSkylMatrix"){
              TestingAdd<TPZSkylMatrix<TVar>, TVar>(dim, dim, 1);
          }
      }
    }

  template<class TVar>
    void Subtract(){
      for (auto dim = 5; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingSubtract<TPZFMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSFMatrix"){
              TestingSubtract<TPZSFMatrix<TVar>, TVar>(dim, dim, 1);
          }
          SECTION("TPZFBMatrix"){
              TestingSubtract<TPZFBMatrix<TVar>,TVar>(dim, dim, 0);
          }
          SECTION("TPZSBMatrix"){
              TestingSubtract<TPZSBMatrix<TVar>,TVar>(dim, dim,1);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingSubtract<TPZFYsmpMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingSubtract<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, 1);
          }
          // SECTION("TPZSkylNSymMatrix"){
          //     TestingSubtract<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, 0);
          // }
          SECTION("TPZSkylMatrix"){
              TestingSubtract<TPZSkylMatrix<TVar>, TVar>(dim, dim, 1);
          }
      }
    }
template<class TVar>
    void MultiplyByScalar(){
      for (auto dim = 5; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingMultiplyByScalar<TPZFMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSFMatrix"){
              TestingMultiplyByScalar<TPZSFMatrix<TVar>, TVar>(dim, dim, 1);
          }
          SECTION("TPZFBMatrix"){
              TestingMultiplyByScalar<TPZFBMatrix<TVar>,TVar>(dim, dim, 0);
          }
          SECTION("TPZSBMatrix"){
              TestingMultiplyByScalar<TPZSBMatrix<TVar>,TVar>(dim, dim,1);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingMultiplyByScalar<TPZFYsmpMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingMultiplyByScalar<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, 1);
          }
          // SECTION("TPZSkylNSymMatrix"){
          //     TestingMultiplyByScalar<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, 0);
          // }
          SECTION("TPZSkylMatrix"){
              TestingMultiplyByScalar<TPZSkylMatrix<TVar>, TVar>(dim, dim, 1);
          }
      }
    }
template<class TVar>
    void AddOperator(){
      for (auto dim = 5; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingAddOperator<TPZFMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSFMatrix"){
              TestingAddOperator<TPZSFMatrix<TVar>, TVar>(dim, dim, 1);
          }
          SECTION("TPZFBMatrix"){
              TestingAddOperator<TPZFBMatrix<TVar>,TVar>(dim, dim, 0);
          }
          SECTION("TPZSBMatrix"){
              TestingAddOperator<TPZSBMatrix<TVar>,TVar>(dim, dim,1);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingAddOperator<TPZFYsmpMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingAddOperator<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, 1);
          }
          // SECTION("TPZSkylNSymMatrix"){
          //     TestingAddOperator<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, 0);
          // }
          SECTION("TPZSkylMatrix"){
              TestingAddOperator<TPZSkylMatrix<TVar>, TVar>(dim, dim, 1);
          }
      }
    }

  template<class TVar>
    void SubtractOperator(){
      for (auto dim = 5; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingSubtractOperator<TPZFMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSFMatrix"){
              TestingSubtractOperator<TPZSFMatrix<TVar>, TVar>(dim, dim, 1);
          }
          SECTION("TPZFBMatrix"){
              TestingSubtractOperator<TPZFBMatrix<TVar>,TVar>(dim, dim, 0);
          }
          SECTION("TPZSBMatrix"){
              TestingSubtractOperator<TPZSBMatrix<TVar>,TVar>(dim, dim,1);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingSubtractOperator<TPZFYsmpMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingSubtractOperator<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, 1);
          }
          // SECTION("TPZSkylNSymMatrix"){
          //     TestingSubtractOperator<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, 0);
          // }
          SECTION("TPZSkylMatrix"){
              TestingSubtractOperator<TPZSkylMatrix<TVar>, TVar>(dim, dim, 1);
          }
      }
    }
template<class TVar>
    void MultiplyByScalarOperator(){
      for (auto dim = 5; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingMultiplyByScalarOperator<TPZFMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSFMatrix"){
              TestingMultiplyByScalarOperator<TPZSFMatrix<TVar>, TVar>(dim, dim, 1);
          }
          SECTION("TPZFBMatrix"){
              TestingMultiplyByScalarOperator<TPZFBMatrix<TVar>,TVar>(dim, dim, 0);
          }
          SECTION("TPZSBMatrix"){
              TestingMultiplyByScalarOperator<TPZSBMatrix<TVar>,TVar>(dim, dim,1);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingMultiplyByScalarOperator<TPZFYsmpMatrix<TVar>, TVar>(10, dim, 0);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingMultiplyByScalarOperator<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, 1);
          }
          // SECTION("TPZSkylNSymMatrix"){
          //     TestingMultiplyByScalarOperator<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, 0);
          // }
          SECTION("TPZSkylMatrix"){
              TestingMultiplyByScalarOperator<TPZSkylMatrix<TVar>, TVar>(dim, dim, 1);
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

TEMPLATE_TEST_CASE("DotNorm","[matrix_tests]",
                   float,
                   double,
                   long double,
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::DotNorm<TestType>();
}

TEMPLATE_TEST_CASE("Add","[matrix_tests]",
                   float,
                   double,
                   long double,
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::Add<TestType>();
}

TEMPLATE_TEST_CASE("Subtract","[matrix_tests]",
                   float,
                   double,
                   long double,
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::Subtract<TestType>();
}

TEMPLATE_TEST_CASE("MultiplyByScalar","[matrix_tests]",
                   float,
                   double,
                   long double,
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::MultiplyByScalar<TestType>();
}

TEMPLATE_TEST_CASE("AddOperator","[matrix_tests]",
                   float,
                   double,
                   long double,
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::AddOperator<TestType>();
}

TEMPLATE_TEST_CASE("SubtractOperator","[matrix_tests]",
                   float,
                   double,
                   long double,
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::SubtractOperator<TestType>();
}

TEMPLATE_TEST_CASE("MultiplyByScalarOperator","[matrix_tests]",
                   float,
                   double,
                   long double,
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::MultiplyByScalarOperator<TestType>();
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
    testmatrix::BasicEigenTests<TestType>();//just for real types
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
                sum += fabs(matr.Get(i, j));
        }
        CAPTURE(fabs(matr.Get(i,i)));
        CAPTURE(sum);
        REQUIRE(fabs(matr.Get(i, i)) > sum);
        if (!(fabs(matr.Get(i, i)) > sum)) {
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
            CAPTURE(i,j,mat.Get(i,j));
            REQUIRE(mat.Get(i,j).imag() == (RTVar)0);
        }
        j++;
        for(; j < ncols; j++){
            CAPTURE(i,j,mat.Get(i,j));
            REQUIRE(mat.Get(i,j)-std::conj(mat.Get(j,i)) == (TVar)0);
        }
    }
}

template <class matx, class TVar>
void TestingInverseWithAutoFill(int dim, int symmetric, DecomposeType dec) {
  int i, j;
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
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
              fabs(cpma.Get(i, j) - res.Get(i, j));
        
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
    Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
    // Checking whether both matrices are equal
    bool check = true;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            RTVar diff = fabs(square.Get(i, j) - square2.Get(i, j));
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
    Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
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

template <class TVar>
void TestingDotNorm(int dim){
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
    REQUIRE(std::abs(dot1) == Approx(norm1*norm1).margin(tol));
    REQUIRE(std::abs(dot1) == Approx(norm2*norm2).margin(tol));
    REQUIRE(std::abs(dot1-dot2) == Approx(0).margin(tol));
  }
  
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}

  
template <class matx, class TVar>
void TestingAdd(int row, int col, int symmetric) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma1, res;
  ma1.AutoFill(row, col, symmetric);
  //this is to ensure that they have the same sparsity pattern
  matx ma2(ma1);

  
  ma1.Add(ma2,res);
    
  bool check = true;
  constexpr RTVar tol = [](){
    if constexpr (std::is_same_v<RTVar,float>) return (RTVar)10;
    else return (RTVar)1;
  }();
  for (int i = 0; i < col; i++) {
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

template <class matx, class TVar>
void TestingSubtract(int row, int col, int symmetric) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma1, res;
  ma1.AutoFill(row, col, symmetric);
  //this is to ensure that they have the same sparsity pattern
  matx ma2(ma1);

  
  ma1.Subtract(ma2,res);
    
  bool check = true;
  constexpr RTVar tol = [](){
    if constexpr (std::is_same_v<RTVar,float>) return (RTVar)10;
    else return (RTVar)1;
  }();
  for (int i = 0; i < col; i++) {
    for (int j = 0; (j < col) && check; j++) {
      CAPTURE(i,j,res.Get(i,j));
      if (!IsZero((res.Get(i, j)-(ma1.Get(i,j)-ma2.Get(i,j)))/tol)) {
        check = false;
      }
      REQUIRE(check);
    }
  }
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}

template <class matx, class TVar>
void TestingMultiplyByScalar(int row, int col, int symmetric) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma1, res;
  ma1.AutoFill(row, col, symmetric);
  
  const TVar val = [](){
    constexpr RTVar lower_bound = 0;
    constexpr RTVar upper_bound = 1;
    static std::uniform_real_distribution<RTVar> unif(lower_bound,upper_bound);
    static std::default_random_engine re;
    const TVar a_random = (TVar)unif(re);
    return a_random;
  }();

  
  ma1.MultiplyByScalar(val,res);
    
  bool check = true;
  constexpr RTVar tol = [](){
    if constexpr (std::is_same_v<RTVar,float>) return (RTVar)10;
    else return (RTVar)1;
  }();
  for (int i = 0; i < col; i++) {
    for (int j = 0; (j < col) && check; j++) {
      CAPTURE(i,j,res.Get(i,j),ma1.Get(i,j)*val);
      if (!IsZero((res.Get(i, j)-(ma1.Get(i,j)*val))/tol)) {
        check = false;
      }
      REQUIRE(check);
    }
  }
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}

template <class matx, class TVar>
void TestingAddOperator(int row, int col, int symmetric) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma1;
  ma1.AutoFill(row, col, symmetric);
  //this is to ensure that they have the same sparsity pattern
  matx ma2(ma1);

  
  
  TPZMatrix<TVar> &&res = ma1 + ma2;
    
  bool check = true;
  constexpr RTVar tol = [](){
    if constexpr (std::is_same_v<RTVar,float>) return (RTVar)10;
    else return (RTVar)1;
  }();
  for (int i = 0; i < col; i++) {
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

template <class matx, class TVar>
void TestingSubtractOperator(int row, int col, int symmetric) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma1;
  ma1.AutoFill(row, col, symmetric);
  //this is to ensure that they have the same sparsity pattern
  matx ma2(ma1);

  
  TPZMatrix<TVar>&& res = ma1-ma2;
    
  bool check = true;
  constexpr RTVar tol = [](){
    if constexpr (std::is_same_v<RTVar,float>) return (RTVar)10;
    else return (RTVar)1;
  }();
  for (int i = 0; i < col; i++) {
    for (int j = 0; (j < col) && check; j++) {
      CAPTURE(i,j,res.Get(i,j));
      if (!IsZero((res.Get(i, j)-(ma1.Get(i,j)-ma2.Get(i,j)))/tol)) {
        check = false;
      }
      REQUIRE(check);
    }
  }
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}

template <class matx, class TVar>
void TestingMultiplyByScalarOperator(int row, int col, int symmetric) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma1;
  ma1.AutoFill(row, col, symmetric);
  
  const TVar val = [](){
    constexpr RTVar lower_bound = 0;
    constexpr RTVar upper_bound = 1;
    static std::uniform_real_distribution<RTVar> unif(lower_bound,upper_bound);
    static std::default_random_engine re;
    const TVar a_random = (TVar)unif(re);
    return a_random;
  }();

  
  TPZMatrix<TVar>&& res = ma1 * val;
    
  bool check = true;
  constexpr RTVar tol = [](){
    if constexpr (std::is_same_v<RTVar,float>) return (RTVar)10;
    else return (RTVar)1;
  }();
  for (int i = 0; i < col; i++) {
    for (int j = 0; (j < col) && check; j++) {
      CAPTURE(i,j,res.Get(i,j),ma1.Get(i,j)*val);
      if (!IsZero((res.Get(i, j)-(ma1.Get(i,j)*val))/tol)) {
        check = false;
      }
      REQUIRE(check);
    }
  }
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
            REQUIRE(IsZero(ma.Get(i, j) - matransptransp.Get(i, j)));
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
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma,mb;
  ma.AutoFill(dim, dim, symmetric);
  mb.AutoFill(dim, dim, symmetric);

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
  if (sizeof(RTVar) == 4) {
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

  
template <class TVar>
void BasicEigenTests() {
  if constexpr (std::is_same_v<TVar,CTVar>) return;
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
    const auto epsilon = std::numeric_limits<TVar>::epsilon()/
      (10*std::numeric_limits<TVar>::digits10);
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
  
template <class matx, class TVar>
void TestingEigenDecompositionAutoFill(int dim, int symmetric) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma;
  ma.AutoFill(dim, dim, symmetric);

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
  if (sizeof(RTVar) == 4) {
    mult *= 10.;
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
    TMatrix UTU(U.Rows(),U.Cols()); //< U transposed * U
    U.Multiply(U,UTU,true);
	for(int i=0; i<UTU.Rows(); i++){
		for(int j=0; j<UTU.Cols(); j++){
			if(!IsZero(std::abs(UTU(i,j)) - TVar(i==j))) check = false;
		}
	}
	// Test if VT is orthogonal
	inverse.Resize(nrows,nrows);
	TVar detVT = 0.;
	VT.DeterminantInverse(detVT,inverse);
	if(!IsZero(std::abs(detVT) - 1.)) 
		{check = false;}
    TMatrix VVT(VT.Rows(),VT.Cols()); //< V * V transposed
    VT.Multiply(VT,VVT,true);
	for(int i=0; i<VVT.Rows(); i++){
		for(int j=0; j<VVT.Cols(); j++){
			if(!IsZero(std::abs(VVT(i,j)) - TVar(i==j))) check = false;
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
