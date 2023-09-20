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
#include "TPZYSMPMatrix.h"
#include "TPZSYSMPMatrix.h"
#include "pzmatred.h"
#ifdef PZ_USING_MKL
#include "TPZYSMPPardiso.h"
#include "TPZSYSMPPardiso.h"
#endif
#include "tpzsparseblockdiagonal.h"
#include "pzseqsolver.h"
#include "pzstepsolver.h"
#include <random>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>

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
void TestGeneratingDiagonalDominantMatrix(SymProp sp);

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
void TestingInverseWithAutoFill(int dim, SymProp sp, DecomposeType dec);

/**
 * @brief Tests the Inverse method of the matrix to any matrix types.
 * @note The matrix should have been filled.
 * @param mat Matrix to be inverted
 * @param dec Decomposition method to be used (See enum DecomposeType at pzmatrix.h)
 * @note Process: gets square matrix, compute its inverse and the inverse of the inverse.
 * Then, checks whether the first and last matrices are identical.
 */
template <class matx, class TVar>
void TestingInverse(matx mat, DecomposeType dec);
    
/**
 * @brief Tests the operator * of the matrix, using AutoFill to build a square matrix of dimension dim (user defined)
 * @param dim Dimension of the square matrix to be build.
 * @param symmetric Whether to build a symmetric matrix
 * @note Process: build square matrix with randomic values, compute its square twice using * operator with two different copies of itself and checks whether the result is the same
 */
template <class matx, class TVar>
void TestingMultiplyOperatorWithAutoFill(int dim, SymProp sp);
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
void TestingMultiplyWithAutoFill(int dim, SymProp sp);
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
void TestingAdd(int row, int col, SymProp sp);
/*
Test Subtract operation
*/
template <class matx, class TVar>
void TestingSubtract(int row, int col, SymProp sp);
/*
Test MultiplyByScalar operation
*/
template <class matx, class TVar>
void TestingMultiplyByScalar(int row, int col, SymProp sp);

template <class matx, class TVar>
void TestingAddOperator(int row, int col, SymProp sp);
/*
Test Subtract operation
*/
template <class matx, class TVar>
void TestingSubtractOperator(int row, int col, SymProp sp);
/*
Test MultiplyByScalar operation
*/
template <class matx, class TVar>
void TestingMultiplyByScalarOperator(int row, int col, SymProp sp);
/**
 * Check if the result is a identity matrix.
 */
template <class matx, class TVar>
void TestingTransposeMultiply(int row, int col, SymProp sp);
/**
 * @brief Tests the Transpose method of the matrix, using AutoFill to build a matrix of dimension row x cols (user defined)
 * @param rows Number of rows
 * @param cols Number of columns
 * @param symmetric Whether to build a symmetric matrix
 * @note Process: build matrix with randomic values, compute its transpose and transpose again then compare the first and last matrices.
 */
template <class matx, class TVar>
void TestingTransposeWithAutoFill(int rows, int cols, SymProp sp);
/**
 * @brief Tests the MultAdd method of the matrix, that calculates z = beta * y + alpha * A * x , using AutoFill to build a square matrix of dimension dim (user defined)
 * @param dim Dimension of the square matrix to be build.
 * @param symmetric Whether to build a symmetric matrix
 * @param dec Decomposition method to be used (See enum DecomposeType at pzmatrix.h)
 * @note Process: build square matrix with randomic values, compute its inverse and uses MultAdd for calculating I - A*A^-1 and check whether you have any non-zero entries on the result
 */
template <class matx, class TVar>
void TestingMultAdd(int dim, SymProp sp, DecomposeType dec);

/**
 * @brief Tests the addContribution method of the matrix, that adds a block C += alpha * A*B starting at C(i,j), using AutoFill to build a square matrix of dimension dim (user defined)
 * @param nrows Number of rows of the matrix to be build.
 * @param ncols Number of columns of the matrix to be build.
 * @note Process: build a matrix C with randomic values, adds a contribution C += A*B of the same size as C. Compare the results with AddContribution and MultAdd.
 */
template <class TVar>
void TestingAddContribution(int nrows, int ncols, int ntype);

/**
 * @brief Tests the MatRed::SolveDirect() method by comparing the solution obtained using reduced and full matrix
 * @param dim0 dimension of the first submatrix K00
 * @param dim dimension of the full matrix K
 * @param ntype type of the test to be performed
 * @note Process: build a matrix C with randomic values, adds a contribution C += A*B of the same size as C. Compare the results with AddContribution and MultAdd.
 */
template <class TVar>
void TestingMatRedSolver(int dim0, int dim, int ntype);

#ifdef PZ_USING_LAPACK

/**
 * @brief Tests the Eigenvalues/eigenvectors of the generalised eigenproblem Av=wBv to any matrix types. It uses the AutoFill method to create a square matrix with
 * @param dim Dimension of the square matrix to be build.
 * @param symmetric Whether to build a symmetric matrix
 * @note Process: uses LAPACK's routine */
template <class matx, class TVar>
void TestingGeneralisedEigenValuesWithAutoFill(int dim, SymProp sp);

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
void TestingEigenDecompositionAutoFill(int dim, SymProp sp);

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
                  TestingInverseWithAutoFill<TPZSYsmpMatrixPardiso<TVar>,TVar>(dim, SymProp::Herm, ECholesky);
                }
            }
#endif
            SECTION("TPZFMatrix"){
                TestingInverseWithAutoFill<TPZFMatrix<TVar>,TVar>(dim, SymProp::Herm, ECholesky);
            }
            SECTION("TPZSFMatrix"){
                TestingInverseWithAutoFill<TPZSFMatrix<TVar>,TVar>(dim, SymProp::Herm, ECholesky);
            }
//            SECTION("TPZFBMatrix"){
//                TestingInverseWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, SymProp::Herm, ECholesky);
//            }
            SECTION("TPZSBMatrix"){
                TestingInverseWithAutoFill<TPZSBMatrix<TVar>,TVar>(dim, SymProp::Herm,ECholesky);
            }
            SECTION("TPZSkylMatrix"){
                TestingInverseWithAutoFill<TPZSkylMatrix<TVar>,TVar>(dim, SymProp::Herm,ECholesky);
            }
            SECTION("TPZNSymSkylMatrix"){
                TestingInverseWithAutoFill<TPZSkylMatrix<TVar>,TVar>(dim, SymProp::Herm,ECholesky);
            }
            
        }
        SECTION("LDLt"){
#ifdef PZ_USING_MKL
            if constexpr (std::is_same<RTVar,double>::value){
                SECTION("TPZSYsmpMatrix"){
                  TestingInverseWithAutoFill<TPZSYsmpMatrixPardiso<TVar>,TVar>(dim, SymProp::Herm, ELDLt);
                }
            }
#endif
            SECTION("TPZFMatrix"){
                TestingInverseWithAutoFill<TPZFMatrix<TVar>,TVar>(dim, SymProp::Herm, ELDLt);
            }
            SECTION("TPZSFMatrix"){
                TestingInverseWithAutoFill<TPZSFMatrix<TVar>,TVar>(dim, SymProp::Herm, ELDLt);
            }
//            SECTION("TPZFBMatrix"){
//                TestingInverseWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, SymProp::Herm, ECholesky);
//            }
            SECTION("TPZSBMatrix"){
                TestingInverseWithAutoFill<TPZSBMatrix<TVar>,TVar>(dim, SymProp::Herm, ELDLt);
            }
            SECTION("TPZSkylMatrix"){
                TestingInverseWithAutoFill<TPZSkylMatrix<TVar>,TVar>(dim, SymProp::Herm,ELDLt);
            }
            SECTION("TPZNSymSkylMatrix"){
                TestingInverseWithAutoFill<TPZSkylMatrix<TVar>,TVar>(dim, SymProp::Herm,ELDLt);
            }
        }
        SECTION("LU"){
            SECTION("TPZFMatrix"){
                TestingInverseWithAutoFill<TPZFMatrix<TVar>,TVar>(dim, SymProp::NonSym, ELU);
                TestingInverseWithAutoFill<TPZFMatrix<TVar>,TVar>(dim, SymProp::Herm, ELU);
            }
            SECTION("TPZFBMatrix"){
                TestingInverseWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, SymProp::NonSym, ELU);
                TestingInverseWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, SymProp::Herm, ELU);
            }
            SECTION("TPZBlockDiagonal"){
                TestingInverseWithAutoFill<TPZBlockDiagonal<TVar>,TVar>(dim, SymProp::NonSym,ELU);
                TestingInverseWithAutoFill<TPZBlockDiagonal<TVar>,TVar>(dim, SymProp::Herm,ELU);
            }
        
            SECTION("TPZFNMatrix"){
                TestingInverseWithAutoFill<TPZFNMatrix<9,TVar>,TVar>(dim, SymProp::NonSym, ELU);
                TestingInverseWithAutoFill<TPZFNMatrix<9,TVar>,TVar>(dim, SymProp::Herm, ELU);
            }
        
            SECTION("TPZSkylNSymMatrix"){
                TestingInverseWithAutoFill<TPZSkylNSymMatrix<TVar>,TVar>(dim, SymProp::NonSym,ELU);
                TestingInverseWithAutoFill<TPZSkylNSymMatrix<TVar>,TVar>(dim, SymProp::Herm,ELU);
            }
#ifdef PZ_USING_MKL
            if constexpr (std::is_same<RTVar,double>::value){
                SECTION("TPZFYsmpMatrix"){
                    TestingInverseWithAutoFill<TPZFYsmpMatrixPardiso<TVar>,TVar>(dim, SymProp::NonSym, ELU);
                }
                SECTION("TPZFYsmpMatrix"){
                    TestingInverseWithAutoFill<TPZFYsmpMatrixPardiso<TVar>,TVar>(dim, SymProp::Herm, ELU);
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
              TestingAdd<TPZFMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSFMatrix"){
              TestingAdd<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingAdd<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          SECTION("TPZFBMatrix"){
              TestingAdd<TPZFBMatrix<TVar>,TVar>(dim, dim, SymProp::NonSym);
          }
          SECTION("TPZSBMatrix"){
              TestingAdd<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Herm);
              TestingAdd<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Sym);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingAdd<TPZFYsmpMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingAdd<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingAdd<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          // SECTION("TPZSkylNSymMatrix"){
          //     TestingAdd<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, SymProp::NonSym);
          // }
          SECTION("TPZSkylMatrix"){
              TestingAdd<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingAdd<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
      }
    }

  template<class TVar>
    void Subtract(){
      for (auto dim = 5; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingSubtract<TPZFMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSFMatrix"){
              TestingSubtract<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingSubtract<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          SECTION("TPZFBMatrix"){
              TestingSubtract<TPZFBMatrix<TVar>,TVar>(dim, dim, SymProp::NonSym);
          }
          SECTION("TPZSBMatrix"){
              TestingSubtract<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Herm);
              TestingSubtract<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Sym);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingSubtract<TPZFYsmpMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingSubtract<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingSubtract<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          // SECTION("TPZSkylNSymMatrix"){
          //     TestingSubtract<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, SymProp::NonSym);
          // }
          SECTION("TPZSkylMatrix"){
              TestingSubtract<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingSubtract<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
      }
    }
template<class TVar>
    void MultiplyByScalar(){
      for (auto dim = 5; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingMultiplyByScalar<TPZFMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSFMatrix"){
              TestingMultiplyByScalar<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingMultiplyByScalar<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          SECTION("TPZFBMatrix"){
              TestingMultiplyByScalar<TPZFBMatrix<TVar>,TVar>(dim, dim, SymProp::NonSym);
          }
          SECTION("TPZSBMatrix"){
              TestingMultiplyByScalar<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Herm);
              TestingMultiplyByScalar<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Sym);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingMultiplyByScalar<TPZFYsmpMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingMultiplyByScalar<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingMultiplyByScalar<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          // SECTION("TPZSkylNSymMatrix"){
          //     TestingMultiplyByScalar<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, SymProp::NonSym);
          // }
          SECTION("TPZSkylMatrix"){
              TestingMultiplyByScalar<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingMultiplyByScalar<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
      }
    }
template<class TVar>
    void AddOperator(){
      for (auto dim = 5; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingAddOperator<TPZFMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSFMatrix"){
              TestingAddOperator<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingAddOperator<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          SECTION("TPZFBMatrix"){
              TestingAddOperator<TPZFBMatrix<TVar>,TVar>(dim, dim, SymProp::NonSym);
          }
          SECTION("TPZSBMatrix"){
              TestingAddOperator<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Herm);
              TestingAddOperator<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Sym);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingAddOperator<TPZFYsmpMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingAddOperator<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingAddOperator<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          // SECTION("TPZSkylNSymMatrix"){
          //     TestingAddOperator<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, SymProp::NonSym);
          // }
          SECTION("TPZSkylMatrix"){
              TestingAddOperator<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingAddOperator<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
      }
    }

  template<class TVar>
    void SubtractOperator(){
      for (auto dim = 5; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingSubtractOperator<TPZFMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSFMatrix"){
              TestingSubtractOperator<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingSubtractOperator<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          SECTION("TPZFBMatrix"){
              TestingSubtractOperator<TPZFBMatrix<TVar>,TVar>(dim, dim, SymProp::NonSym);
          }
          SECTION("TPZSBMatrix"){
              TestingSubtractOperator<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Herm);
              TestingSubtractOperator<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Sym);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingSubtractOperator<TPZFYsmpMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingSubtractOperator<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingSubtractOperator<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          // SECTION("TPZSkylNSymMatrix"){
          //     TestingSubtractOperator<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, SymProp::NonSym);
          // }
          SECTION("TPZSkylMatrix"){
              TestingSubtractOperator<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingSubtractOperator<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
      }
    }
template<class TVar>
    void MultiplyByScalarOperator(){
      for (auto dim = 5; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingMultiplyByScalarOperator<TPZFMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSFMatrix"){
              TestingMultiplyByScalarOperator<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingMultiplyByScalarOperator<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          SECTION("TPZFBMatrix"){
              TestingMultiplyByScalarOperator<TPZFBMatrix<TVar>,TVar>(dim, dim, SymProp::NonSym);
          }
          SECTION("TPZSBMatrix"){
              TestingMultiplyByScalarOperator<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Herm);
              TestingMultiplyByScalarOperator<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Sym);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingMultiplyByScalarOperator<TPZFYsmpMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingMultiplyByScalarOperator<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingMultiplyByScalarOperator<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          // SECTION("TPZSkylNSymMatrix"){
          //     TestingMultiplyByScalarOperator<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, SymProp::NonSym);
          // }
          SECTION("TPZSkylMatrix"){
              TestingMultiplyByScalarOperator<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingMultiplyByScalarOperator<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
      }
    }
  
    template<class TVar>
    void MultiplyTranspose(){
      for (auto dim = 5; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingTransposeMultiply<TPZFMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSFMatrix"){
              TestingTransposeMultiply<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingTransposeMultiply<TPZSFMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          SECTION("TPZFBMatrix"){
              TestingTransposeMultiply<TPZFBMatrix<TVar>,TVar>(dim, dim, SymProp::NonSym);
          }
          SECTION("TPZSBMatrix"){
              TestingTransposeMultiply<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Herm);
              TestingTransposeMultiply<TPZSBMatrix<TVar>,TVar>(dim, dim,SymProp::Sym);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingTransposeMultiply<TPZFYsmpMatrix<TVar>, TVar>(10, dim, SymProp::NonSym);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingTransposeMultiply<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingTransposeMultiply<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          SECTION("TPZSkylNSymMatrix"){
              TestingTransposeMultiply<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim, SymProp::NonSym);
          }
          SECTION("TPZSkylMatrix"){
              TestingTransposeMultiply<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Herm);
              TestingTransposeMultiply<TPZSkylMatrix<TVar>, TVar>(dim, dim, SymProp::Sym);
          }
          SECTION("TPZBlockDiagonal"){
              TestingTransposeMultiply<TPZBlockDiagonal<TVar>, TVar>(dim, dim, SymProp::NonSym);
          }
#ifdef PZ_USING_MKL
          //suported MKL types
          if constexpr ((
                          (std::is_same_v<TVar,float>) ||
                          (std::is_same_v<TVar,double>) ||
                          (std::is_same_v<TVar,std::complex<float>>) ||
                          (std::is_same_v<TVar,std::complex<double>>)
                         )){
            SECTION("TPZFYsmpMatrixPardiso"){
              TestingTransposeMultiply<TPZFYsmpMatrixPardiso<TVar>, TVar>(dim, dim, SymProp::NonSym);
            }
            SECTION("TPZSYsmpMatrixPardiso"){
                TestingTransposeMultiply<TPZSYsmpMatrixPardiso<TVar>, TVar>(dim, dim, SymProp::Herm);
                TestingTransposeMultiply<TPZSYsmpMatrixPardiso<TVar>, TVar>(dim, dim, SymProp::Sym);
            }
          }
#endif
      }
    }

    template<class TVar>
    void Multiply(){
      for (auto dim = 20; dim < 100; dim += 500) {
          SECTION("TPZFMatrix"){
              TestingMultiplyWithAutoFill<TPZFMatrix<TVar>, TVar>(dim, SymProp::NonSym);
          }
          SECTION("TPZSFMatrix"){
              TestingMultiplyWithAutoFill<TPZSFMatrix<TVar>, TVar>(dim, SymProp::Herm);
              TestingMultiplyWithAutoFill<TPZSFMatrix<TVar>, TVar>(dim, SymProp::Sym);
          }
          SECTION("TPZFBMatrix"){
              TestingMultiplyWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, SymProp::NonSym);
          }
          SECTION("TPZSBMatrix"){
              TestingMultiplyWithAutoFill<TPZSBMatrix<TVar>,TVar>(dim,SymProp::Herm);
              TestingMultiplyWithAutoFill<TPZSBMatrix<TVar>,TVar>(dim,SymProp::Sym);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingMultiplyWithAutoFill<TPZFYsmpMatrix<TVar>, TVar>(dim, SymProp::NonSym);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingMultiplyWithAutoFill<TPZSYsmpMatrix<TVar>, TVar>(dim, SymProp::Herm);
              TestingMultiplyWithAutoFill<TPZSYsmpMatrix<TVar>, TVar>(dim, SymProp::Sym);
          }
          SECTION("TPZSkylNSymMatrix"){
              TestingMultiplyWithAutoFill<TPZSkylNSymMatrix<TVar>, TVar>(dim, SymProp::NonSym);
          }
          SECTION("TPZSkylMatrix"){
              TestingMultiplyWithAutoFill<TPZSkylMatrix<TVar>, TVar>(dim, SymProp::Herm);
              TestingMultiplyWithAutoFill<TPZSkylMatrix<TVar>, TVar>(dim, SymProp::Sym);
          }
#ifdef PZ_USING_MKL
          //suported MKL types
          if constexpr ((
                          (std::is_same_v<TVar,float>) ||
                          (std::is_same_v<TVar,double>) ||
                          (std::is_same_v<TVar,std::complex<float>>) ||
                          (std::is_same_v<TVar,std::complex<double>>)
                         )){
            SECTION("TPZFYsmpMatrixPardiso"){
              TestingMultiplyWithAutoFill<TPZFYsmpMatrixPardiso<TVar>, TVar>(dim, SymProp::NonSym);
            }
            SECTION("TPZSYsmpMatrixPardiso"){
                TestingMultiplyWithAutoFill<TPZSYsmpMatrixPardiso<TVar>, TVar>(dim, SymProp::Herm);
                TestingMultiplyWithAutoFill<TPZSYsmpMatrixPardiso<TVar>, TVar>(dim, SymProp::Sym);
            }
          }
#endif
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
        auto sp = SymProp::NonSym;
        SECTION("TPZFMatrix"){            
            TestGeneratingDiagonalDominantMatrix<TPZFMatrix<TVar>>(SymProp::NonSym);
        }
        SECTION("TPZFBMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZFBMatrix<TVar>>(SymProp::NonSym);
        }
        SECTION("TPZSkylNSymMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZSkylNSymMatrix<TVar>>(SymProp::NonSym);
        }
        SECTION("TPZFYsmpMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZFYsmpMatrix<TVar>>(SymProp::NonSym);
        }

        SECTION("TPZSFMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZSFMatrix<TVar>>(SymProp::Sym);
            TestGeneratingDiagonalDominantMatrix<TPZSFMatrix<TVar>>(SymProp::Herm);
        }
        SECTION("TPZSBMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZSBMatrix<TVar>>(SymProp::Sym);
            TestGeneratingDiagonalDominantMatrix<TPZSBMatrix<TVar>>(SymProp::Herm);
        }
        SECTION("TPZSkylMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZSkylMatrix<TVar>>(SymProp::Sym);
            TestGeneratingDiagonalDominantMatrix<TPZSkylMatrix<TVar>>(SymProp::Herm);
        }
        SECTION("TPZSYsmpMatrix"){
            TestGeneratingDiagonalDominantMatrix<TPZSYsmpMatrix<TVar>>(SymProp::Sym);
            TestGeneratingDiagonalDominantMatrix<TPZSYsmpMatrix<TVar>>(SymProp::Herm);
        }
        SECTION("TPZBlockDiagonal"){
            // Unit Test for block diagonal matrix
            TPZVec<int> blocks(13);
            for (int i = 0; i < 13; i++)
                blocks[i] = 15 + (i % 4);
            TPZBlockDiagonal<TVar> mabd(blocks);
            mabd.AutoFill(50, 50, SymProp::NonSym);
            CheckDiagonalDominantMatrix(mabd);
        }
    }
    
    template <class TVar>
    void MultiplyOperatorWithAutoFill() {
      for (int dim = 3; dim < 100; dim += 5) {
        TestingMultiplyOperatorWithAutoFill<TPZFMatrix<TVar>, TVar>(dim, SymProp::NonSym);
        TestingMultiplyOperatorWithAutoFill<TPZFNMatrix<20,TVar>, TVar>(dim, SymProp::NonSym);
      }
    }
    
    template <class TVar>
    void TransposeWithAutoFill() {
      for (int rows = 3; rows < 4; rows += 5) {
        for (int cols = 3; cols < 100; cols += 5) {
          TestingTransposeWithAutoFill<TPZFMatrix<TVar>, TVar>(rows, cols, SymProp::NonSym);
        }
      }
    }

    template<class TVar>
    void TestMultAdd(){
      for (int dim = 5; dim < 6; dim += 10) {
          SECTION("TPZFMatrix"){
              TestingMultAdd<TPZFMatrix<TVar>, TVar>(dim, SymProp::NonSym, ELU);
          }
          SECTION("TPZSFMatrix"){
              TestingMultAdd<TPZSFMatrix<TVar>, TVar>(dim, SymProp::Herm, ECholesky);
          }
          SECTION("TPZFBMatrix"){
              TestingMultAdd<TPZFBMatrix<TVar>, TVar>(dim, SymProp::NonSym, ELU);
          }
          SECTION("TPZSBMatrix"){
              TestingMultAdd<TPZSBMatrix<TVar>, TVar>(dim, SymProp::Herm, ECholesky);
          }
          SECTION("TPZSkylMatrix"){
              TestingMultAdd<TPZSkylMatrix<TVar>, TVar>(dim, SymProp::Herm, ECholesky);
          }
          SECTION("TPZSkylNSymMatrix"){
              TestingMultAdd<TPZSkylNSymMatrix<TVar>, TVar>(dim, SymProp::NonSym, ELU);
          }
          SECTION("TPZSYsmpMatrix"){
              TestingMultAdd<TPZSYsmpMatrix<TVar>, TVar>(dim, SymProp::Herm, ECholesky);
          }
          SECTION("TPZFYsmpMatrix"){
              TestingMultAdd<TPZFYsmpMatrix<TVar>, TVar>(dim, SymProp::NonSym, ELU);
          }
          SECTION("TPZBlockDiagonal"){
              TestingMultAdd<TPZBlockDiagonal<TVar>, TVar>(dim, SymProp::NonSym, ELU);
          }
#ifdef PZ_USING_MKL
          //suported MKL types
          if constexpr ((
                          (std::is_same_v<TVar,float>) ||
                          (std::is_same_v<TVar,double>) ||
                          (std::is_same_v<TVar,std::complex<float>>) ||
                          (std::is_same_v<TVar,std::complex<double>>)
                         )){
            SECTION("TPZFYsmpMatrixPardiso"){
              TestingMultAdd<TPZFYsmpMatrixPardiso<TVar>, TVar>(dim, SymProp::NonSym, ELU);
            }
            SECTION("TPZSYsmpMatrixPardiso"){
              TestingMultAdd<TPZSYsmpMatrixPardiso<TVar>, TVar>(dim, SymProp::Herm, ECholesky);
            }
          }
#endif
      }
    }

    template <class TVar>
    void TestAddContribution()
    {
      SECTION("TPZFMatrix-SQUARE-MULTADD")
      {
          TestingAddContribution<TVar>(10, 10, 0);
      }
#ifdef PZDEBUG
      SECTION("TPZFMatrix-INCOMPATIBLE-DIMENSIONS")
      {
          TestingAddContribution<TVar>(10, 8, 1);
      }
      SECTION("TPZFMatrix-OUT-OF-BOUNDS")
      {
          TestingAddContribution<TVar>(4, 4, 2);
      }
#endif
    }

    template <class TVar>
    void TestMatRedSolver()
    {
      SECTION("TPZMatRedMatrix")
      {
          TestingMatRedSolver<TVar>(2, 4, 0);
      }
    }



#ifdef PZ_USING_LAPACK
    template <class TVar> void GeneralisedEigenvaluesAutoFill() {
        for (int dim = 5; dim < 6; dim += 10) {
            SECTION("TPZFMatrix"){
                TestingGeneralisedEigenValuesWithAutoFill<TPZFMatrix<TVar>, TVar>(dim, SymProp::Sym);
                TestingGeneralisedEigenValuesWithAutoFill<TPZFMatrix<TVar>, TVar>(dim, SymProp::Herm);
            }
            SECTION("TPZSBMatrix"){
                TestingGeneralisedEigenValuesWithAutoFill<TPZSBMatrix<TVar>, TVar>(dim, SymProp::Sym);
                TestingGeneralisedEigenValuesWithAutoFill<TPZSBMatrix<TVar>, TVar>(dim, SymProp::Herm);
            }
        }
    }
    
    template<class TVar>
    void EigenDecompositionAutoFill(){
        for (int dim = 5; dim < 6; dim += 10) {
            
            SECTION("TPZSBMatrix sym"){
                TestingEigenDecompositionAutoFill<TPZSBMatrix<TVar>, TVar>(dim, SymProp::Sym);
                TestingEigenDecompositionAutoFill<TPZSBMatrix<TVar>, TVar>(dim, SymProp::Herm);
            }
            SECTION("TPZFMatrix sym"){
                TestingEigenDecompositionAutoFill<TPZFMatrix<TVar>, TVar>(dim, SymProp::Sym);
                TestingEigenDecompositionAutoFill<TPZFMatrix<TVar>, TVar>(dim, SymProp::Herm);
            }
            SECTION("TPZFMatrix nsym"){
                TestingEigenDecompositionAutoFill<TPZFMatrix<TVar>, TVar>(dim, SymProp::NonSym);
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
    TestingInverse<TPZBlockDiagonal<TVar>, TVar>(blckmat, ELU);
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
    TPZBlockDiagonal<TVar> blckmat(blocksizes,fmat);
    REQUIRE(true);
  }

  template<class TVar>
  void SparseBlockDiagInverse();
  template<class TVar>
  void SparseBlockColorInverse();
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

TEMPLATE_TEST_CASE("AddContribution (REAL)","[matrix_tests]",
                   float,
                   double,
                   long double
                   ) {
    testmatrix::TestAddContribution<TestType>();
}

TEMPLATE_TEST_CASE("AddContribution (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::TestAddContribution<TestType>();
}

TEMPLATE_TEST_CASE("MatRed SolveDirect (REAL)","[matrix_tests]",
                   float,
                   double,
                   long double
                   ) {
    testmatrix::TestMatRedSolver<TestType>();
}

TEMPLATE_TEST_CASE("MatRed SolveDirect (CPLX)","[matrix_tests]",
                   std::complex<float>,
                   std::complex<double>,
                   std::complex<long double>
                   ) {
    testmatrix::TestMatRedSolver<TestType>();
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

TEMPLATE_TEST_CASE("Additional Block Diagonal tests","[matrix_tests]",
                   float,
                   double,
                   std::complex<float>,
                   std::complex<double>
                   ) {
  SECTION("LU Pivot"){
    testmatrix::BlockDiagLUPivot<TestType>();
  }
  SECTION("BlockDiag block with zero size"){
    testmatrix::BlockDiagZeroSizedBlock<TestType>();
  }
  SECTION("Sparse Block Inverse"){
    testmatrix::SparseBlockDiagInverse<TestType>();
    testmatrix::SparseBlockColorInverse<TestType>();
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
    constexpr bool validType{false};
    REQUIRE(validType);
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
void TestGeneratingDiagonalDominantMatrix(SymProp sp){
    for(int dim = 3; dim<100;dim+=7){
        const auto nrows = dim;
        const auto ncols = dim;
        matx mat;
        mat.AutoFill(nrows,ncols,sp);
        CheckDiagonalDominantMatrix(mat);
    }
}
    
template <class matx, class TVar>
void TestGeneratingHermitianMatrix() {
    const auto nrows = 10;
    const auto ncols = 10;
    matx mat;
    mat.AutoFill(nrows,ncols,SymProp::Herm);
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
void TestingInverseWithAutoFill(int dim, SymProp sp, DecomposeType dec) {
  int i, j;
  matx ma;
  ma.AutoFill(dim, dim, sp);
  TestingInverse<matx,TVar>(ma,dec);
}

template <class matx, class TVar>
void TestingInverse(matx ma, DecomposeType dec) {
  const auto symmetric = SymPropName(ma.GetSymmetry());
  if(ma.Rows() != ma.Cols()){
    const bool is_square_mat{false};
    REQUIRE(is_square_mat);
  }
  // Making ma copy because ma is modified by Inverse method (it's decomposed)
  matx cpma(ma);
  const int dim = ma.Rows();
  TPZFMatrix<TVar> inv(dim, dim), invkeep;
  TPZFMatrix<TVar> res(inv);
  // getting inverse twice
  ma.Inverse(inv, dec);
  invkeep = inv;
  inv.Inverse(res, dec);
  //    ma.Print("skyl2 =",std::cout,EMathematicaInput);
  bool check = true;
  /// Checking whether the res matrix is identical to m1 matrix
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
          RTVar diff =
              fabs(cpma.Get(i, j) - res.Get(i, j));
        
          bool loccheck = IsZero(diff / (RTVar) 100.);
          if (loccheck == false) {
              CAPTURE(i, j, diff);
              std::cout << " i " << i << " j " << j << "diff " << diff << std::endl;
          }
          check &= loccheck;
      }
  }
  CAPTURE(dim,symmetric,dec);
    if(!check) {
        cpma.Print(" block mat ",std::cout);
        inv.Print(" inv mat ", std::cout);
        res.Print(" inv inv mat ",std::cout);
    }
  REQUIRE(check);
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}

template <class matx, class TVar>
void TestingMultiplyOperatorWithAutoFill(int dim, SymProp sp) {
    // ma times inv must to be a identity matrix
    matx ma;
    ma.AutoFill(dim, dim, sp);

    TPZFNMatrix<10,TVar> duplicate(ma);
    TPZFNMatrix<10,TVar> square, square2;

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
void TestingMultiplyWithAutoFill(int dim, SymProp sp) {
    const auto sp_str = SymPropName(sp);
    // ma times inv must to be a identity matrix
    matx ma;
    ma.AutoFill(dim, dim, sp);

    TPZFMatrix<TVar> duplicate(ma), square, square2;

    //    ma.Print("SkylineNS");
    //    duplicate.Print("FullMat");
    ma.Multiply(duplicate, square);
    duplicate.Multiply(duplicate, square2);
    auto oldPrecision = Catch::StringMaker<RTVar>::precision;
    Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
    // Checking whether result matrix is the same
    bool check = true;
    constexpr RTVar tol =std::numeric_limits<RTVar>::epsilon()*2000;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            const auto diff = fabs(square(i, j) - square2(i, j));
            CAPTURE(i,j,square(i,j), square2(i,j), diff, tol, dim, sp_str, typeid(matx).name());
            REQUIRE((diff == Catch::Approx(0).margin(tol)));
            if(diff > tol){
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
    REQUIRE(std::abs(dot1) == Catch::Approx(norm1*norm1).margin(tol));
    REQUIRE(std::abs(dot1) == Catch::Approx(norm2*norm2).margin(tol));
    REQUIRE(std::abs(dot1-dot2) == Catch::Approx(0).margin(tol));
  }
  
  Catch::StringMaker<RTVar>::precision = oldPrecision;
}

  
template <class matx, class TVar>
void TestingAdd(int row, int col, SymProp symmetric) {
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
void TestingSubtract(int row, int col, SymProp sp) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma1, res;
  ma1.AutoFill(row, col, sp);
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
void TestingMultiplyByScalar(int row, int col, SymProp sp) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma1, res;
  ma1.AutoFill(row, col, sp);
  
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
    else return (RTVar)10;
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
void TestingAddOperator(int row, int col, SymProp sp) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma1;
  ma1.AutoFill(row, col, sp);
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
void TestingSubtractOperator(int row, int col, SymProp sp) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma1;
  ma1.AutoFill(row, col, sp);
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
void TestingMultiplyByScalarOperator(int row, int col, SymProp sp) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma1;
  ma1.AutoFill(row, col, sp);
  
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
    else return (RTVar)10;
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
void TestingTransposeMultiply(int row, int col, SymProp sp) {
    matx ma;
    ma.AutoFill(row, col, sp);

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
    constexpr RTVar tol =std::numeric_limits<RTVar>::epsilon()*1000;
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < col; j++) {
            const auto diff = fabs(square(i, j) - square2(i, j));
            if(diff > tol){
              check = false;
              ma.Print("ma");
              duplicate.Print("ma2");
              square.Print("sq");
              square2.Print("sq2");
            }
            CAPTURE(i,j,square(i,j), square2(i,j), diff);
            REQUIRE((diff == Catch::Approx(0).margin(tol)));
        }
    }
    if(!check){
        std::cout<<typeid(matx).name()<<std::endl;
        std::cout << __PRETTY_FUNCTION__ << "failed \n";
    }

    REQUIRE(check);

}

template <class matx, class TVar>
void TestingTransposeWithAutoFill(int rows, int cols, SymProp sp) {
    int i, j;

    matx ma;
    ma.AutoFill(rows, cols, sp);

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
void TestingMultAdd(int dim, SymProp sp, DecomposeType dec) {
    const auto str_sp = SymPropName(sp);
    
    int i, j;

    matx ma;
    ma.AutoFill(dim, dim, sp);

    TPZFMatrix<TVar> cpy(ma), cpy2(ma);
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
        else if constexpr (std::is_same_v<RTVar,long double>) return (RTVar)10;
        else return (RTVar)1;
    }();
    /// Checking whether the res matrix is identical to m1 matrix
    bool check = true;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            TVar zval = z(i, j);
            if (!IsZero(zval/tol)) {
                CAPTURE(dim,str_sp,dec);
                CAPTURE(zval,tol);
                std::cout << "i " << i << " j " << j << " zval " << zval << std::endl;
                if(check) {
                    cpy2.Print("FullMat = ",std::cout,EMathematicaInput);
                    ma.Print("SparseMat = ",std::cout,EMathematicaInput);
                    inv.Print("Inv = ",std::cout,EMathematicaInput);
                    z.Print("zmat = ",std::cout,EMathematicaInput);
                }
                check = false;
            }
        }
    }
    REQUIRE(check);
}

template <class TVar>
void TestingAddContribution(int nrows, int ncols, int ntype)
{
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
              return (RTVar)100;
            else if constexpr (std::is_same_v<RTVar, long double>)
              return (RTVar)10;
            else
              return (RTVar)1;
          }();

          bool check = true;

          for (int i = 0; i < nrows; i++)
          {
              for (int j = 0; j < ncols; j++)
              {
                  TVar diff = C1(i, j) - C2(i, j);
                  if (!IsZero(diff / tol))
                  {
                      CAPTURE(nrows, ncols);
                      CAPTURE(C1(i, j), C2(i, j));
                      std::cout << "i " << i << " j " << j << " C1 " << C1(i, j) << " C2 " << C2(i, j) << std::endl;
                      if (check)
                      {
                        A.Print("A = ", std::cout, EMathematicaInput);
                        B.Print("B = ", std::cout, EMathematicaInput);
                        BT.Print("BT = ", std::cout, EMathematicaInput);
                      }
                      check = false;
                  }
              }
          }

          REQUIRE(check);
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
}

template <class TVar>
void TestingMatRedSolver(int dim0, int dim, int ntype)
{
  TPZMatRed<TVar, TPZFMatrix<TVar>> matred;
  matred.Redim(dim, dim0);

  TPZAutoPointer<TPZMatrix<TVar>> k00 = new TPZFMatrix<TVar>(dim0, dim0, 0.);
  k00->SetSymmetry(SymProp::Sym);
  matred.SetK00(k00);

  matred.K00()->AutoFill(dim0, dim0, SymProp::Sym);
  matred.K11().AutoFill(dim-dim0, dim-dim0, SymProp::Sym);
  matred.K01().AutoFill(dim0, dim-dim0, SymProp::Sym);

  // matred(0, 0) = (TVar)(3);
  // matred(0, 1) = (TVar)(1);
  // matred(0, 2) = (TVar)(4);
  // matred(0, 3) = (TVar)(5);
  // matred(1, 0) = (TVar)(1);
  // matred(1, 1) = (TVar)(2);
  // matred(1, 2) = (TVar)(5);
  // matred(1, 3) = (TVar)(6);
  // matred(2, 2) = (TVar)(7);
  // matred(2, 3) = (TVar)(9);
  // matred(3, 3) = (TVar)(8);

  matred.SimetrizeMatRed();

  TPZFMatrix<TVar> fmatrix(dim, dim, 0.);
  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      fmatrix(i, j) = matred(i, j);
    }
  }

  TPZFMatrix<TVar> F(dim, 1);
  F.AutoFill(dim,1,SymProp::NonSym);
  // F(0, 0) = (TVar)(1.);
  // F(1, 0) = (TVar)(2.);
  // F(2, 0) = (TVar)(3.);
  // F(3, 0) = (TVar)(4.);

  matred.SetF(F);
  matred.SetReduced();

  TPZFMatrix<TVar> u(dim, 1, 0.);
  matred.SolveDirect(u, ELU);

  TPZFMatrix<TVar> u2(dim, 1, 0.);
  u2 = F;
  fmatrix.SolveDirect(u2, ELU);

  constexpr RTVar tol = []()
  {
    if constexpr (std::is_same_v<RTVar, float>)
      return (RTVar)100;
    else if constexpr (std::is_same_v<RTVar, long double>)
      return (RTVar)10;
    else
      return (RTVar)1;
  }();

  bool check = true;

  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < u.Cols(); j++)
    {
      TVar diff = u(i, j) - u2(i, j);
      if (!IsZero(diff / tol))
      {
        CAPTURE(dim, u.Cols());
        CAPTURE(u(i, j), u2(i, j));
        std::cout << "i " << i << " j " << j << " u " << u(i, j) << " u2 " << u2(i, j) << std::endl;
        check = false;
      }
    }
  }

  REQUIRE(check);
}

#ifdef PZ_USING_LAPACK

template <class matx, class TVar>
void TestingGeneralisedEigenValuesWithAutoFill(int dim, SymProp sp) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma,mb;
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
void TestingEigenDecompositionAutoFill(int dim, SymProp sp) {
  auto oldPrecision = Catch::StringMaker<RTVar>::precision;
  Catch::StringMaker<RTVar>::precision = std::numeric_limits<RTVar>::max_digits10;
  matx ma;
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
  if (sizeof(RTVar) == 4) {
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
	TMatrix inverse(nrows,nrows);
	TVar detU = 0.;
	U.DeterminantInverse(detU,inverse);
	if(!IsZero(std::abs(detU) - 1.)) 
		{
            U.DeterminantInverse(detU,inverse);
            check = false;
        }
    TMatrix UTU(U.Rows(),U.Cols()); //< U transposed * U
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
    TMatrix VVT(VT.Rows(),VT.Cols()); //< V * V transposed
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
      CAPTURE(res);
      REQUIRE((std::abs(res.GetVal(ieq,0)) == Catch::Approx(0).margin(tol)));
    }
  }

  SECTION("Residual (block mat)"){
    blcp.Residual(du,rhs,resblck);
    for(auto ieq : blockgraph){
      CAPTURE(ieq);
      CAPTURE(res);
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
};
