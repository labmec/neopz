/**
 * @file
 * @brief Contains Unit Tests for methods of the matrices classes.
 */

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzblockdiag.h"
#include "pzbndmat.h"
#include "pzespmat.h"
#include "pzsbndmat.h"
#include "pzsfulmat.h"
#include "pzskylnsymmat.h"
#include "pzskylmat.h"
#include "pzysmp.h"
#include "pzsysmp.h"

#ifdef PZ_USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz matrix tests

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#endif


namespace testmatrix{
/**
 * @brief Tests whether a matrix matr is diagonally dominant. For fixed i, checks the condition |Aii| > sum_j(|Aij|) on j!=i.
 * @param matr Matrix to check, it can to be non square matrix.
 * @note matx is a class of the type matrix matr.
 */
template <class matx>
int TestingGeneratingDiagonalDominantMatrix(matx &matr);

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
    for (auto dim = 9; dim < 10; dim += 5) {
#ifdef PZ_USING_MKL
      if constexpr (std::is_same<TVar,double>::value){
        TestingInverseWithAutoFill<TPZSYsmpMatrix<TVar>,TVar>(dim, 1, ECholesky);
        TestingInverseWithAutoFill<TPZSYsmpMatrix<TVar>,TVar>(dim, 1, ELDLt);
      }
#endif
      TestingInverseWithAutoFill<TPZSBMatrix<TVar>,TVar>(dim, 1, ELDLt);
      TestingInverseWithAutoFill<TPZSBMatrix<TVar>,TVar>(dim, 1,ECholesky);
      TestingInverseWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, 0, ELU);
      TestingInverseWithAutoFill<TPZFMatrix<TVar>,TVar>(dim, 1,ECholesky);
      TestingInverseWithAutoFill<TPZFMatrix<TVar>,TVar>(dim, 0, ELU);
      TestingInverseWithAutoFill<TPZFMatrix<TVar>,TVar>(dim, 1, ELDLt);
      TestingInverseWithAutoFill<TPZFBMatrix<TVar>,TVar>(dim, 0, ELU);
      TestingInverseWithAutoFill<TPZBlockDiagonal<TVar>,TVar>(dim, 0,ELU);
      TestingInverseWithAutoFill<TPZSkylMatrix<TVar>,TVar>(dim, 1,ELDLt);
      TestingInverseWithAutoFill<TPZFNMatrix<9,TVar>,TVar>(dim, 0, ELU);
      TestingInverseWithAutoFill<TPZSFMatrix<TVar>,TVar>(dim, 1, ELDLt);
      TestingInverseWithAutoFill<TPZSkylNSymMatrix<TVar>,TVar>(dim, 0,ELU);
    }
  }

    template<class TVar>
    void MultiplyTranspose(){
      for (auto dim = 5; dim < 100; dim += 500) {
        TestingTransposeMultiply<TPZFYsmpMatrix<TVar>, TVar>(10, dim, 0);
        TestingTransposeMultiply<TPZSYsmpMatrix<TVar>, TVar>(dim, dim, 1);
        TestingTransposeMultiply<TPZFMatrix<TVar>, TVar>(10, dim, 0);
        TestingTransposeMultiply<TPZSkylNSymMatrix<TVar>, TVar>(dim, dim,
                                                                    0);
      }
    }

    template<class TVar>
    void Multiply(){
      for (auto dim = 20; dim < 100; dim += 500) {
        TestingMultiplyWithAutoFill<TPZFYsmpMatrix<TVar>, TVar>(dim, 0);
        TestingMultiplyWithAutoFill<TPZSkylNSymMatrix<TVar>, TVar>(dim, 0);
        TestingMultiplyWithAutoFill<TPZFMatrix<TVar>, TVar>(dim, 0);
      }
    }

    template<class TVar>
    void DiagonalDominant(){
        // Unit Test for full matrix
        for (int dim = 3; dim < 100; dim += 7) {
          TPZFMatrix<TVar> ma;
          ma.AutoFill(dim, dim, 0);
          BOOST_CHECK_EQUAL(
              TestingGeneratingDiagonalDominantMatrix<TPZFMatrix<TVar>>(ma), 1);
        }

        // Unit Test for block diagonal matrix
        TPZVec<int> blocks(13);
        for (int i = 0; i < 13; i++)
          blocks[i] = 15 + (i % 4);
        TPZBlockDiagonal<TVar> mabd(blocks);
        mabd.AutoFill(50, 50, 0);
        BOOST_CHECK_EQUAL(
            TestingGeneratingDiagonalDominantMatrix<TPZBlockDiagonal<TVar>>(
                mabd),
            1);

        // Unit Test No Symmetric Banded matrix
        TPZFBMatrix<TVar> mafb(17, 5);
        mafb.AutoFill(17, 17, 0);
        BOOST_CHECK_EQUAL(
            TestingGeneratingDiagonalDominantMatrix<TPZFBMatrix<TVar>>(mafb),
            1);
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
        TestingMultAdd<TPZFMatrix<float>, float>(dim, 1, ECholesky);
        TestingMultAdd<TPZSkylMatrix<float>, float>(dim, 1, ECholesky);
        TestingMultAdd<TPZSkylNSymMatrix<float>, float>(dim, 1, ECholesky);
      }
    }
#ifdef PZ_USING_LAPACK
    template <class TVar> void GeneralisedEigenvaluesAutoFill() {
    for (int dim = 5; dim < 6; dim += 10) {
      TestingGeneralisedEigenValuesWithAutoFill<TPZFMatrix<TVar>, TVar>(
          dim, 1);
      TestingGeneralisedEigenValuesWithAutoFill<TPZSBMatrix<TVar>, TVar>(
          dim, 1);
    }
  }
    template<class TVar>
    void EigenDecompositionAutoFill(){
      for (int dim = 5; dim < 6; dim += 10) {
        TestingEigenDecompositionAutoFill<TPZSBMatrix<TVar>, TVar>(dim, 1);
        TestingEigenDecompositionAutoFill<TPZFMatrix<TVar>, TVar>(dim, 1);
        TestingEigenDecompositionAutoFill<TPZFMatrix<TVar>, TVar>(dim, 0);
      }
    }
#endif
};

BOOST_AUTO_TEST_SUITE(matrix_tests)

BOOST_AUTO_TEST_CASE(inverse_tests_real) {
    // testmatrix::Inverse<float>();//FAILING
    testmatrix::Inverse<double>();
    testmatrix::Inverse<long double>();
}
BOOST_AUTO_TEST_CASE(inverse_tests_complex) {
    // testmatrix::Inverse<std::complex<float>>();
    // testmatrix::Inverse<std::complex<double>>();
    // testmatrix::Inverse<std::complex<long double>>();
    //FAILING
}

BOOST_AUTO_TEST_CASE(multiplytranspose_tests_real) {
    testmatrix::MultiplyTranspose<float>();
    testmatrix::MultiplyTranspose<double>();
    testmatrix::MultiplyTranspose<long double>();
}
BOOST_AUTO_TEST_CASE(multiplytranspose_tests_complex){
    testmatrix::MultiplyTranspose<std::complex<float>>();
    testmatrix::MultiplyTranspose<std::complex<double>>();
    testmatrix::MultiplyTranspose<std::complex<long double>>();
}

BOOST_AUTO_TEST_CASE(multiply_tests_real) {
    // testmatrix::Multiply<float>();//FAILING
    testmatrix::Multiply<double>();
    testmatrix::Multiply<long double>();
}
BOOST_AUTO_TEST_CASE(multiply_tests_complex){
    testmatrix::Multiply<std::complex<float>>();
    testmatrix::Multiply<std::complex<double>>();
    testmatrix::Multiply<std::complex<long double>>();
}

BOOST_AUTO_TEST_CASE(diagonaldominant_tests) {
    testmatrix::DiagonalDominant<float>();
    testmatrix::DiagonalDominant<double>();
    testmatrix::DiagonalDominant<long double>();
    testmatrix::DiagonalDominant<std::complex<float>>();
    testmatrix::DiagonalDominant<std::complex<double>>();
    testmatrix::DiagonalDominant<std::complex<long double>>();
}


BOOST_AUTO_TEST_CASE(multiplyoperator_tests) {
    testmatrix::MultiplyOperatorWithAutoFill<float>();
    testmatrix::MultiplyOperatorWithAutoFill<double>();
    testmatrix::MultiplyOperatorWithAutoFill<long double >();
    testmatrix::MultiplyOperatorWithAutoFill<std::complex<float>>();
    testmatrix::MultiplyOperatorWithAutoFill<std::complex<double>>();
    testmatrix::MultiplyOperatorWithAutoFill<std::complex<long double >>();
}

BOOST_AUTO_TEST_CASE(transpose_tests) {
    testmatrix::TransposeWithAutoFill<float>();
    testmatrix::TransposeWithAutoFill<double>();
    testmatrix::TransposeWithAutoFill<long double>();
    testmatrix::TransposeWithAutoFill<std::complex<float>>();
    testmatrix::TransposeWithAutoFill<std::complex<double>>();
    testmatrix::TransposeWithAutoFill<std::complex<long double>>();
}
BOOST_AUTO_TEST_CASE(multadd_tests) {
    testmatrix::TestMultAdd<float>();
    testmatrix::TestMultAdd<double>();
    testmatrix::TestMultAdd<long double>();
    testmatrix::TestMultAdd<std::complex<float>>();
    testmatrix::TestMultAdd<std::complex<double>>();
    testmatrix::TestMultAdd<std::complex<long double>>();
}

#ifdef PZ_USING_LAPACK

BOOST_AUTO_TEST_CASE(eigenvalue_tests) {
    
    testmatrix::EigenDecompositionAutoFill<float>();
    testmatrix::EigenDecompositionAutoFill<double>();
    //THERE IS NO LAPACK INTERFACE FOR LONG DOUBLE IN OUR CODE
    // testmatrix::EigenDecompositionAutoFill<long double>();
    testmatrix::EigenDecompositionAutoFill<std::complex<float>>();
    testmatrix::EigenDecompositionAutoFill<std::complex<double>>();
    //THERE IS NO LAPACK INTERFACE FOR LONG DOUBLE IN OUR CODE
    // testmatrix::EigenDecompositionAutoFill<std::complex<long double>>();
}

BOOST_AUTO_TEST_CASE(generalized_eigenvalue_tests) {
    testmatrix::GeneralisedEigenvaluesAutoFill<float>();
    testmatrix::GeneralisedEigenvaluesAutoFill<double>();
    //THERE IS NO LAPACK INTERFACE FOR LONG DOUBLE IN OUR CODE
    //testmatrix::GeneralisedEigenvaluesAutoFill<long double>();
    testmatrix::GeneralisedEigenvaluesAutoFill<std::complex<float>>();
    testmatrix::GeneralisedEigenvaluesAutoFill<std::complex<double>>();
    //THERE IS NO LAPACK INTERFACE FOR LONG DOUBLE IN OUR CODE
    //testmatrix::GeneralisedEigenvaluesAutoFill<std::complex<long double>>();
}
#endif

BOOST_AUTO_TEST_SUITE_END()

//BOOST_AUTO_TEST_CASE(nonsingular_test)
//{
//	TestingDiagonalDominant<TPZFBMatrix>();
//}



namespace testmatrix{

template <class matx>
int TestingGeneratingDiagonalDominantMatrix(matx &matr) {
    REAL sum;
    for (int i = 0; i < matr.Rows(); i++) {
        sum = 0.0;
        for (int j = 0; j < matr.Cols(); j++) {
            if (i != j)
                sum += fabs(matr.GetVal(i, j));
        }
        if (!(fabs(matr.GetVal(i, i)) > sum)) {
            std::cout << "line i " << i << " failed\n";
            matr.Print("matrix = ", std::cout, EMathematicaInput);
            return 0;
        }
    }
    return 1;
}


#ifdef PZ_USING_BOOST

template <class matx, class TVar>
void TestingInverseWithAutoFill(int dim, int symmetric, DecomposeType dec) {
  int i, j;

  matx ma;
  ma.AutoFill(dim, dim, symmetric);
  const bool firsttest = TestingGeneratingDiagonalDominantMatrix<matx>(ma);
  if(!firsttest){
      std::cout << __PRETTY_FUNCTION__;
      std::cout << " failed to generate matrix\n";
      std::cout << std::flush;
  }
  //    ma.Print("skyl =",std::cout,EMathematicaInput);
  BOOST_CHECK_EQUAL(firsttest, 1);
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
      bool loccheck = IsZero(diff / 10.);
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
        cpma.Print("Matrix = ", std::cout, EMathematicaInput);
        invkeep.Print("Inv = ", std::cout, EMathematicaInput);
    }
    BOOST_CHECK(check);
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
    BOOST_CHECK(check);
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
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            if (!IsZero(fabs(square(i, j) - square2(i, j)))) {
                check = false;
            }
        }
    }
    BOOST_CHECK(check);

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
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < col; j++) {
            if (!IsZero(fabs(square(i, j) - square2(i, j)))) {
                check = false;
            }
        }
    }
    if (!check) {
        std::cout << __PRETTY_FUNCTION__ << "failed \n";
        square2.Print("Full");
        square.Print("Matrix");
        square2 -= square;
        square2.Print("Diff");
    }
    BOOST_CHECK(check);

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
            BOOST_CHECK(IsZero(ma.GetVal(i, j) - matransptransp.GetVal(i, j)));
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
    BOOST_CHECK(check);
}

#ifdef PZ_USING_LAPACK

template <class matx, class TVar>
void TestingGeneralisedEigenValuesWithAutoFill(int dim, int symmetric) {

    matx ma, mb;
    ma.AutoFill(dim, dim, symmetric);
    mb.AutoFill(dim, dim, symmetric);

    BOOST_CHECK_EQUAL(TestingGeneratingDiagonalDominantMatrix<matx>(ma), 1);
    BOOST_CHECK_EQUAL(TestingGeneratingDiagonalDominantMatrix<matx>(mb), 1);

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
    BOOST_CHECK(check);
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
    BOOST_CHECK(check);
}


#endif
};
#endif
