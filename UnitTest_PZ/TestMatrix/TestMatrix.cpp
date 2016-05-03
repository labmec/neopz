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

#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz matrix tests

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#endif

/** 
* @brief Tests wheter a matrix matr is diagonally dominant. Checks Aii > Soma(Aij) on j!=i, for all i.
* @param matr Matrix to check, it can to be non square matrix.
* @note matx is a class of the type matrix matr.
*/
template <class matx>
int TestingGeneratingDiagonalDominantMatrix(matx &matr) {
	REAL sum;
	for(int i=0;i<matr.Rows();i++) {
		sum = 0.0;
		for(int j=0;j<matr.Cols();j++) {
			if(i!=j)
				sum += fabs(matr.GetVal(i,j));
		}
		if(!(fabs(matr.GetVal(i,i)) > sum))
        {
            std::cout << "line i " << i << " failed\n";
            matr.Print("matrix = " ,std::cout,EMathematicaInput);
            return 0;
        }
	}
	return 1;
}


#ifdef USING_BOOST

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
void TestingInverseWithAutoFill(int dim, int symmetric, DecomposeType dec) {
	int i, j;
	
	matx ma;
	ma.AutoFill(dim,dim,symmetric);
    
//    ma.Print("skyl =",std::cout,EMathematicaInput);
	BOOST_CHECK_EQUAL(TestingGeneratingDiagonalDominantMatrix<matx>(ma),1);
	// Making ma copy because ma is modified by Inverse method (it's decomposed)
	matx cpma(ma);
	TPZFMatrix<TVar> inv(dim,dim), invkeep;
	TPZFMatrix<TVar> res(inv);
	// getting inverse twice
	ma.Inverse(inv,dec);
  invkeep = inv;
	inv.Inverse(res,dec);
  //    ma.Print("skyl2 =",std::cout,EMathematicaInput);
  bool check = true;
	/// Checking whether the res matrix is identical to m1 matrix
	for(i=0;i<dim;i++)
    {
      for(j=0;j<dim;j++)
      {
          TVar diff = cpma.GetVal(i,j) - res.GetVal(i,j);
          bool loccheck = IsZero(TVar(diff/10.));
          if (loccheck == false) {
              std::cout << "diff " << diff << std::endl;
          }
          check &= loccheck;
      }
    }
    if (!check) {
        cpma.Print("Matrix = ",std::cout,EMathematicaInput);
        invkeep.Print("Inv = ",std::cout,EMathematicaInput);
    }
    BOOST_CHECK(check);
}

/**
* @brief Tests the operator * of the matrix, using AutoFill to build a square matrix of dimension dim (user defined)
* @param dim Dimension of the square matrix to be build.
* @param symmetric Whether to build a symmetric matrix
* @note Process: build square matrix with randomic values, compute its square twice using * operator with two different copies of itself and checks whether the result is the same
*/
template <class matx, class TVar>
void TestingMultiplyOperatorWithAutoFill(int dim, int symmetric) {
	// ma times inv must to be a identity matrix
	matx ma;
	ma.AutoFill(dim,dim,symmetric);

	TPZFMatrix<TVar> duplicate(ma);
	TPZFMatrix<TVar> square,square2;

	square2 = duplicate*duplicate;
	square = ma*duplicate;
	// Checking whether both matrices are equal
  bool check = true;
	for(int i=0;i<dim;i++) {
		for(int j=0;j<dim;j++) {
      TVar diff = fabs(square.GetVal(i, j)- square2.GetVal(i, j));
      if (!IsZero(diff)) {
          check = false;
      }
		}
	}
  BOOST_CHECK(check);
}

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
void TestingMultiplyWithAutoFill(int dim, int symmetric) {
	// ma times inv must to be a identity matrix
	matx ma;
	ma.AutoFill(dim,dim,symmetric);

	TPZFMatrix<TVar> duplicate(ma), square,square2;

//    ma.Print("SkylineNS");
//    duplicate.Print("FullMat");
	ma.Multiply(duplicate,square);
  duplicate.Multiply(duplicate, square2);
	// Checking whether result matrix is the identity matrix
    bool check = true;
	for(int i=0;i<dim;i++) {
		for(int j=0;j<dim;j++) {
      if(!IsZero(fabs(square(i,j)-square2(i,j))))
      {
          check = false;
      }
		}
	}
  BOOST_CHECK(check);

}

/**
* @brief Tests the Transpose method of the matrix, using AutoFill to build a matrix of dimension row x cols (user defined)
* @param rows Number of rows
* @param cols Number of columns
* @param symmetric Whether to build a symmetric matrix
* @note Process: build matrix with randomic values, compute its transpose and transpose again then compare the first and last matrices.
*/
template <class matx, class TVar>
void TestingTransposeWithAutoFill(int rows,int cols, int symmetric) {
	int i, j;
	
	matx ma;
	ma.AutoFill(rows,cols,symmetric);

	matx matransp(cols,rows);
	matx matransptransp(ma);

	// getting inverse twice
	matransptransp.Transpose(&matransp);
	matransp.Transpose(&matransptransp);
	
	/// Checking whether the res matrix is identical to m1 matrix
	for(i=0;i<rows;i++)
		for(j=0;j<cols;j++)
			BOOST_CHECK(IsZero(ma.GetVal(i,j) - matransptransp.GetVal(i,j)));
}

/**
 * @brief Tests the MultAdd method of the matrix, that calculates z = beta * y + alpha * A * x , using AutoFill to build a square matrix of dimension dim (user defined)
 * @param dim Dimension of the square matrix to be build.
 * @param symmetric Whether to build a symmetric matrix
 * @param dec Decomposition method to be used (See enum DecomposeType at pzmatrix.h)
 * @note Process: build square matrix with randomic values, compute its inverse and uses MultAdd for calculating I - A*A^-1 and check whether you have any non-zero entries on the result
 */
template <class matx, class TVar>
void TestingMultAdd(int dim, int symmetric, DecomposeType dec) {
    int i, j;
    
    matx ma;
    ma.AutoFill(dim,dim,symmetric);
    
    TPZFMatrix<TVar> cpy(ma);
    TPZFMatrix<TVar> inv(dim,dim);
    TPZFMatrix<TVar> y(ma),z;
    y.Identity();
    // getting inverse twice
    cpy.Inverse(inv,dec);

//    virtual void MultAdd(const TPZFMatrix<TVar> & x,const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z,
//                         const TVar alpha=1., const TVar beta = 0., const int opt = 0) const;


    ma.MultAdd(inv,y,z,1.,-1.);
    /// Checking whether the res matrix is identical to m1 matrix
    bool check = true;
    for(i=0;i<dim;i++)
    {
        for(j=0;j<dim;j++)
        {
            TVar zval = z(i,j);
            if(!IsZero(zval))
            {
                check = false;
            }
        }
    }
    BOOST_CHECK(check);
}

BOOST_AUTO_TEST_SUITE(matrix_tests)

BOOST_AUTO_TEST_CASE(diagonaldominant_tests) {
	int dim, i;
	
	// Unit Test for full matrix
	for(dim=3;dim<100;dim+=7) {
		TPZFMatrix<REAL> ma;
		ma.AutoFill(dim,dim,0);
		BOOST_CHECK_EQUAL(TestingGeneratingDiagonalDominantMatrix<TPZFMatrix<REAL> >(ma),1);
	}
	
	// Unit Test for block diagonal matrix
	TPZVec<int> blocks(13);
	for(i=0;i<13;i++)
		blocks[i] = 15+(i%4);
	TPZBlockDiagonal<REAL> mabd(blocks);
	mabd.AutoFill(50,50,0);
	BOOST_CHECK_EQUAL(TestingGeneratingDiagonalDominantMatrix<TPZBlockDiagonal<REAL> >(mabd),1);
	
	// Unit Test No Symmetric Banded matrix
	TPZFBMatrix<REAL> mafb(17,5);
	mafb.AutoFill(17,17,0);
	BOOST_CHECK_EQUAL(TestingGeneratingDiagonalDominantMatrix<TPZFBMatrix<REAL> >(mafb),1);
}

BOOST_AUTO_TEST_CASE(inverse_tests)
{
	int dim;
	for(dim = 9;dim < 10; dim+=5) {
    TestingInverseWithAutoFill<TPZSBMatrix<float>, float >(dim,1,ELDLt);
    TestingInverseWithAutoFill<TPZSBMatrix<float>, float >(dim,1,ECholesky);
    TestingInverseWithAutoFill<TPZFBMatrix<double>,double >(dim,0,ELU);
    TestingInverseWithAutoFill<TPZFBMatrix<float>,float >(dim,0,ELU);
    TestingInverseWithAutoFill<TPZFMatrix<float>, float>(dim,1,ECholesky);
    TestingInverseWithAutoFill<TPZFMatrix<float>, float>(dim,0,ELU);
    TestingInverseWithAutoFill<TPZFMatrix<double>, double>(dim,1,ECholesky);
    TestingInverseWithAutoFill<TPZFMatrix<float>, float>(dim,1,ELDLt);
    TestingInverseWithAutoFill<TPZFMatrix<double>, double>(dim,1,ELDLt);
    TestingInverseWithAutoFill<TPZBlockDiagonal<float>,float >(dim,0,ELU);
    TestingInverseWithAutoFill<TPZSkylMatrix<float>, float >(dim,1,ELDLt);
    TestingInverseWithAutoFill<TPZFNMatrix<9,float>, float >(dim,0,ELU);
    TestingInverseWithAutoFill<TPZSFMatrix<float>, float >(dim,1,ELDLt);
    TestingInverseWithAutoFill<TPZSkylNSymMatrix<float>, float >(dim,0,ELU);
    
	}
}

BOOST_AUTO_TEST_CASE(multiply_tests)
{
	int dim;
	for(dim = 3;dim < 100; dim+=5) {
		TestingMultiplyWithAutoFill<TPZSkylNSymMatrix<REAL> ,REAL>(dim,0);
		TestingMultiplyWithAutoFill<TPZFMatrix<REAL>, REAL >(dim,0);
	}
}
BOOST_AUTO_TEST_CASE(multiplyoperator_tests)
{
	int dim;
	for(dim = 3;dim < 100; dim+=5) {
		TestingMultiplyOperatorWithAutoFill<TPZFMatrix<float>, float >(dim,0);
    TestingMultiplyOperatorWithAutoFill<TPZFMatrix<double>, double >(dim,0);
	}
}
BOOST_AUTO_TEST_CASE(transpose_tests)
{
	int rows, cols;
	for(rows = 3;rows < 4; rows+=5) {
		for(cols = 3;cols < 100;cols+=5) {
			TestingTransposeWithAutoFill<TPZFMatrix<float>, float >(rows,cols,0);
      TestingTransposeWithAutoFill<TPZFMatrix<double>, double >(rows,cols,0);
		}
	}
}

BOOST_AUTO_TEST_CASE(multadd_tests)
{
    for (int dim = 5; dim < 6; dim += 10) {
      TestingMultAdd<TPZFMatrix<float>, float>(dim, 1, ECholesky);
      TestingMultAdd<TPZFMatrix<double>, double>(dim, 1, ECholesky);
      TestingMultAdd<TPZSkylMatrix<float>, float>(dim, 1, ECholesky);
      TestingMultAdd<TPZSkylMatrix<double>, double>(dim, 1, ECholesky);
      TestingMultAdd<TPZSkylNSymMatrix<float>, float>(dim, 1, ECholesky);
      TestingMultAdd<TPZSkylNSymMatrix<double>, double>(dim, 1, ECholesky);
    }
}
BOOST_AUTO_TEST_SUITE_END()

//BOOST_AUTO_TEST_CASE(nonsingular_test)
//{
//	TestingDiagonalDominant<TPZFBMatrix>();
//}

#endif