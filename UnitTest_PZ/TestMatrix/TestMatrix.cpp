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

#ifdef USING_BOOST

#define BOOST_TEST_DYN_LINK
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
		if(!fabs(matr.GetVal(i,i)) > sum) return 0;
	}
	return 1;
}


#ifdef USING_BOOST

/**
* @brief Tests the Inverse method of the matrix to any matrix types. It uses the AutoFill method to create a square matrix with 
* @param dim Dimension of the square matrix to be build.
* @note Process: build square matrix with randomic values, compute its inverse and the inverse of the inverse.
* Then, checks whether the first and last matrices are identical.
*/
/** Thanks it:
 * NOTE1: I knowed that Inverse function change data of the current matrix. It stores some decomposition (LU or Cholesky ...)
 * NOTE2: If a matrix is decomposed, all times when you call Inverse, it isn't calculated because the matrix has decomposed flag as TRUE 
 * NOTE3: The first function is for square matrices and the second function is for rectangular matrices.
 */
template <class matx>
void TestingInverseWithAutoFill(int dim) {
	int i, j;
	
	matx ma(dim,dim);
	ma.AutoFill();
	BOOST_CHECK_EQUAL(TestingGeneratingDiagonalDominantMatrix<matx>(ma),1);
	// Making ma copy because ma is modified by Inverse method (it's decomposed)
	matx cpy(ma);
	TPZFMatrix<REAL> inv(dim,dim);
	TPZFMatrix<REAL> res(inv);
	// getting inverse twice
	cpy.Inverse(inv);
	inv.Inverse(res);
	
	/// Checking whether the res matrix is identical to m1 matrix
	for(i=0;i<dim;i++)
		for(j=0;j<dim;j++)
			BOOST_CHECK(IsZero(ma.GetVal(i,j) - res.GetVal(i,j)));
}
template <class matx>
void TestingInverseWithAutoFill(int rows,int bnd) {
	int i, j;
//	int columns = rows;
	
	matx ma(rows,bnd);
	ma.AutoFill();
	BOOST_CHECK_EQUAL(TestingGeneratingDiagonalDominantMatrix<matx>(ma),1);
	// Making ma copy because ma is modified by Inverse method (it's decomposed)
	matx cpy(ma);
	TPZFMatrix<REAL> inv(rows,rows);
	TPZFMatrix<REAL> res(rows,rows);
	// getting inverse twice
	cpy.Inverse(inv);
	inv.Inverse(res);
	
	/// Checking whether the res matrix is identical to m1 matrix
	for(i=0;i<rows;i++)
		for(j=0;j<rows;j++)
			BOOST_CHECK(IsZero(ma.GetVal(i,j) - res.GetVal(i,j)));
}

/**
* @brief Tests the operator * of the matrix, using AutoFill to build a square matrix of dimension dim (user defined)
* @param dim Dimension of the square matrix to be build.
* @note Process: build square matrix with randomic values, compute the product using * operator and then checks whether it is a identity matrix.
*/
template <class matx>
void TestingMultiplyOperatorWithAutoFill(int dim) {
	int i, j;
	REAL val;
	// ma times inv must to be a identity matrix
	matx ma(dim,dim);
	ma.AutoFill();

	matx inv(dim,dim);
	matx macpy(ma);

	macpy.Inverse(inv);
	macpy = ma*inv;
	// Checking whether result matrix is the identity matrix
	for(i=0;i<dim;i++) {
		for(j=0;j<dim;j++) {
			val = fabs(macpy.GetVal(i,j));
			if(i==j) BOOST_CHECK(IsZero(val - 1.0));
			else BOOST_CHECK(IsZero(val));
		}
	}
}

/**
* @brief Tests the Multiply method of the matrix, using AutoFill to build a square matrix of dimension dim (user defined)
* @param dim Dimension of the square matrix to be build.
* @note Process: build square matrix with randomic values, compute its inverse then calculate the product of the first and inverse.
*/
/**
 * Check if the result is a identity matrix.
*/
template <class matx>
void TestingMultiplyWithAutoFill(int dim) {
	int i, j;
	REAL val;
	// ma times inv must to be a identity matrix
	matx ma(dim,dim);
	ma.AutoFill();

	matx inv(dim,dim);
	matx res(ma);

	res.Inverse(inv);
	ma.Multiply(inv,res);
	// Checking whether result matrix is the identity matrix
	for(i=0;i<dim;i++) {
		for(j=0;j<dim;j++) {
			val = fabs(res.GetVal(i,j));
			if(i==j) BOOST_CHECK(IsZero(val - 1.0));
			else BOOST_CHECK(IsZero(val));
		}
	}
}

/**
* @brief Tests the Transpose method of the matrix, using AutoFill to build a square matrix of dimension dim (user defined)
* @param dim Dimension of the square matrix to be build.
* @note Process: build square matrix with randomic values, compute its transpose and transpose again then compare the first and last matrices.
*/
template <class matx>
void TestingTransposeWithAutoFill(int rows,int cols) {
	int i, j;
	
	matx ma(rows,cols);
	ma.AutoFill();

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

BOOST_AUTO_TEST_SUITE(matrix_tests)

BOOST_AUTO_TEST_CASE(diagonaldominant_tests) {
	int dim, i;
	
	// Unit Test for full matrix
	for(dim=3;dim<100;dim+=7) {
		TPZFMatrix<REAL> ma(dim,dim);
		ma.AutoFill();
		BOOST_CHECK_EQUAL(TestingGeneratingDiagonalDominantMatrix<TPZFMatrix<REAL> >(ma),1);
	}
	
	// Unit Test for block diagonal matrix
	TPZVec<int> blocks(13);
	for(i=0;i<13;i++)
		blocks[i] = 15+(i%4);
	TPZBlockDiagonal<REAL> mabd(blocks);
	mabd.AutoFill();
	BOOST_CHECK_EQUAL(TestingGeneratingDiagonalDominantMatrix<TPZBlockDiagonal<REAL> >(mabd),1);
	
	// Unit Test No Symmetric Banded matrix
	TPZFBMatrix<REAL> mafb(17,5);
	mafb.AutoFill();
	BOOST_CHECK_EQUAL(TestingGeneratingDiagonalDominantMatrix<TPZFBMatrix<REAL> >(mafb),1);
}

BOOST_AUTO_TEST_CASE(inverse_tests)
{
	int dim;
	for(dim = 9;dim < 100; dim+=5) {
		TestingInverseWithAutoFill<TPZFMatrix<REAL> >(dim);
	//	TestingInverseWithAutoFill<TPZBlockDiagonal>(dim);
		TestingInverseWithAutoFill<TPZFBMatrix<REAL> >(dim,7);
		TestingInverseWithAutoFill<TPZSpMatrix<REAL> >(dim,dim);
		TestingInverseWithAutoFill<TPZFNMatrix<9,REAL> >(dim,dim);
	//	TestingInverseWithAutoFill<TPZSBMatrix>(dim);
	//	TestingInverseWithAutoFill<TPZSFMatrix>(dim);
	}
}

BOOST_AUTO_TEST_CASE(multiply_tests)
{
	int dim;
	for(dim = 3;dim < 100; dim+=5) {
		TestingMultiplyWithAutoFill<TPZFMatrix<REAL> >(dim);
	}
}
BOOST_AUTO_TEST_CASE(multiplyoperator_tests)
{
	int dim;
	for(dim = 3;dim < 100; dim+=5) {
		TestingMultiplyOperatorWithAutoFill<TPZFMatrix<REAL> >(dim);
	}
}
BOOST_AUTO_TEST_CASE(transpose_tests)
{
	int rows, cols;
	for(rows = 3;rows < 4; rows+=5) {
		for(cols = 3;cols < 100;cols+=5) {
			TestingTransposeWithAutoFill<TPZFMatrix<REAL> >(rows,cols);
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()

//BOOST_AUTO_TEST_CASE(nonsingular_test)
//{
//	TestingDiagonalDominant<TPZFBMatrix>();
//}

#endif