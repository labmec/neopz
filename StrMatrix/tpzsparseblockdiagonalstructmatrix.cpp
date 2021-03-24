/**
 * @file
 * @brief Contains the implementation of the TPZSparseBlockDiagonalStructMatrix methods. 
 */

#include "tpzsparseblockdiagonalstructmatrix.h"

TPZSparseBlockDiagonalStructMatrix::TPZSparseBlockDiagonalStructMatrix() :
TPZStructMatrix()
{
}

TPZSparseBlockDiagonalStructMatrix::TPZSparseBlockDiagonalStructMatrix(TPZCompMesh *mesh) :
TPZStructMatrix(mesh)
{
}

TPZSparseBlockDiagonalStructMatrix::~TPZSparseBlockDiagonalStructMatrix()
{
}

TPZStructMatrix* TPZSparseBlockDiagonalStructMatrix::Clone()
{
    return TPZStructMatrix::Clone();
}

TPZMatrix<STATE> * TPZSparseBlockDiagonalStructMatrix::Create()
{
	//extract the structure of the matrix from the mesh
	//create an object of type SparseBlockDiagonal
	
	// now, the values of the matrix cannot be initialized till the 
	// global matrix has been assembled
	// what about colors and overlapping techniques?
	
	// we need to create a structure which contains a vector of these
	// somehow creating solvers and matrix objects need to interact
	// a non overlapping preconditioner is a sequence solver
	// an overlapping preconditioner is a step solver
	// this preconditioner type should be created by the analysis class
	// it is too complicated to require the user to coordinate this within
	// the current class structure....
	return 0;
}

/*!
 \fn TPZSparseBlockDiagonalStructMatrix::NumColors()
 */
int TPZSparseBlockDiagonalStructMatrix::NumColors()
{
    return 0;
}
