/**
 * @file
 * @brief Contains the implementation of the TPZSparseBlockDiagonalStructMatrix methods. 
 */
//
// C++ Implementation: tpzsparseblockdiagonalstructmatrix
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpzsparseblockdiagonalstructmatrix.h"

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

TPZMatrix<REAL> * TPZSparseBlockDiagonalStructMatrix::Create()
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
