/**
 * @file
 * @brief Contains the implementation of the TPZSparseBlockDiagonalStructMatrix methods. 
 */

#include "tpzsparseblockdiagonalstructmatrix.h"
#include "pzcmesh.h"


TPZSparseBlockDiagonalStructMatrix::TPZSparseBlockDiagonalStructMatrix(TPZCompMesh *mesh) :
TPZStructMatrix(mesh)
{
}

TPZSparseBlockDiagonalStructMatrix::TPZSparseBlockDiagonalStructMatrix(TPZAutoPointer<TPZCompMesh>mesh) :
TPZStructMatrix(mesh)
{
}

TPZStructMatrix* TPZSparseBlockDiagonalStructMatrix::Clone()
{
    DebugStop();
	return nullptr;
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
	DebugStop();
	return nullptr;
}

/*!
 \fn TPZSparseBlockDiagonalStructMatrix::NumColors()
 */
int TPZSparseBlockDiagonalStructMatrix::NumColors()
{
    return 0;
}

int TPZSparseBlockDiagonalStructMatrix::ClassId() const{
	return Hash("TPZSparseBlockDiagonalStructMatrix") ^
		TPZStructMatrix::ClassId() << 1;
}

void TPZSparseBlockDiagonalStructMatrix::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPZStructMatrixOR::Read(buf,context);
}

void TPZSparseBlockDiagonalStructMatrix::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPZStructMatrixOR::Write(buf,withclassid);
}