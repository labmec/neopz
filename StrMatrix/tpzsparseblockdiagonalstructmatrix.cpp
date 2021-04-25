/**
 * @file
 * @brief Contains the implementation of the TPZSparseBlockDiagonalStructMatrix methods. 
 */

#include "tpzsparseblockdiagonalstructmatrix.h"
#include "pzcmesh.h"

template<class TVar, class TPar>
TPZStructMatrix* TPZSparseBlockDiagonalStructMatrix<TVar,TPar>::Clone()
{
    DebugStop();
	return nullptr;
}

template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSparseBlockDiagonalStructMatrix<TVar,TPar>::Create()
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
template<class TVar, class TPar>
 \fn TPZSparseBlockDiagonalStructMatrix<TVar,TPar>::NumColors()
 */
template<class TVar, class TPar>
int TPZSparseBlockDiagonalStructMatrix<TVar,TPar>::NumColors()
{
    return 0;
}

template<class TVar, class TPar>
int TPZSparseBlockDiagonalStructMatrix<TVar,TPar>::ClassId() const{
	return Hash("TPZSparseBlockDiagonalStructMatrix") ^
		TPZStructMatrix::ClassId() << 1 ^
        TPar::ClassId() << 2;
}

template<class TVar, class TPar>
void TPZSparseBlockDiagonalStructMatrix<TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPar::Read(buf,context);
}

template<class TVar, class TPar>
void TPZSparseBlockDiagonalStructMatrix<TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPar::Write(buf,withclassid);
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZSparseBlockDiagonalStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZSparseBlockDiagonalStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZSparseBlockDiagonalStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;
