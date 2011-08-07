/**
 * @file
 * @brief Contains the TPZBlockDiagonalStructMatrix class which implements Block Diagonal Structural Matrices.
 */

#ifndef TPZBLOCKDIAGONALSTRUCTMATRIX_H
#define TPZBLOCKDIAGONALSTRUCTMATRIX_H

#include "pzstrmatrix.h"

class TPZCompMesh;
class TPZFMatrix;
class TPZMatrix;
class TPZBlockDiagonal;
template <class T>
class TPZVec;

/**
 * @brief Implements Block Diagonal Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZBlockDiagonalStructMatrix : public TPZStructMatrix {
public:    
	
	enum MBlockStructure {ENodeBased, EVertexBased, EElementBased};
	
	TPZBlockDiagonalStructMatrix(TPZCompMesh *);
	
	~TPZBlockDiagonalStructMatrix();
	
	TPZBlockDiagonalStructMatrix(const TPZBlockDiagonalStructMatrix &copy) : TPZStructMatrix(copy),
	fBlockStructure(copy.fBlockStructure),fOverlap(copy.fOverlap)
	{
	}
	
	// @brief Creates a sparse blockdiagonal matrix, overlapping should be assumed
	virtual TPZMatrix * Create();
    
	virtual TPZMatrix * CreateAssemble(TPZFMatrix &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	virtual TPZStructMatrix * Clone();    
	
public:
	
	void AssembleBlockDiagonal(TPZBlockDiagonal & block);
private:
	
    void BlockSizes(TPZVec < int > & blocksizes);
    
    MBlockStructure fBlockStructure;
    int fOverlap;
    
    
};

#endif //TPZBLOCKDIAGONALSTRUCTMATRIX_H
