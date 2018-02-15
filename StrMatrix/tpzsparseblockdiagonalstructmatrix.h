/**
 * @file
 * @brief Contains the TPZSparseBlockDiagonalStructMatrix class which builds a sparse block diagonal preconditioner.
 */

#ifndef TPZSPARSEBLOCKDIAGONALSTRUCTMATRIX_H
#define TPZSPARSEBLOCKDIAGONALSTRUCTMATRIX_H

#include "pzstrmatrix.h"

/**
 * @brief It will build a sparse block diagonal preconditioner with a structure determined by the parameters passed to it. \ref structural "Structural Matrix"
 * @ingroup structural
 * @author Philippe R. B. Devloo
 */
class TPZSparseBlockDiagonalStructMatrix : public TPZStructMatrix
{
public:
	/** @brief Constructor for computational mesh */
    TPZSparseBlockDiagonalStructMatrix(TPZCompMesh *mesh);
	/** @brief Destructor */
    ~TPZSparseBlockDiagonalStructMatrix();
	
	virtual TPZMatrix<STATE> * Create();
	
    virtual TPZStructMatrix* Clone();
    int NumColors();
    
private:
    TPZSparseBlockDiagonalStructMatrix();

    friend TPZPersistenceManager;
};

#endif
