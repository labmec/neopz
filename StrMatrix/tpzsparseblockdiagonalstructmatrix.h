/**
 * @file
 * @brief Contains the TPZSparseBlockDiagonalStructMatrix class which builds a sparse block diagonal preconditioner.
 */

#ifndef TPZSPARSEBLOCKDIAGONALSTRUCTMATRIX_H
#define TPZSPARSEBLOCKDIAGONALSTRUCTMATRIX_H

#include "TPZStructMatrix.h"
#include "pzstrmatrixor.h"

/**
 * @brief It will build a sparse block diagonal preconditioner with a structure determined by the parameters passed to it. \ref structural "Structural Matrix"
 * @ingroup structural
 * @author Philippe R. B. Devloo
 */
class TPZSparseBlockDiagonalStructMatrix : public TPZStructMatrix, public TPZStructMatrixOR
{
public:
	/** @brief Constructor for computational mesh */
    TPZSparseBlockDiagonalStructMatrix(TPZCompMesh *mesh);

    TPZSparseBlockDiagonalStructMatrix(TPZAutoPointer<TPZCompMesh>mesh);
	
	TPZMatrix<STATE> * Create() override;
	
    TPZStructMatrix* Clone() override;
    
    int NumColors();

    //@{
    //!Read and Write methods
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
    
private:
    TPZSparseBlockDiagonalStructMatrix() = default;

    friend TPZPersistenceManager;
};

#endif
