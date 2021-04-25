/**
 * @file
 * @brief Contains the TPZSparseBlockDiagonalStructMatrix class which builds a sparse block diagonal preconditioner.
 */

#ifndef TPZSPARSEBLOCKDIAGONALSTRUCTMATRIX_H
#define TPZSPARSEBLOCKDIAGONALSTRUCTMATRIX_H

#include "TPZStructMatrixT.h"
#include "pzstrmatrixor.h"

/**
 * @brief It will build a sparse block diagonal preconditioner with a structure determined by the parameters passed to it. \ref structural "Structural Matrix"
 * @ingroup structural
 * @author Philippe R. B. Devloo
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSparseBlockDiagonalStructMatrix : public TPZStructMatrixT<TVar>,
                                  public TPar
{
public:
    using TPZStructMatrixT<TVar>::TPZStructMatrixT;
	
	TPZMatrix<TVar> * Create() override;
	
    TPZStructMatrix* Clone() override;
    
    int NumColors();

    //@{
    //!Read and Write methods
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}

    friend TPZPersistenceManager;
};

#endif
