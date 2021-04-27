/**
 * @file
 * @brief Contains the TPZSFStructMatrix class which implements Full Symmetric Structural Matrices.
 */

#ifndef TPZSFSTRUCTMATRIX_H
#define TPZSFSTRUCTMATRIX_H

#include "TPZStructMatrixT.h"
#include "pzstrmatrixor.h"

/**
 * @brief Implements a full symmetric structural matrix
 * using TPZSFMatrix as a storage format.
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSFStructMatrix : public TPZStructMatrixT<TVar>,
                                  public TPar
{
public:    	
    using TPZStructMatrixT<TVar>::TPZStructMatrixT;
    TPZMatrix<TVar>* Create() override;
	
    TPZStructMatrix * Clone() override;
    //@{
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
    
    friend TPZPersistenceManager;
};

#endif //TPZFSTRUCTMATRIX_H
