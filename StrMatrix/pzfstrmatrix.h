/**
 * @file
 * @brief Contains the TPZFStructMatrix class which implements Full Structural Matrices.
 */

#ifndef TPZFSTRUCTMATRIX_H
#define TPZFSTRUCTMATRIX_H

#include "TPZStructMatrixT.h"
#include "pzstrmatrixor.h"

/**
 * @brief Implements a full structural matrix using TPZFMatrix as a storage format.
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZFStructMatrix : public TPZStructMatrixT<TVar>,
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
