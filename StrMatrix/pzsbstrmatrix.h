/**
 * @file
 * @brief Contains the TPZSBandStructMatrix class which implements Symmetric Banded Structural Matrices.
 */

#ifndef TPZSBANDSTRUCTMATRIX_H
#define TPZSBANDSTRUCTMATRIX_H

#include "TPZStructMatrixT.h"
#include "pzstrmatrixor.h"


template<class T>
class TPZStructMatrixOR;
/**
 * @brief Implements a symmetric banded structural matrix using TPZSBMatrix as a storage format.
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSBandStructMatrix : public TPZStructMatrixT<TVar>,
                             public TPar
{
	using TPZStructMatrixT<TVar>::TPZStructMatrixT;
    
    TPZMatrix<TVar> * Create() override;
	
    TPZStructMatrix * Clone() override;
    //@{
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
private :
	
    friend TPZPersistenceManager;
};

#endif //TPZSBANDSTRUCTMATRIX_H
