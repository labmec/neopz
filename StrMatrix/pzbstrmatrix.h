/**
 * @file
 * @brief Contains the TPZBandStructMatrix class which implements Banded Structural Matrices.
 */

#ifndef TPZBANDSTRUCTMATRIX_H
#define TPZBANDSTRUCTMATRIX_H
#include "TPZStructMatrixT.h"
#include "pzstrmatrixor.h"
/**
 * @brief Implements a banded structural matrix using TPZFBMatrix as a storage format.
 * @ingroup structural
 */
template<class TVar = STATE,
         class TPar = TPZStructMatrixOR<TVar>>
class TPZBandStructMatrix : public TPZStructMatrixT<TVar>,
                            public TPar{
public:    
	using TPZStructMatrixT<TVar>::TPZStructMatrixT;
	
    TPZMatrix<TVar> * Create() override;
	TPZStructMatrix * Clone() override;

	//@{
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
    friend TPZPersistenceManager;
};

#endif //TPZBANDSTRUCTMATRIX_H
