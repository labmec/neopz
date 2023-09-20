/**
 * @file
 * @brief Contains the TPZMatRedStructMatrix class which implements nested reduced structural matrices.
 * @author Giovane Avancini
 * @date 20/09/2023
 */

#ifndef TPZMatRedStructMatrix_H
#define TPZMatRedStructMatrix_H

#include "TPZStructMatrixT.h"
#include "pzstack.h"

#include "pzstrmatrixor.h"
/**
 * @brief Implements a sparse symmetric structural matrix using TPZMatRedMatrix as a storage format.
 * @ingroup structural
 */
template <class TVar = STATE, class TPar = TPZStructMatrixOR<TVar>>
class TPZMatRedStructMatrix : public TPZStructMatrixT<TVar>,
                              public TPar
{

public:
    TPZMatrix<TVar> *Create() override;

    TPZStructMatrix *Clone() override;

    int ClassId() const override;

    void Read(TPZStream &buf, void *context) override;

    void Write(TPZStream &buf, int withclassid) const override;
};

#endif // TPZMatRedStructMatrix_H
