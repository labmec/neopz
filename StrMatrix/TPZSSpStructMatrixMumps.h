/**
 * @file
 * @brief Contains the TPZSSpStructMatrixMumps class which implements sparse structural matrices.
 */

#ifndef TPZSSpStructMatrixMumps_H
#define TPZSSpStructMatrixMumps_H

#include "TPZSSpStructMatrix.h"
#include "pzstrmatrixor.h"

/**
 * @brief Sparse symmetric structural matrix using MUMPS as the direct solver.
 * Inherits all sparsity pattern setup from TPZSSpStructMatrix; only the
 * matrix object creation and post-setup steps are overridden.
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSSpStructMatrixMumps : public TPZSSpStructMatrix<TVar,TPar> {
public:
    using TPZSSpStructMatrix<TVar,TPar>::TPZSSpStructMatrix;

    TPZStructMatrix * Clone() override;

protected:
    TPZSYsmpMatrix<TVar> * NewSparseMatrix(int64_t neq) const override;
};

#endif //TPZSSpStructMatrixMumps_H