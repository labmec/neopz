/**
 * @file
 * @brief Contains the TPZSpStructMatrixMumps class which implements sparse structural matrices.
 */

#ifndef TPZSpStructMatrixMumps_H
#define TPZSpStructMatrixMumps_H

#include "TPZSpStructMatrix.h"
#include "pzstrmatrixor.h"

/**
 * @brief Sparse structural matrix using MUMPS as the direct solver.
 * Inherits all sparsity pattern setup from TPZSpStructMatrix; only the
 * matrix object creation and post-setup steps are overridden.
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSpStructMatrixMumps : public TPZSpStructMatrix<TVar,TPar> {
public:
    using TPZSpStructMatrix<TVar,TPar>::TPZSpStructMatrix;

    TPZStructMatrix * Clone() override;

protected:
    TPZFYsmpMatrix<TVar> * NewSparseMatrix(int64_t neq) const override;
};

#endif //TPZSpStructMatrixMumps_H
