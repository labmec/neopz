/**
 * @file
 * @brief Implementation of TPZSSpStructMatrixMumps.
 * Only the matrix-creation and post-setup steps differ from the base class;
 * all sparsity-pattern logic is inherited from TPZSSpStructMatrix.
 */

#include "TPZSSpStructMatrixMumps.h"
#include "TPZSYSMPMumps.h"
#include "TPZStructMatrixOMPorTBB.h"
#include "pzstrmatrixflowtbb.h"
#include "pzstrmatrixot.h"

template<class TVar, class TPar>
TPZStructMatrix * TPZSSpStructMatrixMumps<TVar,TPar>::Clone(){
    return new TPZSSpStructMatrixMumps(*this);
}

template<class TVar, class TPar>
TPZSYsmpMatrix<TVar> * TPZSSpStructMatrixMumps<TVar,TPar>::NewSparseMatrix(const int64_t neq) const {
    return new TPZSYsmpMatrixMumps<TVar>(neq, neq);
}

template class TPZSSpStructMatrixMumps<STATE, TPZStructMatrixOR<STATE>>;
template class TPZSSpStructMatrixMumps<STATE, TPZStructMatrixOT<STATE>>;
template class TPZSSpStructMatrixMumps<STATE, TPZStructMatrixTBBFlow<STATE>>;
template class TPZSSpStructMatrixMumps<STATE, TPZStructMatrixOMPorTBB<STATE>>;