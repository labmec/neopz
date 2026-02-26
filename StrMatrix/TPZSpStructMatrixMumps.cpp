/**
 * @file
 * @brief Implementation of TPZSpStructMatrixMumps.
 * Only the matrix-creation and post-setup steps differ from the base class;
 * all sparsity-pattern logic is inherited from TPZSpStructMatrix.
 */

#include "TPZSpStructMatrixMumps.h"
#include "TPZYSMPMumps.h"
#include "TPZStructMatrixOMPorTBB.h"
#include "pzstrmatrixflowtbb.h"
#include "pzstrmatrixot.h"

template<class TVar, class TPar>
TPZStructMatrix * TPZSpStructMatrixMumps<TVar,TPar>::Clone(){
    return new TPZSpStructMatrixMumps(*this);
}

template<class TVar, class TPar>
TPZFYsmpMatrix<TVar> * TPZSpStructMatrixMumps<TVar,TPar>::NewSparseMatrix(const int64_t neq) const {
    return new TPZFYsmpMatrixMumps<TVar>(neq, neq);
}

template class TPZSpStructMatrixMumps<STATE, TPZStructMatrixOR<STATE>>;
template class TPZSpStructMatrixMumps<STATE, TPZStructMatrixOT<STATE>>;
template class TPZSpStructMatrixMumps<STATE, TPZStructMatrixTBBFlow<STATE>>;
template class TPZSpStructMatrixMumps<STATE, TPZStructMatrixOMPorTBB<STATE>>;

template class TPZSpStructMatrixMumps<CSTATE, TPZStructMatrixOR<CSTATE>>;
template class TPZSpStructMatrixMumps<CSTATE, TPZStructMatrixOT<CSTATE>>;
template class TPZSpStructMatrixMumps<CSTATE, TPZStructMatrixTBBFlow<CSTATE>>;
template class TPZSpStructMatrixMumps<CSTATE, TPZStructMatrixOMPorTBB<CSTATE>>;
