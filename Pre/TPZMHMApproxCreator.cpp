#include "TPZMHMApproxCreator.h"


TPZMHMApproxCreator::TPZMHMApproxCreator() {
}

TPZMHMApproxCreator::TPZMHMApproxCreator(TPZVec<int64_t>& elementPartition) {
    fElementPartition = elementPartition;
}

TPZMHMApproxCreator::~TPZMHMApproxCreator() {
}
