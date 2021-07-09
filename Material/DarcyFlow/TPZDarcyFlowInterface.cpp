//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZDarcyFlowInterface.h"

void TPZDarcyFlowInterface::SetPermeabilityFunction(const REAL constant) {
    auto perm_func = [constant](const TPZVec<REAL> &coord, TPZMatrix<REAL> &K, TPZMatrix<REAL> &InvK) {
        for (int i = 0; i < 3; i++) {
            K(i, i) = constant;
            InvK(i, i) = 1 / constant;
        }
    };
    fPermeabilityFunction = perm_func;
}

void TPZDarcyFlowInterface::SetPermeabilityFunction(PermeabilityFunctionType &perm_function) {
    fPermeabilityFunction = perm_function;
}

int TPZDarcyFlowInterface::ClassId() const {
    return Hash("TPZDarcyFlowInterface");
}
