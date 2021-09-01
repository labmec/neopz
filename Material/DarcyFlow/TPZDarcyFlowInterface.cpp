//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZDarcyFlowInterface.h"

void TPZDarcyFlowInterface::SetConstantPermeability(const STATE constant) {
    fConstantPermeability = constant;
}

void TPZDarcyFlowInterface::SetPermeabilityFunction(PermeabilityFunctionType &perm_function) {
    fPermeabilityFunction = perm_function;
}

int TPZDarcyFlowInterface::ClassId() const {
    return Hash("TPZDarcyFlowInterface");
}
