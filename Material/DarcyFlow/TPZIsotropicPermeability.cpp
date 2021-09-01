//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZIsotropicPermeability.h"

void TPZIsotropicPermeability::SetConstantPermeability(const STATE constant) {
    fConstantPermeability = constant;
}

void TPZIsotropicPermeability::SetPermeabilityFunction(PermeabilityFunctionType &perm_function) {
    fPermeabilityFunction = perm_function;
}

STATE TPZIsotropicPermeability::GetPermeability(const TPZVec<REAL> &coord) {
    return fPermeabilityFunction ? fPermeabilityFunction(coord) : fConstantPermeability;
}

int TPZIsotropicPermeability::ClassId() const {
    return Hash("TPZIsotropicPermeability");
}
