//
// Created by Gustavo Batistela on 5/13/21.
//

#ifndef TPZISOTROPICDARCYFLOWINTERFACE_H
#define TPZISOTROPICDARCYFLOWINTERFACE_H

#include <functional>
#include "pzreal.h"
#include "pzvec.h"
#include "pzfmatrix.h"

using PermeabilityFunctionType = std::function<void(const TPZVec<REAL> &coord,
                                                    TPZMatrix<REAL> &K,
                                                    TPZMatrix<REAL> &InvK)>;

class TPZDarcyFlowInterface {

public:
    void SetPermeabilityFunction(REAL constant);

    void SetPermeabilityFunction(PermeabilityFunctionType &perm_function);

protected:
    PermeabilityFunctionType fPermeabilityFunction;

};

#endif //TPZISOTROPICDARCYFLOWINTERFACE_H
