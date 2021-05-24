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

class TPZDarcyFlowInterfaceBC;

class TPZDarcyFlowInterface : public TPZSavable {

public:
    using TInterfaceBC = TPZDarcyFlowInterfaceBC;

    void SetPermeabilityFunction(REAL constant);

    void SetPermeabilityFunction(PermeabilityFunctionType &perm_function);

    [[nodiscard]] int ClassId() const override;

    void Read(TPZStream &buf, void *context) override {};

    void Write(TPZStream &buf, int withclassid) const override {};

protected:
    PermeabilityFunctionType fPermeabilityFunction;
};

class TPZMaterial;

class TPZDarcyFlowInterfaceBC : public TPZDarcyFlowInterface {
protected:
    static void SetMaterialImpl(TPZMaterial *) {}
};

#endif //TPZISOTROPICDARCYFLOWINTERFACE_H
