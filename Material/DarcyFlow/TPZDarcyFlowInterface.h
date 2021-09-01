//
// Created by Gustavo Batistela on 5/13/21.
//

#ifndef TPZDARCYFLOWINTERFACE_H
#define TPZDARCYFLOWINTERFACE_H

#include <functional>
#include "pzreal.h"
#include "pzvec.h"
#include "pzfmatrix.h"

// Alias to improve readability of the permeability function type
using PermeabilityFunctionType = std::function<STATE(const TPZVec<REAL> &coord)>;

// Forward declaration of dummy BC interface class
class TPZDarcyFlowInterfaceBC;

/**
 * @brief  This class implements the interface with the methods required to
 * handle the permeability tensor of an isotropic material.
 */
class TPZDarcyFlowInterface : public TPZSavable {

public:
    using TInterfaceBC = TPZDarcyFlowInterfaceBC;

    /**
     * @brief Set a constant permeability to the material
     * @param [in] constant permeability value
     */
    void SetConstantPermeability(STATE constant);

    /**
     * @brief Set a varying permeability field to the material
     * @param [in] perm_function function that describes the permeability field
     */
    void SetPermeabilityFunction(PermeabilityFunctionType &perm_function);

    [[nodiscard]] int ClassId() const override;

    void Read(TPZStream &buf, void *context) override {};

    void Write(TPZStream &buf, int withclassid) const override {};

protected:

    STATE fConstantPermeability{};
    PermeabilityFunctionType fPermeabilityFunction{};
};

// Dummy BC interface class
class TPZMaterial;
class TPZDarcyFlowInterfaceBC : public TPZDarcyFlowInterface {
protected:
    static void SetMaterialImpl(TPZMaterial *) {}
};

#endif //TPZISOTROPICDARCYFLOWINTERFACE_H
