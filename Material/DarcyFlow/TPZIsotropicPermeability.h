//
// Created by Gustavo Batistela on 5/13/21.
//

#ifndef TPZISOTROPICPERMEABILITY_H
#define TPZISOTROPICPERMEABILITY_H

#include <functional>
#include "pzreal.h"
#include "pzvec.h"
#include "pzfmatrix.h"

// Alias to improve readability of the permeability function type
using PermeabilityFunctionType = std::function<STATE(const TPZVec<REAL> &coord)>;

// Forward declaration of dummy BC interface class
class TPZIsotropicPermeabilityBC;

/**
 * @brief  This class implements the interface with the methods required to
 * handle the permeability field of an isotropic material.
 */
class TPZIsotropicPermeability : public TPZSavable {

public:
    using TInterfaceBC = TPZIsotropicPermeabilityBC;

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

    /**
     * @brief Return the permeability value at a coordinate
     * @param [in] coord coordinate of interest
     */
    STATE GetPermeability(const TPZVec<REAL> &coord);

    [[nodiscard]] int ClassId() const override;

    void Read(TPZStream &buf, void *context) override {};

    void Write(TPZStream &buf, int withclassid) const override {};

private:

    // Member variable to describe a constant permeability field
    STATE fConstantPermeability = 1.;

    // Member variable to describe a varying permeability field
    PermeabilityFunctionType fPermeabilityFunction{};
};

// Dummy BC interface class
class TPZMaterial;
class TPZIsotropicPermeabilityBC : public TPZIsotropicPermeability {
protected:
    // this method is your chance to verify if the material to which this
    // BC interface applies is compatible with this boundary interface
    // it is called in the method SetMaterial of class TPZBndCondBase
    static void SetMaterialImpl(TPZMaterial *) {}
};

#endif //TPZISOTROPICPERMEABILITY_H
