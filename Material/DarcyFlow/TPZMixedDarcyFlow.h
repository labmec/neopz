//
// Created by Gustavo Batistela on 5/13/21.
//

#ifndef TPZMIXEDDARCYFLOW_H
#define TPZMIXEDDARCYFLOW_H

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "TPZIsotropicPermeability.h"

/**
 * @ingroup material
 * @brief This class implements a mixed approximation for the Darcy flow equation for isotropic materials.
 *
 * The Darcy flow equation is given by: \f[\nabla \cdot \boldsymbol{\sigma} = f,\f]
 * where \f$\boldsymbol{\sigma} = -K \nabla u\f$ and \f$u\f$ are the flux and pressure field to be solved respectively,
 * \f$K\f$ is the permeability tensor and \f$f\f$ is the source term.
 *
 * @see TPZDarcyFlow For an H1-conforming approximation.
 * @see TPZIsotropicPermeability For setting the permeability field.
 */

class TPZMixedDarcyFlow : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
        TPZMatErrorCombinedSpaces<STATE>, TPZIsotropicPermeability> {

    // type alias to improve constructor readability
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
            TPZMatErrorCombinedSpaces<STATE>, TPZIsotropicPermeability>;

public:
    /**
     * @brief Default constructor
     */
    TPZMixedDarcyFlow();

    /**
	 * @brief Class constructor
	 * @param [in] id material id
	 * @param [in] dim problem dimension
	 */
    [[maybe_unused]] TPZMixedDarcyFlow(int id, int dim);

    /**
             copy constructor
     */
    TPZMixedDarcyFlow(const TPZMixedDarcyFlow &copy);
    /**
             copy constructor
     */
    TPZMixedDarcyFlow &operator=(const TPZMixedDarcyFlow &copy);
    /**
	 * @brief Returns a 'std::string' with the name of the material
	 */
    [[nodiscard]] std::string Name() const override { return "TPZMixedDarcyFlow"; }

    /**
	 * @brief Returns the problem dimension
	 */
    [[nodiscard]] int Dimension() const override { return this->fDim; }

    /**
	 * @brief Returns the number of state variables
	 */
    [[nodiscard]] int NStateVariables() const override { return 1; }

    /**
	 * @brief Returns the number of errors to be evaluated
     *
     * Returns the number of errors to be evaluated, that is, the number of error norms associated
     * with the problem.
     */
    int NEvalErrors() const override { return 5; }

    /**
     * @brief Sets problem dimension
     */
    virtual void SetDimension(int dim);

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    /**
     * @brief Returns an integer associated with a post-processing variable name
     * @param [in] name string containing the name of the post-processing variable. Ex: "Pressure".
     */
    [[nodiscard]] int VariableIndex(const std::string &name) const override;

    /**
     * @brief Returns an integer with the dimension of a post-processing variable
     * @param [in] var index of the post-processing variable, according to TPZDarcyFlow::VariableIndex method.
     */
    [[nodiscard]] int NSolutionVariables(int var) const override;

    /**
     * @brief Returns the solution associated with the var index based on the
     * finite element approximation at a point
     * @param [in] datavec material data associated with a given integration point
     * @param [in] var index of the variable to be calculated
     * @param [out] solOut vector to store the solution
     */
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut) override;

    /**
     * @brief Calculates the approximation error at a point
     * @param [in] data material data of the integration point
     * @param [out] errors calculated errors
     */
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;

    /*
     * @brief Fill requirements for volumetric contribute
     */
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    /*
     * @brief Fill requirements for boundary contribute
     */
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    /**
     * @brief Returns an unique class identifier
     */
    [[nodiscard]] int ClassId() const override;

    /**
     * @brief Creates another material of the same type
     */
    [[nodiscard]] TPZMaterial *NewMaterial() const override;

    /**
     * @brief Prints data associated with the material.
     */
    void Print(std::ostream & out) const override;

protected:
    /**
     * @brief Problem dimension
     */
    int fDim;

};

#endif //TPZMIXEDDARCYFLOW_H
