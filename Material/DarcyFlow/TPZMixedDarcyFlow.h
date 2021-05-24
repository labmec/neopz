//
// Created by Gustavo Batistela on 5/13/21.
//

#ifndef TPZMIXEDDARCYFLOW_H
#define TPZMIXEDDARCYFLOW_H

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "TPZDarcyFlowInterface.h"

// TODO add doc

class TPZMixedDarcyFlow : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
        TPZMatErrorCombinedSpaces<STATE>, TPZDarcyFlowInterface> {

    // type alias to improve constructor readability
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
            TPZMatErrorCombinedSpaces<STATE>, TPZDarcyFlowInterface>;

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
    int NEvalErrors() override { return 5; }

    /**
     * @brief Sets problem dimension
     */
    virtual void SetDimension(int dim);

    // TODO add doc
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    // TODO add doc
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

    // TODO add doc
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut) override;

    // TODO add doc
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;

    // TODO add doc
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    // TODO add doc
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
