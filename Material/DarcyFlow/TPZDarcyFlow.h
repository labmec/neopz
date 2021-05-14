//
// Created by Gustavo Batistela on 5/13/21.
//

#ifndef TPZDARCYFLOW_H
#define TPZDARCYFLOW_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatErrorSingleSpace.h"
#include "TPZDarcyFlowInterface.h"

// TODO: document class and differential equation

class TPZDarcyFlow : public TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>, TPZMatErrorSingleSpace<STATE>>,
                     public TPZDarcyFlowInterface {

    // type alias to improve constructor readability
    using TBase = TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>, TPZMatErrorSingleSpace<STATE>>;

public:
    /**
     * @brief Default constructor
     */
    TPZDarcyFlow();

    /**
	 * @brief Class constructor
	 * @param [in] id material id
	 * @param [in] dim problem dimension
	 */
    TPZDarcyFlow(int id, int dim);

    /**
	 * @brief Returns a 'std::string' with the name of the material
	 */
    [[nodiscard]] std::string Name() const override { return "TPZDarcyFlow"; }

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
    int NEvalErrors() override { return 6; }

    /**
     * @brief Sets problem dimension
     */
    virtual void SetDimension(int dim);

    // TODO document method
    void Contribute(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    // TODO document method
    void ContributeBC(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;

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

    // TODO document method
    void Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<STATE> &solOut) override;

    // TODO document method
    void GetSolDimensions(uint64_t &u_len, uint64_t &du_row, uint64_t &du_col) const override;

    // TODO document method
    void Errors(const TPZVec<REAL> &x, const TPZVec<STATE> &sol, const TPZFMatrix<STATE> &dsol,
                const TPZFMatrix<REAL> &axes, TPZVec<REAL> &errors) override;

    // TODO add doc
    void FillDataRequirements(TPZMaterialData &data) const override;

    // TODO add doc
    void FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data) const override;

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

#endif //TPZDARCYFLOW_H
