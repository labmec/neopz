//
// Created by Gustavo A. Batistela on 01/06/2021.
//

#ifndef TPZNULLMATERIALCS_H
#define TPZNULLMATERIALCS_H

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"

/**
   @brief Material used in atomic meshes when working with combined approximation spaces.
   * It does not compute any contribution or shape functions.
   */
template<class TVar=STATE>
class TPZNullMaterialCS : public TPZMatBase<TVar, TPZMatCombinedSpacesT<TVar>> {
    using TBase = TPZMatBase<TVar, TPZMatCombinedSpacesT<TVar>>;
public:

    explicit TPZNullMaterialCS(int matid)
            : TPZRegisterClassId(&TPZNullMaterialCS::ClassId), TBase(matid) {
        fDim = 1;
        fNState = 1;
    }

    TPZNullMaterialCS(int matid, int dimension, int nstate)
            : TPZRegisterClassId(&TPZNullMaterialCS::ClassId),
              TBase(matid) {
#ifdef PZDEBUG
        if (dimension < 1 || dimension >3) {
            DebugStop();
        }
#endif
        fDim = dimension;
        fNState = nstate;
    }

    /** @brief Default constructor */
    TPZNullMaterialCS() = default;

    /** @brief Returns the name of the material */
    [[nodiscard]] std::string
    Name() const override { return "TPZNullMaterialCS"; }

    /** @brief Returns the integrable dimension of the material */
    [[nodiscard]] int Dimension() const override { return fDim; }

    /** @brief Sets the material dimension */
    void SetDimension(int dim);

    /** @brief Returns the number of state variables associated with the material */
    [[nodiscard]] int NStateVariables() const override {
        return fNState;
    }

    /** @brief Sets the number of state variables */
    void SetNStateVariables(int nstate) {
        fNState = nstate;
    }

    void FillDataRequirements(TPZVec<TPZMaterialDataT<TVar>> &datavec) const override;

    void Solution(const TPZVec<TPZMaterialDataT<TVar>> &datavec, int var,
                  TPZVec<TVar> &sol) override {}

    void Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec, REAL weight,
                    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override {}

    void ContributeBC(const TPZVec<TPZMaterialDataT<TVar>> &datavec, REAL weight,
                      TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc) override {}

    /** @brief Prints out the data associated with the material */
    void Print(std::ostream &out) const override;

    /** @brief Returns the variable index associated with the name */
    [[nodiscard]] int
    VariableIndex(const std::string &name) const override { return 0; }

    /**
	 * @brief Returns the number of variables associated with the variable indexed by var.
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
    [[nodiscard]] int NSolutionVariables(int var) const override { return 0; }

    /** @brief Unique identifier for serialization purposes */
    [[nodiscard]] int ClassId() const override;

    /** @brief Saves the element data to a stream */
    void Write(TPZStream &buf, int withclassid) const override;

    /** @brief Reads the element data from a stream */
    void Read(TPZStream &buf, void *context) override;

    [[nodiscard]] TPZMaterial *NewMaterial() const override {
        return new TPZNullMaterialCS(*this);
    }

protected:


    /** @brief Problem dimension */
    int fDim = -1;

    /// Number of state variables
    int fNState = 1;
};

extern template class TPZNullMaterialCS<STATE>;
extern template class TPZNullMaterialCS<CSTATE>;

#endif //TPZNULLMATERIALCS_H
