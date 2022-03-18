//
// Created by Gustavo Batistela on 5/4/21.
//

#ifndef TPZNULLMATERIAL_H
#define TPZNULLMATERIAL_H


#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"

/**
   @brief Material used in atomic meshes when working with combined approximation spaces.
   * It does not compute any contribution or shape functions.
   */
template<class TVar=STATE>
class TPZNullMaterial : public TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>> {
    using TBase = TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>>;
public:

    explicit TPZNullMaterial(int matid) : TPZRegisterClassId(&TPZNullMaterial::ClassId),
                                          TBase(matid) {
        fDim = 1;
        fNState = 1;
    }

    TPZNullMaterial(int matid, int dimension, int nstate = 1) : TPZRegisterClassId(&TPZNullMaterial::ClassId),
                                                            TBase(matid) {
#ifdef PZDEBUG
        if (dimension < 0 || dimension > 3) {
            DebugStop();
        }
#endif
        fDim = dimension;
        fNState = nstate;
    }

    /** @brief Default constructor */
    TPZNullMaterial() = default;

    /** @brief Returns the name of the material */
    [[nodiscard]] std::string Name() const override { return "TPZNullMaterial"; }

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

    void FillDataRequirements(TPZMaterialData &) const override;
    
    void GetSolDimensions(uint64_t &u_len,
                          uint64_t &du_row,
                          uint64_t &du_col) const override;

    void Solution(const TPZMaterialDataT<TVar> &data, int var, TPZVec<TVar> &sol) override;

    void Contribute(const TPZMaterialDataT<TVar> &data, REAL weight,
                            TPZFMatrix<TVar> &ek,
                            TPZFMatrix<TVar> &ef) override
    {}
    void ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                       TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                       TPZBndCondT<TVar> &bc) override
    {}
    
    /** @brief Prints out the data associated with the material */
    void Print(std::ostream &out = std::cout) const override;

    /** @brief Returns the variable index associated with the name */
    int VariableIndex(const std::string &name) const override { return 0;}

    /**
	 * @brief Returns the number of variables associated with the variable indexed by var.
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
    [[nodiscard]] int NSolutionVariables(int var) const override;

    /** @brief Unique identifier for serialization purposes */
    [[nodiscard]] int ClassId() const override;

    /** @brief Saves the element data to a stream */
    void Write(TPZStream &buf, int withclassid) const override;

    /** @brief Reads the element data from a stream */
    void Read(TPZStream &buf, void *context) override;

    /**
     * @brief Creates another material of the same type
     */
    [[nodiscard]] TPZMaterial *NewMaterial() const override;

protected:


    /** @brief Problem dimension */
    int fDim;

    /// Number of state variables
    int fNState = 1;
};

extern template class TPZNullMaterial<STATE>;
extern template class TPZNullMaterial<CSTATE>;

#endif //TPZNULLMATERIAL_H
