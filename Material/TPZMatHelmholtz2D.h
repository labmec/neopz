/**
 * @file TPZMatHelmholtz2D.h
 * @brief Header file for class TPZMatHelmholtz2D.\n
 */

#ifndef TPZMATHELMHOLTZ2D_H
#define TPZMATHELMHOLTZ2D_H

#include "TPZVecL2.h"
/**
 * @ingroup material
 * @brief This class implements the inhomogeneous Helmholtz wave equation in 2D.
 */
class TPZMatHelmholtz2D : public TPZVecL2 {
  protected:
    const STATE fC;
    const REAL fScale;
  public:
    TPZMatHelmholtz2D(int id, const STATE &, const REAL &scale = 1.);

    explicit TPZMatHelmholtz2D(int id);

    /** @brief Default constructor */
    TPZMatHelmholtz2D();

    /** @brief Creates a material object based on the referred object and
     * inserts it in the vector of material pointers of the mesh. */
    /** Upon return vectorindex contains the index of the material object within
     * the vector */
    TPZMatHelmholtz2D(const TPZMatHelmholtz2D &mat);
    /** @brief Default destructor */
    virtual ~TPZMatHelmholtz2D();

    /** @brief Returns the name of the material */
    std::string Name() override { return "TPZMatHelmholtz2D"; }

    /** @brief Returns the integrable dimension of the material */
    int Dimension() const override { return 2; }

    /** @brief Returns the number of state variables associated with the
     * material
     */
    int NStateVariables() const override { return 1; }

  public:
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector
     * at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    void Contribute(TPZMaterialData &data, REAL weight,
                            TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector
     * at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    void ContributeBC(TPZMaterialData &data, REAL weight,
                              TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                              TPZBndCond &bc) override;

    /**
     * @brief This method defines which parameters need to be initialized in
     * order to compute the contribution of an element
     * @param datavec [out] vector of TPZMaterialData, each position will
     * specifie the requirements for its correspondent state variable
     */
    void FillDataRequirements(TPZMaterialData &data) override {
        data.SetAllRequirements(false);
    }

    /** @brief This method defines which parameters need to be initialized in
     * order to compute the contribution of the boundary condition */
    void FillBoundaryConditionDataRequirement(int type,
                                                      TPZMaterialData &data) override {
        data.SetAllRequirements(false);
        data.fNeedsNormal = true;
    }

    /** @brief Gets the order of the integration rule necessary to integrate an
     * element with polinomial order p */
    int IntegrationRuleOrder(int elPMaxOrder) const override;

    int VariableIndex(const std::string &name) override;

    /**
     * @brief Returns the number of variables associated with the variable
     * indexed by var.
     * @param var Index variable into the solution, is obtained by calling
     * VariableIndex
     */
    int NSolutionVariables(int var) override;

    /** @brief Returns the solution associated with the var index based on the
     * finite element approximation */
    void Solution(TPZMaterialData &data, int var,
                          TPZVec<STATE> &Solout) override;

    void ErrorsHdiv(TPZMaterialData &data, TPZVec<STATE> &u_exact,
                    TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &values) override;
    int NEvalErrors() override { return 3; } // l2 and hcurl
};

#endif
