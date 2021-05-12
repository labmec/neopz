//
//  TPZDualPoisson.hpp
//  PZ
//
//  Created by omar on 04/07/2016.
//
//

#ifndef TPZDualPoisson_hpp
#define TPZDualPoisson_hpp


#include <stdio.h>
#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"
    
    
/**
 * @ingroup Material that implements Continuous Galerkin Approximation of laplace equation
 * @brief This material consider exactly just laplace equation (i.e. coefficient equal to 1)
 */

/**
 * \f$ q = - grad(u)  \f$
 * \f$ div( q ) =  f  \f$
 */
class TPZDualPoisson : public TPZMaterial {
    
    
public:
    
    /** @brief default constructor  */
    TPZDualPoisson();
    
    /** @brief default desconstructor  */
    ~TPZDualPoisson();
    
    /** @brief constructor based on material id */
    TPZDualPoisson(int mat_id);
    
    /** @brief constructor based on object copy */
    TPZDualPoisson(const TPZDualPoisson &copy);
    
    /** @brief return a material object from a this object */
    TPZMaterial * NewMaterial() override;
    
    /** @brief assignment operator */
    TPZDualPoisson &operator=(const TPZDualPoisson &other);
    
    /**
     * @name Material attributes
     * @{
     */
    
    /** @brief return the euclidean dimension of the weak statement */
    int Dimension() const override;
    
    /** @brief return the number of state variables associated with each trial function */
    virtual int NStateVariables() const override;
    
    /** @brief print all material information */
    void Print(std::ostream & out) override;
    
    /** @brief print material name */
    std::string Name() override;
    
    /** @brief fill requirements for volumetric contribute methods */
    void FillDataRequirements(TPZMaterialData &data) override;
    
    /** @brief fill requirements for bounadry contribute methods */
    void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data) override;
    
    /** @brief fill requirements for volumetric contribute methods multiphsycis mesh */
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec) override;
    
    /** @brief fill requirements for boundary contribute methods multiphsycis mesh */
    void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec) override;
    
    /** @brief unique class identifier */
    public:
    
    virtual int ClassId() const override;

    
    /** @brief write class in disk */
    virtual void Write(TPZStream &buf, int withclassid) const override;
    
    /** @brief write class from disk */
    void Read(TPZStream &buf, void *context) override;
    
    
    /** @} */
    
    /**
     * @name Contribute methods (weak formulation) no multiphysics mesh
     * @{
     */
    
    /** @brief Volumetric contribute without jacobian matrix */
    void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
    
    /** @brief Volumetric contribute */
    void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef) override;
    
    /** @brief Boundary contribute without jacobian matrix */
    void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
    /** @brief Boundary contribute */
    void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
    /** @} */
    
    
    /**
     * @name Contribute methods (weak formulation) multiphysics mesh
     * @{
     */
    
    /** @brief Volumetric contribute without jacobian matrix */
    void Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
    
    /** @brief Volumetric contribute */
    void Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef) override;
    
    /** @brief Boundary contribute without jacobian matrix */
    void ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
    /** @brief Boundary contribute */
    void ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
    /** @} */
    
    /**
     * @name Post-processing methods
     * @{
     */
    
    /** @brief Number of errors being computed = { 0-> h1, 1 ->L2 primal, 2 L2 dual} */
    int NEvalErrors() override;
    
    /** @brief Variable index based on variable naming */
    int VariableIndex(const std::string &name) override;
    
    /** @brief size of the current variable (1 -> scalar, 3-> vector, 9 ->  Tensor ) */
    int NSolutionVariables(int var) override;
    
    /** @brief Postprocess required variables */
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
    
    /** @brief Postprocess required variables multiphysics */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
    
    
    void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors);
    
    /** @} */
    
    
};

#endif /* TPZDualPoisson_hpp */
