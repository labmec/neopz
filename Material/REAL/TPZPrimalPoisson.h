//
//  TPZPrimalPoisson.h
//  PZ
//
//  Created by omar on 04/07/2016.
//
//

#ifndef TPZPrimalPoisson_h
#define TPZPrimalPoisson_h

#include <stdio.h>
#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"


/**
 * @ingroup Material that implements Continuous Galerkin Approximation of laplace equation
 * @brief This material consider exactly just laplace equation (i.e. coefficient equal to 1)
 */

/**
 * \f$ div( - grad (u) ) =  f  \f$
 */
class TPZPrimalPoisson : public TPZMaterial {
    
    
public:
    
    /** @brief default constructor  */
    TPZPrimalPoisson();
    
    /** @brief default desconstructor  */
    ~TPZPrimalPoisson();
    
    /** @brief constructor based on material id */
    TPZPrimalPoisson(int mat_id);
    
    /** @brief constructor based on object copy */
    TPZPrimalPoisson(const TPZPrimalPoisson &copy);
    
    /** @brief return a material object from a this object */
    TPZMaterial * NewMaterial() override;

    /** @brief assignment operator */
    TPZPrimalPoisson &operator=(const TPZPrimalPoisson &other);
    
    /**
     * @name Material attributes
     * @{
     */

    /** @brief return the euclidean dimension of the weak statement */
    int Dimension() const override ;
    
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

    
    /** @brief Volumetric contribute with jacobian matrix */
    void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;

    /** @brief Volumetric contribute */
    void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef) override;

    /** @brief Boundary contribute with jacobian matrix */
    void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;

    /** @brief Boundary contribute */
    void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
    /** @brief Volumetric contribute with jacobian matrix */
    void Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
    
    /** @brief Volumetric contribute */
    void Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef) override;
    
    /** @brief Boundary contribute with jacobian matrix */
    void ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
    /** @brief Boundary contribute */
    void ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;

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

    /** @brief Compute errors, no comments!!! */
protected:
    void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZFMatrix<STATE> &du, TPZFMatrix<REAL> &axes, TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &error) override;
    
};


#endif /* TPZPrimalPoisson_hpp */
