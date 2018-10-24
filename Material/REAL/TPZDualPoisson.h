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
    TPZMaterial * NewMaterial();
    
    /** @brief assignment operator */
    TPZDualPoisson &operator=(const TPZDualPoisson &other);
    
    /**
     * @name Material attributes
     * @{
     */
    
    /** @brief return the euclidean dimension of the weak statement */
    int Dimension() const;
    
    /** @brief return the number of state variables associated with each trial function */
    int NStateVariables();
    
    /** @brief print all material information */
    void Print(std::ostream & out);
    
    /** @brief print material name */
    std::string Name();
    
    /** @brief fill requirements for volumetric contribute methods */
    void FillDataRequirements(TPZMaterialData &data);
    
    /** @brief fill requirements for bounadry contribute methods */
    void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data);
    
    /** @brief fill requirements for volumetric contribute methods multiphsycis mesh */
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec);
    
    /** @brief fill requirements for boundary contribute methods multiphsycis mesh */
    void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec);
    
    /** @brief unique class identifier */
    public:
    
    virtual int ClassId() const;

    
    /** @brief write class in disk */
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /** @brief write class from disk */
    void Read(TPZStream &buf, void *context);
    
    
    /** @} */
    
    /**
     * @name Contribute methods (weak formulation) no multiphysics mesh
     * @{
     */
    
    /** @brief Volumetric contribute without jacobian matrix */
    void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    /** @brief Volumetric contribute */
    void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef);
    
    /** @brief Boundary contribute without jacobian matrix */
    void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /** @brief Boundary contribute */
    void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /** @} */
    
    
    /**
     * @name Contribute methods (weak formulation) multiphysics mesh
     * @{
     */
    
    /** Computes the divergence over the parametric space */
    void ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &DivergenceofU);
    
    /** @brief Volumetric contribute without jacobian matrix */
    void Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    /** @brief Volumetric contribute */
    void Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef);
    
    /** @brief Boundary contribute without jacobian matrix */
    void ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /** @brief Boundary contribute */
    void ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /** @} */
    
    /**
     * @name Post-processing methods
     * @{
     */
    
    /** @brief Number of errors being computed = { 0-> h1, 1 ->L2 primal, 2 L2 dual} */
    int NEvalErrors();
    
    /** @brief Variable index based on variable naming */
    int VariableIndex(const std::string &name);
    
    /** @brief size of the current variable (1 -> scalar, 3-> vector, 9 ->  Tensor ) */
    int NSolutionVariables(int var);
    
    /** @brief Postprocess required variables */
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    
    /** @brief Postprocess required variables multiphysics */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    /** @brief Compute errors, no comments!!! */
    void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZFMatrix<STATE> &du, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &error);
    
    void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors);
    
    /** @} */
    
    
};

#endif /* TPZDualPoisson_hpp */
