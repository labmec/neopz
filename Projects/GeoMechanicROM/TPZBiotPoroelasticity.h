//
//  TPZBiotPoroelasticity.h
//  PZ
//
//  Created by Omar on 2/25/17.
//
//

#ifndef TPZBiotPoroelasticity_h
#define TPZBiotPoroelasticity_h

#include <stdio.h>
#include <iostream>
#include <cmath>
#include "pzmatwithmem.h"
#include "pzdiscgal.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzvec.h"
#include "TPZSimulationData.h"
#include "TPZPoroPermMemory.h"
#include "pzaxestools.h"

class TPZBiotPoroelasticity : public TPZMatWithMem< TPZPoroPermMemory,TPZDiscontinuousGalerkin> {
    
private:
    
    /** @brief define the simulation data */
    TPZSimulationData * fSimulationData;
    
    /** @brief material dimension */
    int fdimension;
    
    /** @brief Material parameter for symetric system mixed */
    REAL fMFsymetric;

    /** @brief Material parameter for symetric coupled system */
    REAL fCsymetric;
    
    /** @brief first Lame Parameter */
    REAL flambda;
    
    /** @brief undrained first Lame Parameter */
    REAL flambdau;
    
    /** @brief Bulk modulus */
    REAL fK;

    /** @brief undrained Bulk modulus */    
    REAL fKu;
    
    /** @brief Second Lame Parameter */
    REAL fmu;
    
    /** @brief constants Biot poroelasticity */
    REAL falpha;
    
    /** @brief Storage coefficient poroelasticity */
    REAL fSe;
    
    /** @brief Intact rock porosity */
    REAL fporosity_0;
    
    /** @brief Permeability of the rock */
    REAL fk;
    
    /** @brief Fluid viscosity */
    REAL feta;

    /** @brief Define the use of mixed formulation */
    bool fIsMixedQ;
    
    /** @brief Define the use of symmetric factors */
    bool fIsSymmetricQ;
    
public:
    
    
    /** @brief Default constructor */
    TPZBiotPoroelasticity();
    
    /** @brief Constructor based on a material id */
    TPZBiotPoroelasticity(int matid, int dimension);
    
    /** @brief Constructor based on a Biot Poroelasticity  object */
    TPZBiotPoroelasticity(const TPZBiotPoroelasticity &mat);
    
    /** @brief Constructor based on a Biot Poroelasticity  object */
    TPZBiotPoroelasticity &operator=(const TPZBiotPoroelasticity &mat)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief Default destructor */
    ~TPZBiotPoroelasticity();
    
    /** @brief Set the required data at each integration point */
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec);
    
    /** @brief Set the required data at each integration point */
    void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec);
    
    /** @brief Returns the name of the material */
    std::string Name() {
        return "TPZBiotPoroelasticity";
    }
    
    /** returns the integrable dimension of the material */
    int Dimension() const {return fdimension;}
    
    /** returns the number of state variables associated with the material */
    int NStateVariables() {return 1;}
    
    virtual TPZMaterial *NewMaterial()
    {
        return new TPZBiotPoroelasticity(*this);
    }
    
    /** print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    /** returns the number of variables associated with the variable
     indexed by var.  var is obtained by calling VariableIndex */
    int NSolutionVariables(int var);
    
    /** returns the solution associated with the var index based on
     * the finite element approximation */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
    // Attribute
    
    void SetPorolasticParameters(REAL l, REAL mu, REAL l_u)
    {
        flambda = l;
        fmu = mu;
        flambdau = l_u;
        fK = flambda + (2.0/3.0)*fmu;
        fKu = flambdau + (2.0/3.0)*fmu;
    }
    
    void SetBiotParameters(REAL alpha, REAL Se)
    {
        if(alpha==0){
            std::cout << "Biot constan should be at leats equal to the intact porosity, alpha = " << alpha  << std::endl;
            DebugStop();
        }
        
        falpha = alpha;
        fSe = Se;
    }
    
    /** @brief Parameters of rock and fluid: */
    void SetFlowParameters(REAL perm, REAL fporosity, REAL eta)
    {
        fk = perm;
        feta = eta;
        fporosity_0 = fporosity;
    }
    
    /** @brief Set the simulation data */
    void SetSimulationData(TPZSimulationData * SimulationData)
    {
        fSimulationData = SimulationData;
    }
    
    void SetMixedFormulation(bool IsMixedQ){
        fIsMixedQ = IsMixedQ;
    }
    
    /** @brief Define the use of symmetric factors */
    void SetSymmetricFormulation(bool IsSymmetricQ){
        fIsSymmetricQ = IsSymmetricQ;
    }
    
    /** @brief Get the space generator */
    TPZSimulationData * SimulationData()
    {
        return fSimulationData;
    }
    
    void Compute_Sigma(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_v);
    
    void Compute_Sigma_fast(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_v);
    
    // Contribute Methods being used
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    // Reduce basis methods
    void ContributeRB(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void ContributeRB_BC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    
    /** Computes the divergence over the parametric space */
    void ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &Divergence_of_q);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    void ContributeMF(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    void ContributeMF(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeMFBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /** returns the solution associated with the var index based on
     * the finite element approximation */
    void SolutionMF(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
    /** @brief Unique identifier for serialization purposes */
    int ClassId() const;
    
    /** @brief Save object data to a stream */
    void Write(TPZStream &buf, int withclassid);
    
    /** @brief Read object data from a stream */
    void Read(TPZStream &buf, void *context);
    
};

#endif /* TPZBiotPoroelasticity_h */
