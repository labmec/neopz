//
//  TPZPoroPermCoupling3D.hpp
//  PZ
//
//  Created by Manouchehr on 4/16/18.
//
//

#ifndef TPZPoroPermCoupling3D_h
#define TPZPoroPermCoupling3D_h

#include <stdio.h>
#include "TPZMaterial.h"
#include "TPZMatWithMem.h"
#include "TPZPoroPermMemory.h"
#include "pzdiscgal.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzvec.h"
#include "TPZSimulationData.h"
#include <iostream>
#include <cmath>
#include "pzlog.h"


class TPZPoroPermCoupling3D : public TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> {
    
protected:
    
    /** @brief define the simulation data */
    TPZSimulationData * m_SimulationData;
    
    /** @brief Problem dimension */
    int m_Dim;
    
    /** @brief solid weight */
    TPZManVector<REAL,3>  m_b;
    
    /** @brief Poison coeficient */
    REAL m_nu;
    REAL m_nuu;
    
    /** @brief first Lame Parameter */
    REAL m_lambda;
    REAL m_lambdau;

    /** @brief Bulk modulus */
    REAL m_K;
    REAL m_Ku;
    
    /** @brief Second Lame Parameter */
    REAL m_mu;
    
    /** @brief constants Biot poroelasticity */
    REAL m_alpha;
    
    /** @brief Storage coefficient poroelasticity */
    REAL m_Se;
    
    /** @brief Intact rock porosity */
    REAL m_porosity_0;
    
    /** @brief Permeability of the rock */
    REAL m_k;
    
    /** @brief Fluid viscosity */
    REAL m_eta;
    
    /** @brief coehsion of the rock */
    REAL m_c;
    
    /** @brief Friction angle */
    REAL m_phi_f;
    
    /** @brief eta parameter for Drucker-Prager model */
    REAL m_eta_dp;
    
    /** @brief xi parameter for Drucker-Prager model  */
    REAL m_xi_dp;
    
    /** @brief permeability coupling model  */
    int m_k_model;
    
    /** @brief Rock density */
    REAL m_rho_s;
    
    /** @brief Fluid density */
    REAL m_rho_f;
    
public:
    
    /** @brief Default constructor */
    TPZPoroPermCoupling3D();
    
     /** @brief Constructor based on a material id */
    TPZPoroPermCoupling3D(int matid, int dim);
    
    /** @brief Default desconstructor */
    ~TPZPoroPermCoupling3D();
    
    /** @brief Copy constructor $ */
    TPZPoroPermCoupling3D(const TPZPoroPermCoupling3D& other);
    
    /** @brief Copy assignemnt operator $ */
    TPZPoroPermCoupling3D & operator = (const TPZPoroPermCoupling3D& other);
    
    /** @brief Set the required data at each integration point */
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    /** @brief Set the required data at each integration point; BoundaryConditionData */
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    /** @brief Returns the name of the material */
    std::string Name() {
        return "TPZPoroPermCoupling3D";
    }
    
    /** returns the integrable dimension of the material */
    int Dimension() const {return m_Dim;}
    
    
    /** returns the number of state variables associated with the material */
    virtual int NStateVariables();
    
    
    /** @brief permeability model */
    REAL k_permeability(REAL &phi, REAL &k);
    
    
    /** @brief Poroelastic porosity correction */
    REAL porosoty_corrected(TPZVec<TPZMaterialData> &datavec);
    
    
    /** @brief the print function */
    void Print(std::ostream & out);
    
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    
    /** returns the number of variables associated with the variable indexed */
    int NSolutionVariables(int var);
    
    
    /** @brief Compute effective stress */
    void Compute_Sigma(REAL & l, REAL & mu, TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_u);
    
    // @brief the inner product function
    REAL Inner_Product(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & T);
    
    
    /** returns the solution associated with the var index based on
     * the finite element approximation */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right)
    {
        DebugStop();
    }
    
    
    /** @brief Parameters of rock and fluid: */
    void SetDimension(int dimension)
    {
        m_Dim = dimension;
    }
    
    
    /** @brief Parameters of rock and fluid: */
    void SetParameters(REAL perm, REAL m_porosity, REAL eta)
    {
        m_k = perm;
        m_eta = eta;
        m_porosity_0 = m_porosity;
    }
    
    
    /** @brief Set the simulation data */
    void SetSimulationData(TPZSimulationData * SimulationData)
    {
        m_SimulationData = SimulationData;
    }
    
    
    /** @brief Get the space generator */
    TPZSimulationData * SimulationData()
    {
        return m_SimulationData;
    }
    
    
    /** @brief Set the Porolastic data */
    void SetPorolasticParameters(REAL l, REAL mu, REAL l_u)
    {
        m_lambda = l;
        m_mu = mu;
        m_lambdau = l_u;
        m_K = m_lambda + (2.0/3.0)*m_mu;
        m_Ku = m_lambdau + (2.0/3.0)*m_mu;
    }
    
    
    /** @brief Set the Porolastic Engineer data */
    void SetPorolasticParametersEngineer(REAL Ey, REAL nu)
    {
        
        m_lambda = (Ey*nu)/((1.0+nu)*(1.0-2.0*nu));
        m_mu = (Ey)/(2.0*(1.0+nu));
        m_lambdau = (Ey*nu)/((1.0+nu)*(1.0-2.0*nu));
        m_K = m_lambda + (2.0/3.0)*m_mu;
        m_Ku = m_lambdau + (2.0/3.0)*m_mu;
    }
    
    
    /** @brief Set the Biot data */
    void SetBiotParameters(REAL alpha, REAL Se)
    {
        if(alpha==0){
            std::cout << "Biot constan should be at leats equal to the intact porosity, alpha = " << alpha  << std::endl;
            DebugStop();
        }
        m_alpha = alpha;
        m_Se = Se;
    }
    
   
    /** @brief Set the type Permeability model: */
    void SetKModel(int model)
    {
        m_k_model = model;
    }
    
    
    /** @brief Get the type Permeability model: */
    int KModel()
    {
        return m_k_model;
    }
    
    
    // Contribute Methods being used
    void Contribute_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void ContributeInterface(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
                             REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
                                     REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                               REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        DebugStop();
    }
   
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
    {
        DebugStop();
    }
    
    

    
public:
        
        
    
};


#endif /* TPZPoroPermCoupling_hpp */
