//
//  TPZPoroPermCoupling.hpp
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#ifndef TPZPoroPermCoupling_h
#define TPZPoroPermCoupling_h

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


class TPZPoroPermCoupling : public TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> {
    
protected:
    
    /** @brief define the simulation data */
    TPZSimulationData * m_SimulationData;
    
    /** @brief Problem dimension */
    int m_Dim;
    
    /** @brief solid weight */
    TPZManVector<REAL,2>  m_b;
    
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
    
    /** @Drucker Prager property */
    REAL m_eta_dp;
    REAL m_xi_dp;
    
    
    /** @brief permeability coupling model  */
    int m_k_model;
    
    /** @brief Uses plain stress
     * @note \f$m_PlaneStress = 1\f$ => Plain stress state
     * @note \f$m_PlaneStress != 1\f$ => Plain Strain state
     */
    int m_PlaneStress;
    
    
    /** @brief Rock density */
    REAL m_rho_s;
    
    /** @brief Fluid density */
    REAL m_rho_f;
    

    
public:
    
    TPZPoroPermCoupling();
    
    TPZPoroPermCoupling(int matid, int dim);
    
    ~TPZPoroPermCoupling();
    
    /** @brief Copy constructor $ */
    TPZPoroPermCoupling(const TPZPoroPermCoupling& other);
    
    /** @brief Copy assignemnt operator $ */
    TPZPoroPermCoupling & operator = (const TPZPoroPermCoupling& other);
    
    
    
    
    void Print(std::ostream & out);
    
    std::string Name() { return "TPZPoroPermCoupling"; }
    
    int Dimension() const {return m_Dim;}
    
    virtual int NStateVariables();
    
    /** @brief dimension of the model: */
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
    
    /** @brief Get the simulation data */
    TPZSimulationData * SimulationData()
    {
        return m_SimulationData;
    }
    
    /** @brief Set the porolastic parameters data */
    void SetPorolasticParameters(REAL l, REAL mu, REAL l_u)
    {
        m_lambda = l;
        m_mu = mu;
        m_lambdau = l_u;
        m_K = m_lambda + (2.0/3.0)*m_mu;
        m_Ku = m_lambdau + (2.0/3.0)*m_mu;
    }
    
    /** @brief Set the porolastic engineer parameters data */
    void SetPorolasticParametersEngineer(REAL Ey, REAL nu)
    {
        
        m_lambda = (Ey*nu)/((1.0+nu)*(1.0-2.0*nu));
        m_mu = (Ey)/(2.0*(1.0+nu));
        m_lambdau = (Ey*nu)/((1.0+nu)*(1.0-2.0*nu));
        m_K = m_lambda + (2.0/3.0)*m_mu;
        m_Ku = m_lambdau + (2.0/3.0)*m_mu;
    }
    
    /** @brief Set the Biot parameters data */
    void SetBiotParameters(REAL alpha, REAL Se)
    {
        if(alpha==0){
            std::cout << "Biot constan should be at leats equal to the intact porosity, alpha = " << alpha  << std::endl;
            DebugStop();
        }
        m_alpha = alpha;
        m_Se = Se;
    }
    

    
    
    /** @Drucker Prager parameters */
    void SetDruckerPragerParameters(REAL phi_f, REAL c)
    {
        m_phi_f = phi_f;
        m_c = c;
        
        // Outer edges condition
        m_eta_dp = 6.0*(sin(m_phi_f))/(sqrt(3.0)*(3.0-sin(m_phi_f)));
        m_xi_dp = 6.0*(cos(m_phi_f))/(sqrt(3.0)*(3.0-sin(m_phi_f)));
        
//        // Inner edges condition
//        m_eta_dp = 6.0*(sin(m_phi_f))/(sqrt(3.0)*(3.0+sin(m_phi_f)));
//        m_xi_dp = 6.0*(cos(m_phi_f))/(sqrt(3.0)*(3.0+sin(m_phi_f)));
//        
//        // Plane strain condition
//        m_eta_dp = 3.0*(tan(m_phi_f))/(sqrt(9.0+12.0*tan(m_phi_f)*tan(m_phi_f)));
//        m_xi_dp = 3.0/(sqrt(9.0+12.0*tan(m_phi_f)*tan(m_phi_f)));
    }

    
    
    /** @brief Set plane problem
     * planestress = 1 => Plain stress state
     * planestress != 1 => Plain Strain state
     */
    void SetPlaneProblem(int planestress)
    {
        m_PlaneStress = planestress;
    }
    
    
    /** @brief set the peremability models: */
    void SetKModel(int model)
    {
        m_k_model = model;
    }
    
    /** @brief return the peremability models: */
    int KModel()
    {
        return m_k_model;
    }
    
    /** @brief permeability correction model */
    REAL k_permeability(REAL &phi, REAL &k);
    
    /** @brief porosity correction model */
    REAL porosity_corrected(TPZVec<TPZMaterialData> &datavec);
    
    
    
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);

    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    void ContributeBC_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    void ContributeBC_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    void Compute_Sigma(TPZFMatrix<REAL> & S_eff,TPZFMatrix<REAL> & Grad_u,  REAL p_ex);
    void Compute_Sigma(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_v);
    
    void Compute_Sigma(REAL & l, REAL & mu, TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_u);

    
    REAL Inner_Product(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & T);
    
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right)
    {
        DebugStop();
    }

    void ContributeInterface(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
                             REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
                                     REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
        DebugStop();
    }
    
    void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef) {
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                               REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
   
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    

    
public:
    
    /** @brief mean stress */
    REAL p(TPZFMatrix<REAL> T);
    
    /** @brief mean stress */
    TPZFMatrix<REAL> s(TPZFMatrix<REAL> T);
    
    /** @brief J2 invariant stress */
    REAL J2(TPZFMatrix<REAL> T);
    
    /** @brief J3 invariant stress */
    REAL J3(TPZFMatrix<REAL> T);
    
    /** @brief theta */
    REAL theta(TPZFMatrix<REAL> T);
    
    /** @brief Phi Mohr-Coulomb */
    REAL Phi_MC(TPZFMatrix<REAL> T);
    
    /** @brief Phi Drucker-Prager */
    REAL Phi_DP(TPZFMatrix<REAL> T);
    
    /** @brief plasticity multiplier delta_gamma */
    REAL Phi_tilde_DP(TPZFMatrix<REAL> T, REAL d_gamma_guest);
    
    /** @brief plasticity multiplier delta_gamma */
    REAL Phi_tilde_DP_delta_gamma(TPZFMatrix<REAL> T, REAL d_gamma_guest);
    
    /** @brief plasticity multiplier delta_gamma using newton iterations */
    REAL delta_gamma_finder(TPZFMatrix<REAL> T, REAL d_gamma_guest);
    
    /** @brief Drucker prager strain update */
    TPZFMatrix<REAL> strain_DP(TPZFMatrix<REAL> T);
    
    /** @brief Drucker prager stress update */
    TPZFMatrix<REAL> stress_DP(TPZFMatrix<REAL> T);
    
    /** @brief Principal Stress */
    void Principal_Stress(TPZFMatrix<REAL> T, TPZFMatrix<REAL> & S);
       
    /** @brief Drucker prager elastoplastic corrector  */
    void corrector_DP(TPZFMatrix<REAL> Grad_u_n, TPZFMatrix<REAL> Grad_u, TPZFMatrix<REAL> &e_e, TPZFMatrix<REAL> &e_p, TPZFMatrix<REAL> &S);
    
    
        
    
};


#endif /* TPZPoroPermCoupling_hpp */
