//
//  TPZLinearElliptic.hpp
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#ifndef TPZLinearElliptic_h
#define TPZLinearElliptic_h

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
#include <algorithm>    // std::max
#include "pzlog.h"


class TPZLinearElliptic : public TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> {
    
protected:
    
    /** @brief define the simulation data */
    TPZSimulationData * fSimulationData;
    
    /** @brief Problem dimension */
    int fDim;
    
    /** @brief solid weight */
    TPZManVector<REAL,2>  fb;
    
    /** @brief Poison coeficient */
    REAL fnu;
    REAL fnuu;
    
    /** @brief first Lame Parameter */
    REAL flambda;
    REAL flambdau;
    
    /** @brief Bulk modulus */
    REAL fK;
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
    
    /** @brief coehsion of the rock */
    REAL fc;
    
    /** @brief Friction angle */
    REAL fphi_f;
    
    /** @brief eta parameter for Drucker-Prager model */
    REAL feta_dp;
    
    /** @brief xi parameter for Drucker-Prager model  */
    REAL fxi_dp;
    
    /** @brief permeability coupling model  */
    int fk_model;
    
    /** @brief Uses plain stress
     * @note \f$fPlaneStress = 1\f$ => Plain stress state
     * @note \f$fPlaneStress != 1\f$ => Plain Strain state
     */
    int fPlaneStress;
    
    /** @brief Rock density */
    REAL frho_s;
    
    /** @brief Fluid density */
    REAL frho_f;
    
public:
    
    TPZLinearElliptic();
    
    TPZLinearElliptic(int matid, int dim);
    
    ~TPZLinearElliptic();
    
    void Print(std::ostream & out);
    
    std::string Name() { return "TPZLinearElliptic"; }
    
    int Dimension() const {return fDim;}
    
    virtual int NStateVariables();
    
    /** @brief Parameters of rock and fluid: */
    void SetDimension(int dimension)
    {
        fDim = dimension;
    }
    
    
    /** @brief Parameters of rock and fluid: */
    void SetParameters(REAL perm, REAL fporosity, REAL eta)
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
    
    /** @brief Get the space generator */
    TPZSimulationData * SimulationData()
    {
        return fSimulationData;
    }
    
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
    
    void SetDruckerPragerParameters(REAL phi_f, REAL c)
    {
        fphi_f = phi_f;
        fc = c;
        feta_dp = 6.0*(sin(fphi_f))/(sqrt(3.0)*(3.0-sin(fphi_f)));
        fxi_dp = 6.0*(cos(fphi_f))/(sqrt(3.0)*(3.0-sin(fphi_f)));
    }
    
    /** @brief Set plane problem
     * planestress = 1 => Plain stress state
     * planestress != 1 => Plain Strain state
     */
    void SetPlaneProblem(int planestress)
    {
        fPlaneStress = planestress;
    }
    
    /** @brief Parameters of rock and fluid: */
    void SetKModel(int model)
    {
        fk_model = model;
    }
    
    /** @brief Parameters of rock and fluid: */
    int KModel()
    {
        return fk_model;
    }
    
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    // Scalar shape functions
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

    // vectorial shape functions
    void ContributeVec(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void ContributeVec(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    void ContributeVecBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    void Compute_Sigma(TPZFMatrix<REAL> & S_eff,TPZFMatrix<REAL> & Grad_u,  REAL p_ex);
    void Compute_Sigma(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_v);
    
    REAL Inner_Product(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & T);
    
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right) {
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
    
    
    REAL k_permeability(REAL &phi, REAL &k);
    
    /** @brief Poroelastic porosity correction */
    REAL porosoty_corrected(TPZVec<TPZMaterialData> &datavec);
    

    
    
};

#endif /* TPZLinearElliptic_hpp */
