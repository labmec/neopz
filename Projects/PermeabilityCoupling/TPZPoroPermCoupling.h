//
//  TPZPoroPermCoupling.hpp
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#ifndef TPZPoroPermCoupling_h
#define TPZPoroPermCoupling_h

#include <stdio.h>
#include "pzmaterial.h"
#include "pzdiscgal.h"
#include "pzvec.h"
#include "TPZSimulationData.h"
#include <iostream>


class TPZPoroPermCoupling : public TPZDiscontinuousGalerkin {
    
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
    
    /** @brief Second Lame Parameter */
    REAL fmu;
    
    /** @brief constants Biot poroelasticity */
    REAL falpha;
    
    /** @brief Storage coefficient poroelasticity */
    REAL fSe;
    
    /** @brief Permeability of the rock */
    REAL fk;
    
    /** @brief Fluid viscosity */
    REAL fvisc;
    
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
    
    TPZPoroPermCoupling();
    
    TPZPoroPermCoupling(int matid, int dim);
    
    ~TPZPoroPermCoupling();
    
    void Print(std::ostream & out);
    
    std::string Name() { return "TPZPoroPermCoupling"; }
    
    int Dimension() const {return fDim;}
    
    virtual int NStateVariables();
    
    /** @brief Parameters of rock and fluid: */
    void SetDimension(int dimension)
    {
        fDim = dimension;
    }
    
    
    /** @brief Parameters of rock and fluid: */
    void SetParameters(REAL perm, REAL visc)
    {
        fk = perm;
        fvisc = visc;
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
    
    void SetParameters(REAL l, REAL mu, REAL l_u)
    {
//        fE = (mu*(3.0*Lambda+2.0*mu))/(Lambda+mu);
//        fnu = (Lambda)/(2*(Lambda+mu));
//        
//        fEu = (mu*(3.0*Lambdau+2.0*mu))/(Lambdau+mu);
//        fnuu = (Lambdau)/(2*(Lambdau+mu));
        
        flambdau = l_u;
        flambda = l;
        fmu = mu;
    }
    
    void SetBiotParameters(REAL alpha, REAL Se)
    {
        falpha = alpha;
        fSe = Se;
    }
    
    /** @brief Set plane problem
     * planestress = 1 => Plain stress state
     * planestress != 1 => Plain Strain state
     */
    void SetPlaneProblem(int planestress)
    {
        fPlaneStress = planestress;
    }
    
    
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);

    void ContributeII(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    void ContributeBCII(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    void Compute_Sigma(TPZFMatrix<REAL> & S_eff,TPZFMatrix<REAL> & Grad_u);
    
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
};


#endif /* TPZPoroPermCoupling_hpp */
