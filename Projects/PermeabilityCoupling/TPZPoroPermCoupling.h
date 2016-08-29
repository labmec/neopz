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
    
    /** @brief Intact rock porosity */
    REAL fporosity_0;
    
    /** @brief Permeability of the rock */
    REAL fk;
    
    /** @brief Fluid viscosity */
    REAL feta;
    
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

    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
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
    
    /** @brief Rudnicki diffusion coefficient */
    /** J. W. Rudnicki. Fluid mass sources and point forces in linear elastic di usive solids. Journal of Mechanics of Materials, 5:383–393, 1986. */
    REAL c_diffusion(REAL phi);
    
    /** @brief Poroelastic porosity correction */
    REAL porosoty_corrected(TPZVec<TPZMaterialData> &datavec);
    
};


#endif /* TPZPoroPermCoupling_hpp */
