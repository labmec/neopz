//
//  TPZNonLinearElliptic.hpp
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#ifndef TPZNonLinearElliptic_h
#define TPZNonLinearElliptic_h

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


class TPZNonLinearElliptic : public TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> {
    
protected:
    
    /** @brief define the simulation data */
    TPZSimulationData * fSimulationData;
    
    /** @brief Problem dimension */
    int fDim;

    /** @brief Problem dimension nstate variables */
    int fnstate;

    /** @brief Parameter 0 */
    REAL fmu_0;
    
    /** @brief Parameter 1 */
    REAL fmu_1;
    
    /** @brief Parameter 2 */
    REAL fmu_2;
    
    
public:
    
    TPZNonLinearElliptic();
    
    TPZNonLinearElliptic(int matid, int dim);
    
    ~TPZNonLinearElliptic();
    
    void Print(std::ostream & out);
    
    std::string Name() { return "TPZNonLinearElliptic"; }
    
    int Dimension() const {return fDim;}
    
    virtual int NStateVariables();
    
    /** @brief Dimension of the problem */
    void SetDimension(int dimension)
    {
        fDim = dimension;
    }
    
    /** @brief Ste n state variables */
    void SetNState(int nstate)
    {
        fnstate = nstate;
    }
    
    /** @brief set parameters  */
    void SetParameters(REAL mu_0, REAL mu_1, REAL mu_2)
    {
        fmu_0 = mu_0;
        fmu_1 = mu_1;
        fmu_2 = mu_2;
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
    
    
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);

    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
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


#endif /* TPZNonLinearElliptic_hpp */
