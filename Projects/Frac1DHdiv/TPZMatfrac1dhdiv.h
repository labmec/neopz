//
//  pzmultiphase.h
//  PZ
//
//  Created by Omar Duran on 19/08/2013.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef PZ_pzmultiphase_h
#define PZ_pzmultiphase_h

#include "pzmaterial.h"
#include "pzdiscgal.h"
#ifdef _AUTODIFF
#include "fad.h"
#endif
#include <iostream>
#include <fstream>
#include <string>

/**
 * @ingroup material
 * @author Omar Duran
 * @since 19/08/2013
 * @brief Material to solve a 2d multiphase transport problem by multiphysics simulation
 * @brief Here is used L, Hdiv ... spaces and first order upwind scheme
 */


class TPZMatfrac1dhdiv : public TPZDiscontinuousGalerkin {
    
protected:
    
    /** @brief Problem dimension */
    int fDim;
    
    /** @brief State: Stiffness or Mass Matrix Calculations */
    enum EState { ELastState = 0, ECurrentState = 1 };
    EState gState;
  
    
public:
    
    TPZMatfrac1dhdiv();
    
    TPZMatfrac1dhdiv(int matid, int dim, REAL viscosity);
    
    virtual ~TPZMatfrac1dhdiv();
    
    /** @brief copy constructor */
    TPZMatfrac1dhdiv(const TPZMatfrac1dhdiv &copy);
    
    TPZMatfrac1dhdiv &operator=(const TPZMatfrac1dhdiv &copy);
    
    virtual void Print(std::ostream & out);
    
    virtual std::string Name() { return "TPZMatfrac1dhdiv"; }
    
    virtual int Dimension() const;
    
    virtual int NStateVariables();
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
    {
        DebugStop(); // Should never be called for this material
    }
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        DebugStop();
    }
    
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    
    virtual int VariableIndex(const std::string &name);
    
    virtual int NSolutionVariables(int var);
    
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    
private:
    
    /** @brief Fluid Viscosity - Pa.s */
    REAL fmu;
    
    /** @brief Simulation time step */
    REAL fDeltaT;
    
    /** @brief Simulation current time */
    REAL fTime;
    
    /** @brief Parameter representing temporal scheme for transport equation */
    REAL fTheta;
    
public:
    
    /** @brief Set fluid viscosity. */
    void SetViscosity(REAL mu){this->fmu = mu;}
    
    /** @brief Defines simulation time step. */
    void SetTimeStep(REAL timestep){ this->fDeltaT = timestep;}
    
    /** @brief Defines simulation time step. */
    void SetTime(REAL time){ this->fTime = time;}
    
    /** @brief Defines stemporal scheme. */
    void SetTScheme(REAL timetheta){ this->fTheta = timetheta;}

    /** @brief Set evaluating step n */
    void SetLastState(){ gState = ELastState;}

    /** @brief Set evaluating step n + 1 */
    void SetCurrentState(){ gState = ECurrentState;}
    
    /** @brief Apply dirichlet condition on flux strongly */
    void ApplyQnD(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
    
    /** @brief Apply neumann condition on pressure weakly */
    void ApplyPN(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
    
    
};

#endif