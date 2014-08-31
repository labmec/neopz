//
//  pzmultiphase.h
//  PZ
//
//  Created by Omar Duran and Nathan Shauer on 19/08/2014.
//  Copyright (c) 2014 LabMec-Unicamp. All rights reserved.
//

#ifndef PZ_pzmultiphase_h
#define PZ_pzmultiphase_h

#include "pzmaterial.h"
#ifdef _AUTODIFF
#include "fad.h"
#endif
#include <iostream>
#include <fstream>
#include <string>

/**
 * @ingroup material
 * @author Omar Duran and Nathan Shauer
 * @since 19/08/2013
 * @brief Material to solve a 1d mixed formulation for fracture opening
 * @brief Here is used H1 for flux and L2 for pressure
 */
class TPZMatfrac1dhdiv : public TPZMaterial {
  
protected:
  
  /** @brief Problem dimension */
  int fDim;
  
  /** @brief State: Stiffness or Mass Matrix Calculations */
  enum EState { ELastState = 0, ECurrentState = 1 };
  EState gState;
  
  
public:
  
  /** @brief Default Constructor */
  TPZMatfrac1dhdiv();
  
  /** @brief Constructor with matid */
  TPZMatfrac1dhdiv(int matid);

  /** @brief Destructor */  
  virtual ~TPZMatfrac1dhdiv();
  
  /** @brief copy constructor */
  TPZMatfrac1dhdiv(const TPZMatfrac1dhdiv &copy);

  /** @brief operator equal */
  TPZMatfrac1dhdiv &operator=(const TPZMatfrac1dhdiv &copy);
  
  /** @brief Print Method */
  virtual void Print(std::ostream & out);
  
  /** @brief Name of the material */
  virtual std::string Name() { return "TPZMatfrac1dhdiv"; }
  
  /** @brief Returns the integrable dimension */
  virtual int Dimension() const;
  
  /** @brief Return the number of state variables */
  virtual int NStateVariables();
  
  /** @brief Contribute method for not multiphysics materials */
  virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
  {
    DebugStop();
  }
  
  /** @brief Contribute method beeing used */
  virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
  
  /** @brief ContributeBC method for not multiphysics materials */
  virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
  {
    DebugStop();
  }

  /** @brief ContributeBC method beeing used */
  virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
  
  /** @brief Fill material data parameter with necessary requirements for the Contribute method*/
  virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);

  /** @brief Fill material data parameter with necessary requirements for the ContributeBC method*/
  virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
  
  /** @brief Returns the variable index associated with the name */
  virtual int VariableIndex(const std::string &name);
  
  /** @brief Returns the number of variables associated with the variable indexed by var */
  virtual int NSolutionVariables(int var);
  
  /** @brief Calculates a solution given datavec*/
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
  
};

#endif