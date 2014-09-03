#ifndef TPZFRACDATA_H
#define TPZFRACDATA_H

#include "pzvec.h"

/**
 * @author Omar Duran and Nathan Shauer
 * @since 19/08/2014
 * @brief Material to store data of fracture propagation simulation
 */
class TPZFracData {
  
public:
  
  /** @brief Default Constructor */
  TPZFracData();

  /** @brief Destructor */  
  ~TPZFracData();
  
  /** @brief copy constructor */
  TPZFracData(const TPZFracData &copy);

  /** @brief operator equal */
  TPZFracData &operator=(const TPZFracData &copy);
  
private:
  
  /** @brief State: Stiffness or Mass Matrix Calculations */
  enum EState { ELastState = 0, ECurrentState = 1 };
  EState gState;
  
  /** @brief Fluid Viscosity - Pa.s */
  REAL fmu;
  
  /** @brief Simulation time step */
  REAL fDeltaT;
  
  /** @brief Simulation current time */
  REAL fTime;
  
  /** @brief Fracture length */
  REAL fLfrac;
  
  /** @brief P order of pressure (p) analysis for fracturing simulation */
  REAL fPorderPressure;

  /** @brief P order of flow (Q) analysis for fracturing simulation */
  REAL fPorderFlow;

public:
  
  /** @brief Set fluid viscosity. */
  void SetViscosity(REAL mu){this->fmu = mu;}
  
  /** @brief Returns fluid viscosity. */
  REAL Viscosity(){return this->fmu;}
  
  /** @brief Defines simulation time step. */
  void SetTimeStep(REAL timestep){ this->fDeltaT = timestep;}

  /** @brief Returns simulation time step. */
  REAL TimeStep(){return this->fDeltaT;}
  
  /** @brief Defines simulation time */
  void SetTime(REAL time){ this->fTime = time;}
  
  /** @brief Returns simulation time */
  REAL Time(){return this->fTime;}

  /** @brief Defines simulation fracture length */
  void SetLfrac(REAL Lfrac){ this->fLfrac = Lfrac;}
  
  /** @brief Returns simulation fracture length */
  REAL Lfrac(){return this->fLfrac;}

  /** @brief Defines p order of the pressure in L2 space */
  void SetPorderPressure(int PorderPressure){ this->fPorderPressure = PorderPressure;}
  
  /** @brief Returns p order of the pressure */
  int PorderPressure(){return this->fPorderPressure;}
  
  /** @brief Defines p order of the flow in H1 space */
  void SetPorderFlow(int PorderFlow){ this->fPorderFlow = PorderFlow;}
  
  /** @brief Returns p order of the flow */
  int PorderFlow(){return this->fPorderFlow;}

  /** @brief Set evaluating step n */
  void SetLastState(){ gState = ELastState;}
  
  /** @brief Set evaluating step n + 1 */
  void SetCurrentState(){ gState = ECurrentState;}
  
  bool IsLastState() { return gState == ELastState;}
  
};

#endif