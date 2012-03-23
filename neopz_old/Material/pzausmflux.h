//$Id: pzausmflux.h,v 1.2 2009-08-28 19:43:43 fortiago Exp $

#ifndef TPZAUSMFLUX_H
#define TPZAUSMFLUX_H

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzreal.h"

class TPZAUSMFlux{

  private:

   /**
   * ratio between specific heat is constant and the specific heat the constant
   * volume of a polytropic gas
   */
  REAL fGamma;

  /** Method constants
   */
  REAL fAlpha, fBeta;

  public:

  TPZAUSMFlux(REAL gamma);

  TPZAUSMFlux(const TPZAUSMFlux &cp);

  /** Computes numerical flux
   */
  void ComputeFlux(TPZVec<REAL> &solL, TPZVec<REAL> &solR, TPZVec<REAL> &normal, TPZVec<REAL> & F);

  private:

  /** Returns sound speed
   */
  REAL SoundSpeed(TPZVec<REAL> &sol, REAL press);

  /** Returns pressure values
   */
  REAL Pressure(TPZVec<REAL> &sol);

  /** Returns speed
   */
  REAL Speed(TPZVec<REAL> &sol, TPZVec<REAL> &normal, REAL &NormalSpeed);

  /** Returns enthalpy
   */
  REAL Enthalpy(REAL soundSpeed, REAL speed);

  /** Auxiliar method only
  */
  void ComputeInitialData(TPZVec<REAL>&sol,TPZVec<REAL> &normal, REAL&soundSpeed, 
                          REAL &Speed, REAL &NormalSpeed, REAL &Enthalpy, REAL &press);

  /** Returns pressure in the face */
  REAL FacePressure(REAL pL, REAL pR, REAL Ml, REAL Mr);

  /** Returns mach number in the face */
  REAL FaceMachNumber(REAL Ml, REAL Mr);

  /** Computes the numerical sound speed at the face
  */
  REAL NumSoundSpeed(REAL LeftSoundSpeed,REAL RightSoundSpeed);

  /** Returns the mass flux
  */
  REAL MassFlux(REAL NumericalSoundSpeed, REAL rhoL, REAL rhoR, REAL FaceMach);


};///class

#endif
