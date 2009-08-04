//$Id: pzeuler.h,v 1.2 2009-08-04 21:37:54 fortiago Exp $

#ifndef PZEULER_H
#define PZEULER_H

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzdiscgal.h"

#include "pzausmflux.h"
#include "pzlog.h"


/** This material implements the weak statement of the compressible euler equations
 * for Olivier Roussel's project. It is for transient analysis, finite volume method and
 * explicit time integrator only.
 */
class TPZEulerEquation : public TPZDiscontinuousGalerkin{

public:

  enum BCType{EFreeSlip = 1};

private:

  /**
   * ratio between specific heat is constant and the specific heat the constant
   * volume of a polytropic gas
   */
  REAL fGamma;

  /** Flux object */
  TPZAUSMFlux fAUSMFlux;

  /** Compute Flux
   */
  void ComputeFlux(TPZVec<REAL> &sol, TPZFMatrix & F);

public:

  /** Class constructor
   */
  TPZEulerEquation(int nummat, REAL gamma);

  /** Default destructor
   */
  ~TPZEulerEquation();

  /** Default constructor
   */
  TPZEulerEquation();

  /** Copy constructor
   */
  TPZEulerEquation(const TPZEulerEquation &cp);

  /** Creates a copy of this
   */
  TPZAutoPointer<TPZMaterial> NewMaterial();

  /**
   * Object-based overload
   */
  int NStateVariables();

  /**
   * Object-based overload
   */
  virtual int Dimension();

  /**
   * Returns the pressure value
   */
  virtual REAL Pressure(TPZVec<REAL> &U);

  /** Computes sound speed 
   */
  REAL cSpeed(TPZVec<REAL> & sol);

  /**
  * Returns u = Sqrt(u2 + v2 + w2);
  */
  REAL uRes(TPZVec<REAL> & sol);

  /**
   * See declaration in base class
   */
  virtual void Print(std::ostream & out);

  /**
   * See declaration in base class
   */
  virtual std::string Name(){return "TPZEulerEquation";}

  /**
   * See declaration in base class
   */
  virtual int VariableIndex(const std::string &name);

  /**
   * See declaration in base class
   */
  virtual int NSolutionVariables(int var);

  /**
   * See declaration in base class
   */
  virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);


  /**
   * See declaration in base class
   */
  virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

  /**
   * See declaration in base class
   * Contributes only to the rhs.
   */
  virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ef);

  /**
   * See declaration in base class
   */
  virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

  /**
   * See declaration in base class
   */
  virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef);

  /**
   * See declaration in base class
   */

  virtual void ContributeBC(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek, TPZFMatrix &ef,
                            TPZBndCond &bc);

  virtual void ContributeBCInterface(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix &ek,TPZFMatrix &ef,
                                     TPZBndCond &bc);

  virtual void ContributeBCInterface(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix &ef,
                                     TPZBndCond &bc);
};

#endif///PZEULER_H

