//$Id: pzeuler.h,v 1.3 2009-08-28 19:43:43 fortiago Exp $

#ifndef PZEULER_H
#define PZEULER_H

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzdiscgal.h"

#include "pzausmflux.h"
#include "pzgradientflux.h"
#include "pzlog.h"

#define LinearConvection

/** This material implements the weak statement of the compressible euler equations
 * for Olivier Roussel's project. It is for transient analysis, finite volume method and
 * explicit time integrator only.
 */
class TPZEulerEquation : public TPZDiscontinuousGalerkin{

public:

  enum BCType{EFreeSlip = 1};

  enum CALCType{EFlux = 1, EGradient = 2};

  static void SetComputeFlux(){
    gType = EFlux;
  }

  static void SetComputeGradient(){
    gType = EGradient;
  }

#ifdef LinearConvection
  static void SetLinearConvection(TPZVec<REAL> &Celerity);
#endif

private:

#ifdef LinearConvection
  static TPZVec<REAL> gCelerity;
#endif

  static CALCType gType;

  /**
   * ratio between specific heat is constant and the specific heat the constant
   * volume of a polytropic gas
   */
  static REAL gGamma;

  /** Convective flux object */
  TPZAUSMFlux fAUSMFlux;

  /** Gradient flux object */
  TPZGradientFlux fGradientFlux;

  /** Compute Euler Flux
   */
  void ComputeEulerFlux(TPZVec<REAL> &sol, TPZFMatrix & F);

public:

  static REAL Gamma(){ return gGamma; }

  /** Convert from primitive to conservative variables */
  static void FromPrimitiveToConservative(TPZVec<REAL> &sol,REAL gamma);

  /** Convert from conservative to primitive variables */
  static void FromConservativeToPrimitive(TPZVec<REAL> &sol,REAL gamma);

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
  static REAL Pressure(TPZVec<REAL> &U, double gamma);

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
   * data.soll and data.solr are expected in primitive variables
   */
  virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

  /**
   * See declaration in base class
   * data.soll and data.solr are expected in primitive variables
   */
  virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef);

  /**
   * See declaration in base class
   */

  virtual void ContributeBC(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek, TPZFMatrix &ef,
                            TPZBndCond &bc);

  /**
   * See declaration in base class
   * data.soll and data.solr are expected in primitive variables
   */
  virtual void ContributeBCInterface(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix &ek,TPZFMatrix &ef,
                                     TPZBndCond &bc);

  /**
   * See declaration in base class
   * data.soll and data.solr are expected in primitive variables
   */
  virtual void ContributeBCInterface(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix &ef,
                                     TPZBndCond &bc);
};

#endif///PZEULER_H

