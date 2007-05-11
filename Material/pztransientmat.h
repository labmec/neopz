// -*- c++ -*-

//$Id: pztransientmat.h,v 1.5 2007-05-11 19:15:18 joao Exp $


#ifndef TRANSIENTMATH
#define TRANSIENTMATH

#include "pzmaterial.h"
// #include "pzvec.h"
// #include "pzfmatrix.h"

/**
 * Class TPZTransient material implements an implicit Euler time integrator.
 * Weak statement is supposed to be Integral[(un+1 - un)/deltaT * v, Omega] + Bilinear Form = Linear Form
 * This class implements only Integral[(un+1 - un)/deltaT * v, Omega]. Bilinear and linear form must be implemented in base class TBASEMAT.
 */
template<class TBASEMAT>
class TPZTransientMaterial : public TBASEMAT {

  public:
  
  TPZTransientMaterial(int nummat, int dim, REAL TimeStep);
  
  ~TPZTransientMaterial();
  
  TPZTransientMaterial(const TPZTransientMaterial &cp);
  
  /** Set integral scheme as an explicit Euler */
  void SetExplicit();
  
  /** Set integral scheme as an implicit Euler */
  void SetImplicit();
  
  virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
                            TPZFMatrix &ef);
                          
  virtual void ContributeBC(TPZMaterialData &data,
                             REAL weight,
                             TPZFMatrix &ek,
                             TPZFMatrix &ef,
                             TPZBndCond &bc);
  
  virtual void ContributeInterface(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix &ek,
                                     TPZFMatrix &ef);

  virtual void ContributeBCInterface(TPZMaterialData &data,
                                       REAL weight,
                                       TPZFMatrix &ek,
                                       TPZFMatrix &ef,
                                       TPZBndCond &bc);

 /**
  * Set material to compute only Integral[- un/deltaT * v, Omega]
  */
  void SetLastState();
  
 /**
  * Set material to compute Integral[un+1/deltaT * v, Omega] + Bilinear Form = Linear Form 
  */
  void SetCurrentState();
  
  /**
  * Set material to compute ek = Integral[phi_i phi_j, Omega]/deltaT
  */
  void SetMassMatrix();
  
  /**
  * Set material to compute ef = Linear Form - Bilinear Form(u) = F -ku
  */ 
  void SetFluxOnly();
  
  /** 
   * Define time step DeltaT
   */
  void SetTimeStep(REAL TimeStep);
  
  /**
   * Returns time step value.
   */
  REAL TimeStep();
  
  /** Indicates if the material requires the solution to compute Contribute
   * By default its value is true, but it can be set as false by derived material classes
   * to increase the performance of method TPZCompEl::CalcStiff
   */
  virtual bool NeedsSolutionToContribute(){
    return true;
  }

  /** Indicates if the material requires the global coordinate X to compute Contribute
   * By default its value is true, but it can be set as false by derived material classes
   * to increase the performance of method TPZCompEl::CalcStiff
   */
  virtual bool NeedsXCoord(){
    return (this->fStep != ELast);
  }
  
  protected:
  
  enum ETemporalScheme{EImplicit = 1, EExplicit = 2};
  
  static int gTemporalIntegrator;
  
  enum STEPS{ENone = -1, ELast = 0, ECurrent = 1, EMassMatrix = 2, EFluxOnly = 3};
  
  STEPS fStep;
  
  REAL fTimeStep;
  
  void ContributeSolutionRhs(TPZVec<REAL> &sol, TPZFMatrix &phi, REAL weight, TPZFMatrix &ef);
  
  void ContributeTangent(TPZFMatrix &phi, REAL weight, TPZFMatrix &ek);
};

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetLastState(){
  this->SetImplicit();
  this->fStep = ELast;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetCurrentState(){
  this->SetImplicit();
  this->fStep = ECurrent;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetMassMatrix(){
  this->SetExplicit();
  this->fStep = EMassMatrix;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetFluxOnly(){
  this->SetExplicit();
  this->fStep = EFluxOnly;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetTimeStep(REAL TimeStep){
  this->fTimeStep = TimeStep;
}
  
template<class TBASEMAT>
inline REAL TPZTransientMaterial< TBASEMAT >::TimeStep(){
  return this->fTimeStep;
}
  
#endif
