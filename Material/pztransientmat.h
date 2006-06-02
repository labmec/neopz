// -*- c++ -*-

//$Id: pztransientmat.h,v 1.1 2006-06-02 17:03:59 tiago Exp $


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
  
  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
                          TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);
                          
  virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
                            TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

 /**
  * Set material to compute only Integral[- un/deltaT * v, Omega]
  */
  void SetLastState();
  
 /**
  * Set material to compute Integral[un+1/deltaT * v, Omega] + Bilinear Form = Linear Form 
  */
  void SetCurrentState();
  
  /** 
   * Define time step DeltaT
   */
  void SetTimeStep(REAL TimeStep);
  
  /**
   * Returns time step value.
   */
  REAL TimeStep();
  
  protected:
  
  enum STEPS{ENone = -1, ELast = 0, ECurrent = 1};
  
  STEPS fStep;
  
  REAL fTimeStep;
  
  void ContributeSolutionRhs(TPZVec<REAL> &sol, TPZFMatrix &phi, REAL weight, TPZFMatrix &ef);
  
  void ContributeTangent(TPZFMatrix &phi, REAL weight, TPZFMatrix &ek);
};

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetLastState(){
  this->fStep = ELast;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetCurrentState(){
  this->fStep = ECurrent;
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
