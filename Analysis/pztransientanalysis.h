// -*- c++ -*-

//$Id: pztransientanalysis.h,v 1.1 2006-06-02 17:03:27 tiago Exp $

#ifndef TRANSIENTANALH
#define TRANSIENTANALH

#include "pzanalysis.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include <iostream>
class TPZCompMesh;
class TPZFMatrix;

class TPZTransientAnalysis : public TPZAnalysis {

public:

  TPZTransientAnalysis(TPZCompMesh *mesh,std::ostream &out = std::cout);
  
  ~TPZTransientAnalysis();

  virtual void Run(std::ostream &out = std::cout);
  
  /** 
   * Defines max number of steps and steady state convergence tolerance.
   */
  void SetConvergence(int niter, REAL eps);

  /** 
   * Defines max number of steps and error convergence tolerance for Newton's method.
   */  
  void SetNewtonConvergence(int niter, REAL eps);

  REAL & TimeStep();
  
  TPZFMatrix & GetSolution(int step);
  
  void SetInitialSolution(TPZFMatrix & InitialSol);
  void SetInitialSolutionAsZero();
  
  int GetCurrentIter();
    
protected:

  REAL fTimeStep;
  
  int fCurrentIter;
  
  int fNIter;
  
  REAL fSteadyTol;
  
  int fNewtonMaxIter;
  
  REAL fNewtonTol;
  
  TPZVec< TPZFMatrix > fAllSolutions;
  
  template<class T> 
  void SetLastState();
  
  void SetCurrentState();
  
  void SetAllMaterialsDeltaT();

};
 
inline void TPZTransientAnalysis::SetConvergence(int niter, REAL eps){
  this->fNIter = niter;
  this->fSteadyTol = eps;
  this->fAllSolutions.Resize(niter+1);
}

inline void TPZTransientAnalysis::SetNewtonConvergence(int niter, REAL eps){
  this->fNewtonMaxIter = niter;
  this->fNewtonTol = eps;
}

inline REAL & TPZTransientAnalysis::TimeStep(){
  return this->fTimeStep;
}

inline TPZFMatrix & TPZTransientAnalysis::GetSolution(int step){
  return this->fAllSolutions[step];
}

inline int TPZTransientAnalysis::GetCurrentIter(){
  return this->fCurrentIter;
}

#endif
