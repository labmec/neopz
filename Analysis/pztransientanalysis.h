// -*- c++ -*-

//$Id: pztransientanalysis.h,v 1.2 2006-07-06 15:59:09 tiago Exp $

#ifndef TRANSIENTANALH
#define TRANSIENTANALH

#include "pzanalysis.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include <iostream>

class TPZCompMesh;
class TPZFMatrix;
class TPZFStructMatrix;

class TPZTransientAnalysis : public TPZAnalysis {

public:

  static double gTime;

  TPZTransientAnalysis(TPZCompMesh *mesh, bool IsLinear = false, std::ostream &out = std::cout);
  
  ~TPZTransientAnalysis();
  
    /**
   *Assemble the stiffness matrix
   **/
  virtual void Assemble();

  virtual void Run(std::ostream &out = std::cout);
  
  virtual void PostProcess(int resolution){ TPZAnalysis::PostProcess(resolution);}
  
  virtual void PostProcess(int resolution, int dimension);
  
  virtual void PostProcess(TPZVec<REAL> &loc, std::ostream &out = std::cout);
  
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

  bool fIsLinearProblem;

  REAL fTimeStep;
  
  int fCurrentIter;
  
  int fNIter;
  
  REAL fSteadyTol;
  
  int fNewtonMaxIter;
  
  REAL fNewtonTol;
  
  TPZVec< TPZFMatrix > fAllSolutions;
  
  TPZFMatrix fLastState;
  
  void SetLastState();
  
  void SetCurrentState();
  
  void SetAllMaterialsDeltaT();
  
  void TangentResidual();
  
  void ComputeLinearTangentMatrix();
  
  

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
