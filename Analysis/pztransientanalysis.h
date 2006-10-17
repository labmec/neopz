// -*- c++ -*-

//$Id: pztransientanalysis.h,v 1.3 2006-10-17 02:00:58 phil Exp $

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

template<class TRANSIENTCLASS>
class TPZTransientAnalysis : public TPZAnalysis {

public:

  static double gTime;

  TPZTransientAnalysis(TPZCompMesh *mesh, bool IsLinear = false, std::ostream &out = std::cout);
  
  ~TPZTransientAnalysis();
  
    /**
   *Assemble the stiffness matrix
   **/
  virtual void Assemble();

  virtual void Run(std::ostream &out = std::cout, bool FromBegining = true);
  virtual void RunExplicit(std::ostream &out = std::cout, bool FromBegining = true);
  
  virtual void PostProcess(int resolution){ TPZAnalysis::PostProcess(resolution);}
  
  virtual void PostProcess(int resolution, int dimension);
  
  virtual void PostProcess(TPZVec<REAL> &loc, std::ostream &out = std::cout);
  
  /** 
   * Defines max number of steps and steady state convergence tolerance.
   */
  void SetConvergence(int niter, REAL eps, bool ForceAllSteps = true);

  /** Defines properties of DX file */
  void SetSaveFrequency(int SaveFrequency, int resolution);

  /** 
   * Defines max number of steps and error convergence tolerance for Newton's method.
   */  
  void SetNewtonConvergence(int niter, REAL eps);

  REAL & TimeStep();
  
  void SetInitialSolution(TPZFMatrix & InitialSol);
  void SetInitialSolutionAsZero();
  
  int GetCurrentIter();
    
protected:

  /** Flag indicating whether the problem is linear or not. 
   * Linear problems require the computation and decompostition of tangent
   * matrix only once.
   */
  bool fIsLinearProblem;

  /** Simulation time step */
  REAL fTimeStep;
  
  /** Current iteration. Variable allowing to restart the simulation. */
  int fCurrentIter;
  
  /** Number of iterations counting from fCurrentIter to fCurrentIter+fNIter */
  int fNIter;
  
  /** Tolerance to consider the problem solution as steady state */
  REAL fSteadyTol;
  
  /** Flag indicating whether all steps must be performed even if tolerance is achieved. */
  bool fForceAllSteps;
  
  /** Frequency which solution must be saved in DX file. */
  int fSaveFrequency;  
  /** Resolution of DX mesh */
  int fDXResolution;
   
  /** Max iteration number of Newton's method */
  int fNewtonMaxIter;
  
  /** Tolerance of Newton's method */
  REAL fNewtonTol;
  
  /** Sel all materials in LastState */
  void SetLastState();
  
  /** Sel all materials in CurrentState */
  void SetCurrentState();
  
  /** Set all materials to compute the mass matrix - used in the explicit scheme */
  void SetMassMatrix();
  
  /** Set all materials to compute only the flux contributions - used in the explicit scheme */
  void SetFluxOnly();
  
  /** Set all materials the time step */
  void SetAllMaterialsDeltaT();
  
  /** Compute linear tangent matrix for linear problems */
  void ComputeLinearTangentMatrix();
  
  /** Compute the mass matrix for the explicit scheme */
  void ComputeMassMatrix();
  
  /** Compute the only the flux contribution for the explicit scheme */
  void ComputeFluxOnly();

};

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis< TRANSIENTCLASS >::SetConvergence(int niter, REAL eps, bool ForceAllSteps){
  this->fNIter = niter;
  this->fSteadyTol = eps;
  this->fForceAllSteps = ForceAllSteps;
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis< TRANSIENTCLASS >::SetSaveFrequency(int SaveFrequency, int resolution){
  this->fSaveFrequency = SaveFrequency;
  this->fDXResolution = resolution;
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis< TRANSIENTCLASS >::SetNewtonConvergence(int niter, REAL eps){
  this->fNewtonMaxIter = niter;
  this->fNewtonTol = eps;
}

template<class TRANSIENTCLASS>
inline REAL & TPZTransientAnalysis< TRANSIENTCLASS >::TimeStep(){
  return this->fTimeStep;
}

template<class TRANSIENTCLASS>
inline int TPZTransientAnalysis< TRANSIENTCLASS >::GetCurrentIter(){
  return this->fCurrentIter;
}

#endif
