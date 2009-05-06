// -*- c++ -*-

//$Id: pztransientanalysis.h,v 1.5 2009-05-06 20:13:37 fortiago Exp $

#ifndef TRANSIENTANALH
#define TRANSIENTANALH

#include "pzanalysis.h"
#include "pznonlinanalysis.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include <iostream>
#include <list>

class TPZCompMesh;
class TPZFMatrix;
class TPZFStructMatrix;

/** Implements a very simple manner to perform transient simulations
 * It uses an implicit or explicit Euler scheme for time derivative
 * It is associated to a TPZTransientMaterial< TRANSIENTCLASS > class
 */
template<class TRANSIENTCLASS>
class TPZTransientAnalysis : public TPZNonLinearAnalysis {

public:

  /** Static attribute storing the current time of simulation
   */
  static double gTime;

  double GetgTime(){ return gTime; }

  /** Class constructor
   * @param mesh for base class
   * @param IsLinear for optimizating the process time, linear problems have flux tangent computed only once
   * @param out for base class
   */
  TPZTransientAnalysis(TPZCompMesh *mesh, bool IsLinear = false, std::ostream &out = std::cout);
  
  /** Default destructor
   */
  ~TPZTransientAnalysis();
  
  /**
   * Assemble flux vector and jacobian matrix
   */
  virtual void Assemble();

  /** 
   * Executes a Newton's method for the solution of the implicit in time equation 
   */
  virtual void Run(std::ostream &out = std::cout, bool FromBegining = true, bool linesearch = true);

  /** 
   * Solves a explicit Euler's scheme in time
   */
  virtual void RunExplicit(std::ostream &out = std::cout, bool FromBegining = true);

  /** See base class for comments
  */  
  virtual void PostProcess(int resolution){ TPZAnalysis::PostProcess(resolution);}
  
  /** See base class for comments
  */
  virtual void PostProcess(int resolution, int dimension);
  
  /** See base class for comments
  */
  virtual void PostProcess(TPZVec<REAL> &loc, std::ostream &out = std::cout);
  
  /** 
   * Defines max number of steps and steady state convergence tolerance.
   */
  void SetConvergence(int niter, REAL eps, bool ForceAllSteps = true);

  /** Defines properties of DX file */
  void SetSaveFrequency(int SaveFrequency, int resolution);

  /** Defines to save solution vector with SaveFrequency frequency.
   * If not set, no solution is kept in the process.
   */
  void SetSaveSolution(int SaveFrequency);

  /** Access to saved solution. Pair of (solution vec, simulation time)
    */
  std::list< std::pair<TPZFMatrix, REAL> > & GetSavedSolutions();

  /** 
   * Defines max number of steps and error convergence tolerance for Newton's method.
   */  
  void SetNewtonConvergence(int niter, REAL eps);

  /** Access to time step attribute
  */
  REAL & TimeStep();

  /** Set problem initial solution
   */
  void SetInitialSolution(TPZFMatrix & InitialSol);

  /** Set problem initial solution as zero
   */
  void SetInitialSolutionAsZero();
  
  /** Returns current iteration
   */
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

  /** Frequency which solution vector must be saved.
   *  Zero (default value) means no solution vector but the current one is saved.
   */
  int fSaveSolutionVecFrequency;

  /** Attribute to store solution vectors during process. Pair of (solution vec, simulation time)
   * This attribute is cleaned every time Run method is called
   */
  std::list< std::pair< TPZFMatrix, REAL> > fSavedSolutionVec;

  /** If fSaveSolutionVecFrequency != 0, save current solution vector in fSavedSolutionVec attribute. */
  void SaveCurrentSolutionVec();

  /** Max iteration number of Newton's method */
  int fNewtonMaxIter;
  
  /** Tolerance of Newton's method */
  REAL fNewtonTol;

  /** Sel all materials in temporal scheme as an implicit Euler */
  void SetImplicit();

  /** Sel all materials in temporal scheme as an explicit Euler */
  void SetExplicit();

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
