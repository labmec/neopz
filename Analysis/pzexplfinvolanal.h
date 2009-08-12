//$Id: pzexplfinvolanal.h,v 1.3 2009-08-12 21:04:57 fortiago Exp $

#ifndef EXPLFINVOLANALH
#define EXPLFINVOLANALH

#include "pzanalysis.h"
#include "pznonlinanalysis.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include <iostream>

class TPZCompMesh;
class TPZFMatrix;
class TPZFStructMatrix;

/** This class implements an explicit finite volume analysis.
 */
class TPZExplFinVolAnal : public TPZAnalysis{

public:

  TPZExplFinVolAnal(TPZCompMesh *mesh, std::ostream &out = std::cout);

  ~TPZExplFinVolAnal();

  /**
   * Assemble fluxes
  **/
  void AssembleFluxes(const TPZFMatrix & Solution, std::set<int> *MaterialIds = NULL);

  /** Computes next solution based on the last
   */
  void UpdateSolution(TPZFMatrix &LastSol, TPZFMatrix &NextSol);

  virtual void Run(std::ostream &out = std::cout);

  virtual void PostProcess(int resolution){ 
    TPZAnalysis::PostProcess(resolution);
  }

  virtual void PostProcess(int resolution, int dimension);

  /**
   * Defines max number of steps and steady state convergence tolerance.
   */
  void Set(REAL timestep, int niter, REAL eps, bool ForceAllSteps = true);

  /** Defines properties of DX file
   */
  void SetSaveFrequency(int SaveFrequency, int resolution);

  REAL TimeStep();

  void SetInitialSolution(TPZFMatrix & InitialSol);

  void SetInitialSolutionAsZero();

protected:

  /** Simulation time step */
  REAL fTimeStep;

  /** Number of iterations counting from fCurrentIter to fCurrentIter+fNIter */
  int fNMaxIter;

  /** Auxiliar variable only for DX post-processing
   */
  REAL fSimulationTime;

  /** Tolerance to consider the problem solution as steady state */
  REAL fSteadyTol;

  /** Flag indicating whether all steps must be performed even if tolerance is achieved. */
  bool fForceAllSteps;

  /** Frequency which solution must be saved in DX file. */
  int fSaveFrequency;

  /** Resolution of DX mesh */
  int fDXResolution;

};

#endif

