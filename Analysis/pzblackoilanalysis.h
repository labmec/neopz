//$Id: pzblackoilanalysis.h,v 1.6 2011-04-05 19:32:55 calle Exp $

#ifndef BLACKOILANALH
#define BLACKOILANALH

#include "pzanalysis.h"
#include "pznonlinanalysis.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include <iostream>

class TPZCompMesh;
class TPZFMatrix;
class TPZFStructMatrix;

#ifdef _AUTODIFF
/** Jorge
 * @brief Derived from TPZNonLinearAnalysis class
 * @ingroup analysis
 */
class TPZBlackOilAnalysis : public TPZNonLinearAnalysis {

private:

  TPZFMatrix fLastState;

public:

  TPZBlackOilAnalysis(TPZCompMesh *mesh, double TimeStep, std::ostream &out = std::cout);

  ~TPZBlackOilAnalysis();

  /**
   * Assemble residual vector and tangent matrix
  */
  virtual void Assemble();

  /**
  * Assemble only the residual vector
  **/
  virtual void AssembleResidual();

  /**
   * Invert the algebraic system
  **/
  virtual void Solve();

  virtual void Run(std::ostream &out = std::cout, bool linesearch = true);

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

  REAL &TimeStep();

  void SetInitialSolution(TPZFMatrix & InitialSol);
  void SetInitialSolutionAsZero();

protected:

  /** Simulation time step */
  REAL fTimeStep;

  /** */
  REAL fSimulationTime;

  /** Number of iterations counting from fCurrentIter to fCurrentIter+fNIter */
  int fNIter;

  /** Local variable indicating the current step of simulation */
  int fCurrentStep;

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

  /** Set all materials the time step */
  void SetAllMaterialsDeltaT();
  
  static double PressaoMedia(TPZBlackOilAnalysis &an, int matid);
  static void Vazao(TPZBlackOilAnalysis &an, int matid, double & VazaoAguaSC, double  & VazaoOleoSC, double & VazaoAguaFundo, double  & VazaoOleoFundo);
};

#endif

#endif
