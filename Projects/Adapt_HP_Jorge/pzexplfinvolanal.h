//$Id: pzexplfinvolanal.h,v 1.8 2009-11-04 14:13:24 fortiago Exp $

#ifndef EXPLFINVOLANALH
#define EXPLFINVOLANALH

#include "pzanalysis.h"
#include "pznonlinanalysis.h"
#include "pzcompel.h"
#include "pzinterpolationspace.h"
#include <map>
#include <list>
#include "pzfmatrix.h"
#include "pzvec.h"
#include <iostream>

class TPZCompMesh;
class TPZFMatrix;
class TPZFStructMatrix;
class TMTFaceData;

/** This class implements an explicit finite volume analysis.
 */
class TPZExplFinVolAnal : public TPZAnalysis{

public:

  TPZExplFinVolAnal(TPZCompMesh *mesh, std::ostream &out = std::cout);

  ~TPZExplFinVolAnal();

  /**
   * Assemble fluxes
  **/
  void AssembleFluxes(const TPZFMatrix & Solution, std::set<int> *MaterialIds = NULL){
    this->AssembleFluxes2ndOrder(Solution);
  }

  int NStateVariables(){
    return 5;
  }

  int Dimension(){
    return 3;
  }

  /** Computes next solution based on the last
   */
  void TimeEvolution(TPZFMatrix &LastSol, TPZFMatrix &NextSol);

  virtual void Run(std::ostream &out = std::cout);

  void MultiResolution(double Epsl, std::ostream &out = std::cout);

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

 /** Computes solution gradients 
  * Input is SolutionConsVars: conservative vars
  * Result is given in fSolution which is a vector of solution and its derivatives
  * in primitive vars
  * <!> fRhs is modified
  */
  void ComputeGradient(const TPZFMatrix & SolutionConsVars);
	
  void ComputeGradientForDetails(const TPZFMatrix & PrimitiveSolution, TPZFMatrix & SolutionWithGrad);

protected:


  void DX(int iter, std::string filename);

  /** divide vec elements by cell volume and multiply by alpha */
  void DivideByVolume(TPZFMatrix &vec, double alpha);

  /** Make loop over interfaces requesting flux computation */
  void ComputeFlux(std::list< TPZInterfaceElement* > &FacePtrList);
  void ParallelComputeFlux(std::list< TPZInterfaceElement* > &FacePtrList);
  static void * ExecuteParallelComputeFlux(void * ExtData);

  void AssembleFluxes2ndOrder(const TPZFMatrix & Solution);


  void GetNeighbourSolution(TPZInterfaceElement *face, TPZVec<REAL> &LeftSol, TPZVec<REAL> &RightSol);
  void GetSol(TPZCompElDisc * disc, TPZVec<REAL> &sol);

  /** Compute the element residual.
   * This special method for finite volume method only
   * does not use shape functions.
   * It assumes the solution vector stores the cell solution and its derivatives
   * For instance: {rho, u, p, drhodx, drhody, dudx, dudy, dpdx,dpdy }
   * Neighbour elements must be TPZCompElDisc, p = 0
   * Implemented for Olivier Roussel
   * @author Tiago Forti
   * @since 2009 Aug 18
   */
  void CalcResidualFiniteVolumeMethod(TPZInterfaceElement *face, TPZElementMatrix &ef, TPZVec<REAL> &LeftSol, TPZVec<REAL> &RightSol);

  void FromConservativeToPrimitiveAndLoad(const TPZFMatrix & Solution);

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

  /** Initialize fFacePtrList, fVolumeData, and fVecFaces class attributes */
  void InitializeAuxiliarVariables();

  /** Clean fFacePtrList and fVolumeData class attributes */
  void CleanAuxiliarVariables();

  /** Stores the list of interface elements to avoid dynamic cast all the time
   */
  std::list< TPZInterfaceElement * > fFacePtrList;

  /** Stores volume data:
   * From its pointer to the pair < cell volume, destination indices >
   */
  std::map< TPZInterpolationSpace*, std::pair< REAL, TPZVec<int> > > fVolumeData;

  /** For parallel computing */
  TPZVec< TMTFaceData * > fVecFaces;

};

struct TMTFaceData{
  TMTFaceData():fFaces(),fAn(NULL){
  }
  ~TMTFaceData(){
    fFaces.clear();
    fAn= NULL;
  }
  std::list< TPZInterfaceElement* > fFaces;
  TPZExplFinVolAnal * fAn;
};

#endif

