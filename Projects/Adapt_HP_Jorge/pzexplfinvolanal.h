/**
 * @file
 * @brief Contains the declaration of TPZExplFinVolAnal class. (Explicit Finite volume method)
 */

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

template<class TVar>
class TPZFMatrix;
class TPZFStructMatrix;
class TMTFaceData;

/**
 * @brief This class implements an explicit finite volume analysis.
 */
class TPZExplFinVolAnal : public TPZAnalysis{

public:

  TPZExplFinVolAnal(TPZCompMesh *mesh, std::ostream &out = std::cout);

  ~TPZExplFinVolAnal();

  /** @brief Assemble fluxes */
  void AssembleFluxes(const TPZFMatrix<REAL> & Solution, std::set<int> *MaterialIds = NULL) {
    this->AssembleFluxes2ndOrder(Solution);
  }
	/** @brief Number of the variables */
  int NStateVariables() {
    return 5;
  }

	/** @brief Model dimension */
  int Dimension() {
    return 3;
  }

  /** @brief Computes next solution based on the last */
  void TimeEvolution(TPZFMatrix<REAL> &LastSol, TPZFMatrix<REAL> &NextSol);

  virtual void Run(std::ostream &out = std::cout);

  void MultiResolution(double Epsl, std::ostream &out = std::cout);

  virtual void PostProcess(int resolution){ 
    TPZAnalysis::PostProcess(resolution);
  }

  virtual void PostProcess(int resolution, int dimension);

  /** @brief Defines max number of steps and steady state convergence tolerance. */
  void Set(REAL timestep, int niter, REAL eps, bool ForceAllSteps = true);

  /** @brief Defines properties of DX file */
  void SetSaveFrequency(int SaveFrequency, int resolution);

  REAL TimeStep();

  void SetInitialSolution(TPZFMatrix<REAL> & InitialSol);

  void SetInitialSolutionAsZero();

 /**
  * @brief Computes solution gradients 
  * @param SolutionConsVars [in] conservative variables
  * @return Result is given in fSolution which is a vector of solution and its derivatives in primitive vars
  * @note fRhs is modified
  */
  void ComputeGradient(const TPZFMatrix<REAL> & SolutionConsVars);
	
  void ComputeGradientForDetails(const TPZFMatrix<REAL> & PrimitiveSolution, TPZFMatrix<REAL> & SolutionWithGrad);

protected:


  void DX(int iter, std::string filename);

  /** @brief Divide vec elements by cell volume and multiply by alpha */
  void DivideByVolume(TPZFMatrix<REAL> &vec, double alpha);

  /** @brief Make loop over interfaces requesting flux computation */
  void ComputeFlux(std::list< TPZInterfaceElement* > &FacePtrList);
  void ParallelComputeFlux(std::list< TPZInterfaceElement* > &FacePtrList);
  static void * ExecuteParallelComputeFlux(void * ExtData);

  void AssembleFluxes2ndOrder(const TPZFMatrix<REAL> & Solution);


  void GetNeighbourSolution(TPZInterfaceElement *face, TPZVec<REAL> &LeftSol, TPZVec<REAL> &RightSol);
  void GetSol(TPZCompElDisc * disc, TPZVec<REAL> &sol);

  /**
   * @brief Compute the element residual. This special method for finite volume method only does not use shape functions.
   * @note It assumes the solution vector stores the cell solution and its derivatives. For instance: {rho, u, p, drhodx, drhody, dudx, dudy, dpdx,dpdy } \n
   * Neighbour elements must be TPZCompElDisc, p = 0
   * @author Olivier Roussel
   * @author Tiago Forti
   * @since 2009 Aug 18
   */
  void CalcResidualFiniteVolumeMethod(TPZInterfaceElement *face, TPZElementMatrix &ef, TPZVec<REAL> &LeftSol, TPZVec<REAL> &RightSol);

  void FromConservativeToPrimitiveAndLoad(const TPZFMatrix<REAL> & Solution);

  /** @brief Simulation time step */
  REAL fTimeStep;

  /** @brief Number of iterations counting from fCurrentIter to fCurrentIter+fNIter */
  int fNMaxIter;

  /** @brief Auxiliar variable only for DX post-processing */
  REAL fSimulationTime;

  /** @brief Tolerance to consider the problem solution as steady state */
  REAL fSteadyTol;

  /** @brief Flag indicating whether all steps must be performed even if tolerance is achieved. */
  bool fForceAllSteps;

  /** @brief Frequency which solution must be saved in DX file. */
  int fSaveFrequency;

  /** @brief Resolution of DX mesh */
  int fDXResolution;

  /** @brief Initialize fFacePtrList, fVolumeData, and fVecFaces class attributes */
  void InitializeAuxiliarVariables();

  /** @brief Clean fFacePtrList and fVolumeData class attributes */
  void CleanAuxiliarVariables();

  /** @brief Stores the list of interface elements to avoid dynamic cast all the time */
  std::list< TPZInterfaceElement * > fFacePtrList;

  /** @brief Stores volume data: From its pointer to the pair < cell volume, destination indices > */
  std::map< TPZInterpolationSpace*, std::pair< REAL, TPZVec<int> > > fVolumeData;

  /** @brief For parallel computing */
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

