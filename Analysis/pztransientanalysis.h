/**
 * @file
 * @brief Contains TPZTransientAnalysis class which implements a simple manner to perform transient simulations.
 */

#ifndef TRANSIENTANALH
#define TRANSIENTANALH

#include "pznonlinanalysis.h"
#include "pzfmatrix.h"
#include <iosfwd>
#include <list>

class TPZCompMesh;

class TPZFStructMatrix;

/** 
 * @note It is associated to a TPZTransientMaterial< TRANSIENTCLASS >
 * @brief Implements a very simple manner to perform transient simulations. \ref analysis "Analysis"
 * @ingroup analysis
 */
/** 
 * It uses an implicit or explicit Euler scheme for time derivative
 */
template<class TRANSIENTCLASS>
class TPZTransientAnalysis : public TPZNonLinearAnalysis {
	
public:
	
	/** @brief Static attribute storing the current time of simulation */
	static double gTime;
	
	/** @brief Method for gTime attribute access */
	double GetgTime(){ return gTime; }
	
	/** @brief Constructor
	 * @param mesh for base class
	 * @param IsLinear for optimizating the process time, linear problems have flux tangent computed only once
	 * @param out for base class
	 */
	TPZTransientAnalysis(TPZCompMesh *mesh, bool IsLinear = false, std::ostream &out = std::cout);
	
	/** @brief Default destructor */
	~TPZTransientAnalysis();
	
	/** @brief Assemble flux vector and jacobian matrix */
	virtual void Assemble() override;
	
	/** @brief Executes a Newton's method for the solution of the implicit in time equation */
	virtual void RunTransient(std::ostream &out = std::cout, bool FromBegining = true, bool linesearch = true);
	
	/** @brief Solves a explicit Euler's scheme in time */
	virtual void RunExplicit(std::ostream &out = std::cout, bool FromBegining = true);
	
	/** @brief See base class for comments */  
	virtual void PostProcess(int resolution) override { TPZLinearAnalysis::PostProcess(resolution);}
	
	/** @brief See base class for comments */
	virtual void PostProcess(int resolution, int dimension) override;
	
	/** @brief See base class for comments */
	virtual void PostProcess(TPZVec<REAL> &loc, std::ostream &out = std::cout) override;
	
	/** @brief Defines max number of steps and steady state convergence tolerance. */
	void SetConvergence(int niter, REAL eps, bool ForceAllSteps = true);
	
	/** @brief Defines properties of DX file */
	void SetSaveFrequency(int SaveFrequency, int resolution);
	
	/** @brief Defines to save solution vector with SaveFrequency frequency. */
	/** If not set, no solution is kept in the process. */
	void SetSaveSolution(int SaveFrequency);
	
	/** @brief Access to saved solution. Pair of (solution vec, simulation time) */
	std::list< std::pair<TPZFMatrix<STATE>, STATE> > & GetSavedSolutions();
	
	/** @brief Defines max number of steps and error convergence tolerance for Newton's method. */  
	void SetNewtonConvergence(int niter, REAL eps);
	
	/** @brief Access to time step attribute */
	REAL & TimeStep();
	
	/** @brief Sets problem initial solution */
	void SetInitialSolution(TPZFMatrix<STATE> & InitialSol);
	
	/** @brief Sets problem initial solution as zero */
	void SetInitialSolutionAsZero();
	
	/** @brief Returns current iteration */
	int GetCurrentIter();
    
protected:
	
	/** @brief Flag indicating whether the problem is linear or not. */
	/** Linear problems require the computation and decompostition of tangent matrix only once. */
	bool fIsLinearProblem;
	
	/** @brief Simulation time step */
	REAL fTimeStep;
	
	/** @brief Current iteration. Variable allowing to restart the simulation. */
	int fCurrentIter;
	
	/** @brief Number of iterations counting from fCurrentIter to fCurrentIter+fNIter */
	int fNIter;
	
	/** @brief Tolerance to consider the problem solution as steady state */
	REAL fSteadyTol;
	
	/** @brief Flag indicating whether all steps must be performed even if tolerance is achieved. */
	bool fForceAllSteps;
	
	/** @brief Frequency which solution must be saved in DX file. */
	int fSaveFrequency;  
	/** @brief Resolution of DX mesh */
	int fDXResolution;
	
	/** @brief Frequency which solution vector must be saved. */
	/**  Zero (default value) means no solution vector but the current one is saved. */
	int fSaveSolutionVecFrequency;
	
	/** @brief Attribute to store solution vectors during process. Pair of (solution vec, simulation time) */
	/** This attribute is cleaned every time Run method is called */
	std::list< std::pair< TPZFMatrix<STATE>, STATE> > fSavedSolutionVec;
	
	/** @brief If fSaveSolutionVecFrequency != 0, save current solution vector in fSavedSolutionVec attribute. */
	void SaveCurrentSolutionVec();
	
	/** @brief Max iteration number of Newton's method */
	int fNewtonMaxIter;
	
	/** @brief Tolerance of Newton's method */
	REAL fNewtonTol;
	
	/** @brief Sets all materials in temporal scheme as an implicit Euler */
	void SetImplicit();
	
	/** @brief Sets all materials in temporal scheme as an explicit Euler */
	void SetExplicit();
	
	/** @brief Sets all materials in LastState */
	void SetLastState();
	
	/** @brief Sets all materials in CurrentState */
	void SetCurrentState();
	
	/** @brief Sets all materials to compute the mass matrix - used in the explicit scheme */
	void SetMassMatrix();
	
	/** @brief Sets all materials to compute only the flux contributions - used in the explicit scheme */
	void SetFluxOnly();
	
	/** @brief Sets all materials the time step */
	void SetAllMaterialsDeltaT();
	
	/** @brief Computes linear tangent matrix for linear problems */
	void ComputeLinearTangentMatrix();
	
	/** @brief Computes the mass matrix for the explicit scheme */
	void ComputeMassMatrix();
	
	/** @brief Computes the only the flux contribution for the explicit scheme */
	void ComputeFluxOnly();

private:
	template<class TVar>
	void AssembleInternal();
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
