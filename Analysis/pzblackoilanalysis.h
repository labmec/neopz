/**
 * @file
 * @brief Contains TPZBlackOilAnalysis class derived from TPZNonLinearAnalysis class.
 */
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

template <class TVar>
class TPZFMatrix;
class TPZFStructMatrix;

#ifdef _AUTODIFF

/**
 * @brief Derived from TPZNonLinearAnalysis class. \ref analysis "Analysis"
 * @ingroup analysis
 */
class TPZBlackOilAnalysis : public TPZNonLinearAnalysis {
	
private:
	
	TPZFMatrix<REAL> fLastState;
	
public:
	
	TPZBlackOilAnalysis(TPZCompMesh *mesh, double TimeStep, std::ostream &out = std::cout);
	
	~TPZBlackOilAnalysis();
	
	/**
	 * @brief Assemble residual vector and tangent matrix
	 */
	virtual void Assemble();
	
	/**
	 * @brief Assemble only the residual vector
	 **/
	virtual void AssembleResidual();
	
	/**
	 * @brief Invert the algebraic system
	 **/
	virtual void Solve();
	
	virtual void Run(std::ostream &out = std::cout, bool linesearch = true);
	
	virtual void PostProcess(int resolution){ TPZAnalysis::PostProcess(resolution);}
	
	virtual void PostProcess(int resolution, int dimension);
	
	virtual void PostProcess(TPZVec<REAL> &loc, std::ostream &out = std::cout);
	
	/**
	 * @brief Defines max number of steps and steady state convergence tolerance.
	 */
	void SetConvergence(int niter, REAL eps, bool ForceAllSteps = true);
	
	/** @brief Defines properties of DX file */
	void SetSaveFrequency(int SaveFrequency, int resolution);
	
	/**
	 * @brief Defines max number of steps and error convergence tolerance for Newton's method.
	 */
	void SetNewtonConvergence(int niter, REAL eps);
	
	REAL &TimeStep();
	
	void SetInitialSolution(TPZFMatrix<REAL> & InitialSol);
	void SetInitialSolutionAsZero();
	
protected:
	
	/** @brief Simulation time step */
	REAL fTimeStep;
	
	/** */
	REAL fSimulationTime;
	
	/** @brief Number of iterations counting from fCurrentIter to fCurrentIter+fNIter */
	int fNIter;
	
	/** @brief Local variable indicating the current step of simulation */
	int fCurrentStep;
	
	/** @brief Tolerance to consider the problem solution as steady state */
	REAL fSteadyTol;
	
	/** @brief Flag indicating whether all steps must be performed even if tolerance is achieved. */
	bool fForceAllSteps;
	
	/** @brief Frequency which solution must be saved in DX file. */
	int fSaveFrequency;
	
	/** @brief Resolution of DX mesh */
	int fDXResolution;
	
	/** @brief Max iteration number of Newton's method */
	int fNewtonMaxIter;
	
	/** @brief Tolerance of Newton's method */
	REAL fNewtonTol;
	
	/** @brief Sets all materials in LastState */
	void SetLastState();
	
	/** @brief Sets all materials in CurrentState */
	void SetCurrentState();
	
	/** @brief Sets all materials the time step */
	void SetAllMaterialsDeltaT();
	
	static double PressaoMedia(TPZBlackOilAnalysis &an, int matid);
	static void Vazao(TPZBlackOilAnalysis &an, int matid, double & VazaoAguaSC, double  & VazaoOleoSC, double & VazaoAguaFundo, double  & VazaoOleoFundo);
};

#endif

#endif
