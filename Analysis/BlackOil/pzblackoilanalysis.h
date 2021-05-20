/**
 * @file
 * @brief Contains TPZBlackOilAnalysis class derived from TPZNonLinearAnalysis class.
 */

#ifndef BLACKOILANALH
#define BLACKOILANALH

#include "pznonlinanalysis.h"
#include <iosfwd>

class TPZCompMesh;
class TPZFStructMatrix;


/**
 * @brief Derived from TPZNonLinearAnalysis class. \ref analysis "Analysis"
 * @ingroup analysis
 */
class TPZBlackOilAnalysis : public TPZNonLinearAnalysis {
	
private:
	/** @brief To store last load vector */
	TPZFMatrix<STATE> fLastState;
	
public:
	/** @brief Constructor for given time step */
	TPZBlackOilAnalysis(TPZCompMesh *mesh, double TimeStep, std::ostream &out = std::cout);
	/** @brief Simple destructor */
	~TPZBlackOilAnalysis();
	
	/** @brief Assemble residual vector and tangent matrix */
	virtual void Assemble() override;
	
	/** @brief Assemble only the residual vector */
	virtual void AssembleResidual() override;
	
	/** @brief Invert the algebraic system */
	virtual void Solve() override;
	
//	virtual void Run(std::ostream &out = std::cout, bool linesearch = true);
    virtual void Run(std::ostream &out,bool linesearch);
    virtual void Run(std::ostream &out = std::cout) override {
        TPZBlackOilAnalysis::Run(out,true);
    }
	
	virtual void PostProcess(int resolution) override { TPZLinearAnalysis::PostProcess(resolution);}
	
	virtual void PostProcess(int resolution, int dimension) override;
	
	virtual void PostProcess(TPZVec<REAL> &loc, std::ostream &out = std::cout) override;
	
	/** @brief Defines max number of steps and steady state convergence tolerance. */
	void SetConvergence(int niter, REAL eps, bool ForceAllSteps = true);
	
	/** @brief Defines properties of DX file */
	void SetSaveFrequency(int SaveFrequency, int resolution);
	
	/** @brief Defines max number of steps and error convergence tolerance for Newton's method. */
	void SetNewtonConvergence(int niter, REAL eps);
	/** @brief Gets time step */
	REAL &TimeStep();
	
	/** @brief Sets a vector InitialSol as initial solution to start iterative process */ 
	void SetInitialSolution(TPZFMatrix<STATE> & InitialSol);
	/** @brief Zeroes solution to start iterative process */
	void SetInitialSolutionAsZero();
	
protected:
	
	/** @brief Simulation time step */
	REAL fTimeStep;
	
	/** @brief Stores the current time of the simulation, to compare with \f$ TotalTime = fTimeSet * fNIter\f$ */
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
	
	/** @brief Computes the average pressure */
	static double AveragePressure(TPZBlackOilAnalysis &an, int matid);
	/** @brief Computes the oil flow */
	static void Flow(TPZBlackOilAnalysis &an, int matid, double & WaterFlowSC, double  & OilFlowSC, double & WaterFlowBottom, double  & OilFlowBottom);
};

#endif

