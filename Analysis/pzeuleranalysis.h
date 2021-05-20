/**
 * @file
 * @brief Contains TPZEulerAnalysis class which implements an analysis procedure for compressible Euler flow simulation.
 */

#ifndef PZEULERANALYSIS_H
#define PZEULERANALYSIS_H

#ifndef STATE_COMPLEX

#include "TPZLinearAnalysis.h"
#include "pzcmesh.h"
#include "pzflowcmesh.h"
#include "pzadmchunk.h"
#include "TPZMaterial.h"
#include "ConsLaw/TPZConsLaw.h"
#include "TPZMatrixSolver.h"
#include "pzdxmesh.h"
#include "pzstepsolver.h"
#include "pzblockdiag.h"

#include <iostream>

/**
 * @brief Implements an analysis procedure for computing the steady state solution of a compressible Euler flow simulation. \ref analysis "Analysis"
 * @ingroup analysis
 * @author Erick Raggio Slis dos Santos
 */
/** 
 * This class implements several tricks to obtain the steady state solution as fast as possible \n
 * It evaluates the CFL of the simulation according to the residual of the system of equations and according to the number of iterations of the Newton method at each step
 */
class TPZEulerAnalysis : public TPZLinearAnalysis
{
	
public:
	/** @brief Default constructor */
	TPZEulerAnalysis();
	/** @brief Constructor over computational mesh. out is not used */
	TPZEulerAnalysis(TPZFlowCompMesh *mesh,std::ostream &out = std::cout);
	/** @brief Simple destructor */
	~TPZEulerAnalysis();

	/** @brief Writes the computational mesh onto disk */
	void WriteCMesh( const char * str);
	
	/** @brief See declaration in the base class. */
	virtual void Run(std::ostream &out, const std::string & dxout, int dxRes);
	
	/** @brief See declaration in the base class. */
	virtual void Run(std::ostream &out) override
	{
		std::cout <<__PRETTY_FUNCTION__ << " should never be called!!!\n";
		TPZLinearAnalysis::Run(out);
	}
	
	/**
	 * @brief Sets the solution vector to be the one
	 * representing the Advanced State (Wn+1) or the advanced implicit solution.
	 */
	void SetAdvancedState();
	
	/**
	 * @brief Sets the solution vector to be the one
	 * representing the Current State (Wn) or the last explicit solution.
	 */
	void SetLastState();
	
	/**
	 * @brief Informs the Analysis class the time at which
	 * the current solution in the computational mesh belongs \n
	 * So that the materials can choose whether to contribute implicitly
	 * or explicitly.
	 */
	void SetContributionTime(TPZContributeTime time);
	
	/** @brief Indicates whether the CFL is to evolute or not */
	void SetEvolCFL(int EvolCFL);
	
	
	/**
	 * @brief Adds deltaSol to the last solution and stores it as the current solution, preparing it to contribute again.
	 * @param deltaSol [in] result of previous call to Solve.
	 * @param epsilon [out] norm of resultant fRhs.
	 */
	/** 
	 * Also updates the CompMesh soution.
	 */
	template<class TVar>
	void UpdateSolAndRhs(TPZFMatrix<TVar> & deltaSol, REAL & epsilon);
	
	/**
	 * @brief After a call to UpdateSolution, this method
	 * shifts the history in one solution,
	 * so that the current solution in the last iteration
	 * becomes the last solution in the newest
	 * iteration. \n
	 * This method does not reset the current solution. \n
	 * To reset it, UpdateSolution with a deltaSol
	 * containing any valid value (zeros, for instance)
	 * must be called.
	 */
	void UpdateHistory();
	
	/** @brief Buffers the assemblage of the rhs with respect to the last state. */
	/** Sets the current solution in the mesh when finished. */
	void BufferLastStateAssemble();
	
	/** @brief Evaluates the flux part of the residual for convergence check. */
	REAL EvaluateFluxEpsilon();
	
	/** @brief Assembles the stiffness matrix */
	/** BufferLastStateAssemble or UpdateHistory must be called first. */
	virtual void Assemble() override;
	
	/** @brief Assembles the right hand side vector. */ 
	/** BufferLastStateAssemble or UpdateHistory must be called first. */
	virtual void AssembleRhs();
	
	/** 
	 * @brief Solves an assembled stiffness matrix. 
	 * @note Not to be misunderstood as a global solver, because it performs only part of one linear iteration.
	 * @param res [out] residual of the linear system solution
	 * @param residual [in] pointer to a temporary matrix storage
	 * @param delSol [out] returns the delta Solution vector.
	 * @return 1 if everything fine, 0 otherwise (no better solution found)
	 */
	int Solve(REAL & res, TPZFMatrix<STATE> * residual, TPZFMatrix<STATE> & delSol);
	
	virtual void Solve() override
	{
		std::cout << __PRETTY_FUNCTION__ << " should never be called\n";
		TPZLinearAnalysis::Solve();
	}
	/**
	 * @brief Implements the Newton's method.
	 * @return Returns 1 if the desired tolerance was reached, 0 if it exceeded the maximal number of iterations.
	 * @param epsilon [in/out] expected tolerance in the Newton's method. Returns the value of the achieved tolerance.
	 * @param numIter [in/out] maximal number of iterations of the Newton's method. 
	 * @note Returns the number of iterations needed by the method.
	 */
	int RunNewton(REAL & epsilon, int & numIter);
	
	/** @brief Defines the time step for each material in the mesh */
	REAL ComputeTimeStep();
	
	/** @brief Settings for the Newton's method. */
	void SetNewtonCriteria(REAL epsilon, int maxIter);
	
	/** @brief Settings for the time integration method. */
	void SetTimeIntCriteria(REAL epsilon, int maxIter);
	
	/** @brief Prepares the DX graph mesh */
	TPZDXGraphMesh * PrepareDXMesh(const std::string &dxout, int dxRes);
	
	/** @brief Informs a block diagonal to be used as preconditioning */
	void SetBlockDiagonalPrecond(TPZBlockDiagonal<STATE> * blockDiag);
	
	/**
	 * @brief This method will search for the \f$ sol0 + \alpha dir\f$ solution which minimizes the residual
	 * @param residual [in/out] initial residual (of sol0) and on exiting final residual
	 * @param sol0 [in/out] solution vector
	 * @param dir [in] search direction
	 * @return 1 if a direction was found, 0 otherwise
	 */
	int LineSearch(REAL &residual, TPZFMatrix<STATE> &sol0, TPZFMatrix<STATE> &dir);
	/** @brief Trying difference between two times. */
	/** @note Actually it do nothing */
    void CompareRhs();
	
	/**
	 * @brief Evaluates the CFL control based on the newest residual norm
	 * (epsilon) and last residual norm (lastEpsilon).
	 * @param lastEpsilon [in/out] norm of the last flux vector, it is set to epsilon after CFL update.
	 * @param epsilon [in] norm of the current flux vector.
	 * @param epsilon_Newton [in] residual of the nonlinear invertion.
	 * @param timeStep [in/out] parameter to accumulate the same scales.
	 */
	/** Scales are also applied to timeStep */
	void CFLControl(REAL & lastEpsilon, REAL & epsilon, REAL & epsilon_Newton, REAL & timeStep);
	/** @brief Set to solve with GMRes algorithm */
    void SetGMResFront(REAL tol, int numiter, int numvectors);
	/** @brief Set to solve with Frontal method */
    void SetFrontalSolver();
	/** @brief Set to solve with GMRes algorithm for block diagonal matrix */
    void SetGMResBlock(REAL tol, int numiter, int numvec);
	
protected:
	
	/** @brief Stores a pointer to the computational flow mesh. */
	/** 
	 * The inherited fCompMesh is also updated with valid values to keep compatibility \n
	 * with inherited methods, but in this class fFlowCompMesh is used whenever a method exclusive \n
	 * from this class is needed.
	 */
	TPZFlowCompMesh * fFlowCompMesh;
	
	/** @brief Vector to hold the contribution of last state to the rhs. */
	TPZFMatrix<STATE> fRhsLast;
	
	/** @brief Stop criteria for the Newton loops. */
	REAL fNewtonEps;
	/** @brief Maxime iterations for iterative Newton loops */
	int fNewtonMaxIter;
	
	/** @brief Stop criteria for time integration loops */
	REAL fTimeIntEps;
	/** @brief Maxime iterations for iterative time integration loops */
	int fTimeIntMaxIter;
	
	/** @brief Indicates whether the CFL is to evolute or not. */
	int fEvolCFL;
	
	/** @brief Preconditioner */
	TPZBlockDiagonal<STATE> * fpBlockDiag;
	
	/** @brief Total number of newton iterations during this run */
	int fTotalNewton;
	
	/** @brief Indication if a frontal matrix is being used as a preconditioner */
	int fHasFrontalPreconditioner;
private:
	template<class TVar>
	void AssembleInternal();
};

#endif

#endif
