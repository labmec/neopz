/**
 * @file
 * @brief Contains TPZNonLinearAnalysis class which implements the non linear analysis.
 */

#ifndef NONLINANALYSISH
#define NONLINANALYSISH

#include "TPZLinearAnalysis.h"
#include <iosfwd>

template<class TVar>
class TPZFMatrix;
/**
 * @brief Derived class from TPZLinearAnalysis implements non linear analysis (Newton's method). \ref analysis "Analysis"
 * @ingroup analysis
 */
class TPZNonLinearAnalysis : public TPZLinearAnalysis {
	
public:
	/** @brief Constructor with computational mesh */
	TPZNonLinearAnalysis(TPZCompMesh *mesh,std::ostream &out);
	/** @brief Default constructor */
	TPZNonLinearAnalysis();
	/** @brief Default destructor */
	virtual ~TPZNonLinearAnalysis();
	
	/** @brief It process a Newton's method to solve the non-linear problem. */
	/** It has the possibility of line search with parameter linesearch = true. */
	virtual void IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch = false, bool checkconv = false);

	/**
	 * @brief Implements a golden section line search.
	 * @param Wn
	 * @param DeltaW must be a copy
	 * @param NextW
	 * @param tol tolerance
	 * @param niter number of iterations
	 * @note Please do not put a &. It is because usually here and in derived classes fSolution was passed
	 * as DeltaW. \n But fSolution changes in the linesearch procedure when LoadSolution
	 * is called before AssembleResidual.
	 */
	REAL LineSearch(const TPZFMatrix<STATE> &Wn, TPZFMatrix<STATE> DeltaW, TPZFMatrix<STATE> &NextW, REAL tol, int niter);
	/** @brief Computes the L2 norm of the solution */
	REAL SolutionNorm();
	
	/** @note Incomplete */
	virtual void ComputeTangent(TPZFMatrix<STATE> &tangent, TPZVec<REAL> &coefs, int icase);
	/** @brief Actually return 1 */
	int NumCases();
	
	virtual void Residual(TPZFMatrix<STATE> &residual, int icase);
	
	virtual void LoadSolution() override;
	
	virtual void LoadSolution(const TPZFMatrix<STATE> &state) override;
	
	/** @brief Load solution with state as solution. But fSolution is not modified */
	void LoadState(TPZFMatrix<STATE> &state);
	
};

#endif
