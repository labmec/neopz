/**
 * @file
 * @brief Contains TPZSequenceSolver class which defines sequence solvers.
 */

#ifndef TPZSEQUENCESOLVER_H
#define TPZSEQUENCESOLVER_H
#include "TPZMatrixSolver.h"
#include "pzstack.h"
#include "pzfmatrix.h"

/**
 * @brief Defines sequence solvers. \ref solver "Solver"
 * @ingroup solver
 */
template <class TVar>
class TPZSequenceSolver : public TPZMatrixSolver<TVar> {
public:
	/**
     * @brief Constructor with initialization parameter
     * @param refmat Sets reference matrix
	 */
	TPZSequenceSolver(TPZMatrix<TVar> *refmat = nullptr);
	/**
     * @brief Copy constructor
     * @param copy Model object to be copied from
	 */
	TPZSequenceSolver(const TPZSequenceSolver<TVar> & copy);
	
	void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual = 0) override;
	
	/** @brief This method will reinitialize the solver object, including the solution procedure */  
	void ResetSolver();
	
	/**
	 * @brief This method will reset the matrix associated with the solver
	 */
	 /** This is useful when the matrix needs to be recomputed in a non linear problem */
	virtual void ResetMatrix() override;
	
	/** @brief Updates the values of the preconditioner based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZBaseMatrix> mat) override;
	
	void AppendSolver(TPZMatrixSolver<TVar>& solve);
	
	virtual TPZSolver * Clone() const override;
	
	/** @brief Saveable specific methods */
	
        int ClassId() const override;
	void Write(TPZStream &buf, int withclassid) const override;
	void Read(TPZStream &buf, void *context) override;
	
	
private:    
	TPZStack < TPZMatrixSolver<TVar> * > fSolvers;
};

template <class TVar>
int TPZSequenceSolver<TVar>::ClassId() const{
    return Hash("TPZSequenceSolver") ^ TPZMatrixSolver<TVar>::ClassId() << 1;
} 

#endif //TPZSEQUENCESOLVER_H
