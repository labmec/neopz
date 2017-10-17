/**
 * @file
 * @brief Contains TPZSequenceSolver class which defines sequence solvers.
 */

#ifndef TPZSEQUENCESOLVER_H
#define TPZSEQUENCESOLVER_H
#include "pzsolve.h"
#include "pzstack.h"
#include "pzfmatrix.h"

/** 
 * @ingroup solver
 * @brief Id for sequence solver
 */
#define TPZSQUENCESOLVER_ID 29281006

/**
 * @brief Defines sequence solvers. \ref solver "Solver"
 * @ingroup solver
 */
template <class TVar>
class TPZSequenceSolver : public TPZMatrixSolver<TVar> {
public:
	/**
     * @brief Constructor with initialization parameter
     * @param refmat Sets reference matrix to NILL
	 */
	TPZSequenceSolver(TPZMatrix<TVar> *refmat = 0);
	/**
     * @brief Copy constructor
     * @param copy Model object to be copied from
	 */
	TPZSequenceSolver(const TPZSequenceSolver<TVar> & copy);
	
	void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual = 0);
	
	/** @brief This method will reinitialize the solver object, including the solution procedure */  
	void ResetSolver();
	
	/**
	 * @brief This method will reset the matrix associated with the solver
	 */
	 /** This is useful when the matrix needs to be recomputed in a non linear problem */
	virtual void ResetMatrix();
	
	/** @brief Updates the values of the preconditioner based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat);
	
	void AppendSolver(TPZMatrixSolver<TVar>& solve);
	
	virtual TPZSolver<TVar> * Clone() const;
	
	/** @brief Saveable specific methods */
	
        virtual int ClassId() const;
	virtual void Write(TPZStream &buf, int withclassid) const;
	virtual void Read(TPZStream &buf, void *context);
	
	
private:    
	TPZStack < TPZMatrixSolver<TVar> * > fSolvers;
};

template <class TVar>
int TPZSequenceSolver<TVar>::ClassId() const{
    return Hash("TPZSequenceSolver") ^ TPZMatrixSolver<TVar>::ClassId() << 1;
} 

#endif //TPZSEQUENCESOLVER_H
