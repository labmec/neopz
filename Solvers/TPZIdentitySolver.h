/**
 * @file
 * @brief Contains TPZIdentitySolver class which acts as an identity matrix preconditioner.
 */

#ifndef TPZIDENTITYSOLVER_H
#define TPZIDENTITYSOLVER_H
#include "TPZMatrixSolver.h"
#include "pzfmatrix.h"


/**
 * @ingroup solver
 * @brief This class behaves as an identity matrix preconditioner. \ref solver "Solver"
 */
template<class TVar>
class TPZIdentitySolver: public TPZMatrixSolver<TVar>
{
public:

	//!Default constructor
	TPZIdentitySolver() = default;
	/**
	 * @brief Solves the system of linear equations stored in current matrix \n
	 * As this class implements only a copy operation, it just copies u to F;
	 * @param F contains Force vector
	 * @param result contains the solution
	 * @param residual [out] residual computed
	 */
	inline void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual) override
	{
		result = F;
	}
	
	/** @brief Clones the current object returning a pointer of type TPZSolver */
	TPZSolver *Clone() const override
	{
		return new TPZIdentitySolver(*this);
	}
};

#endif //TPZIDENTITYSOLVER_H
