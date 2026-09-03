/**
 * @file
 * @brief Contains TPZJacobiPrecond class which implements a Jacobi preconditioner.
 */

#ifndef TPZJACOBIPRECOND_H
#define TPZJACOBIPRECOND_H
#include "TPZMatrixSolver.h"
#include "pzfmatrix.h"


/**
 * @ingroup solver
 * @brief This class behaves as a Jacobi matrix preconditioner. \ref solver "Solver"
 * @note fContainer will only store inverted diagonal entries of reference matrix
 */
template<class TVar>
class TPZJacobiPrecond: public TPZMatrixSolver<TVar>
{
public:

	//!Default constructor
	TPZJacobiPrecond() = default;

  //!Initializes reference matrix and diagonal vector
  TPZJacobiPrecond(TPZAutoPointer<TPZMatrix<TVar> >  refmat);

  //!Updates diagonal vector if matrix is the same as reference matrix
  void UpdateFrom(TPZAutoPointer<TPZBaseMatrix> mat) override;
	/**
	 * @brief Solves the system of linear equations stored in current matrix \n
	 * As this class implements only a copy operation, it just copies u to F;
	 * @param F contains Force vector
	 * @param result contains the solution
	 * @param residual [out] residual computed
	 */
	void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual) override;
	
	/** @brief Clones the current object returning a pointer of type TPZSolver */
	TPZSolver *Clone() const override
	{
		return new TPZJacobiPrecond(*this);
	}
  /** @brief Sets coloring of blocks (for GS-like iteration)*/
  void SetColoring(const TPZVec<int64_t> &colors, const int numcolors);
protected:
  //! Vector storing coloring info
  TPZVec<TPZVec<int64_t>> fColorVec;
};

#endif //TPZJACOBIPRECOND_H