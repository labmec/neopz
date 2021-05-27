/**
 * @file
 * @brief Contains the TPZMGSolver class which represents a solution process in three steps.
 */

#ifndef TPZMGSOLVER_H
#define TPZMGSOLVER_H
#include "TPZMatrixSolver.h"
#include "pzstepsolver.h"
#include "pzfmatrix.h"
#include "pztransfer.h"


/**
 * @ingroup solver
 * @brief Represents a solution process in three steps: transfer of the residual, execute a solver on the coarse mesh, extend the solution. \ref solver "Solver"
 */
template <class TVar>
class TPZMGSolver: public TPZMatrixSolver<TVar>
{
public:
	/** @brief Default constructor */
	TPZMGSolver() : TPZRegisterClassId(&TPZMGSolver::ClassId),TPZMatrixSolver<TVar>() {}
	/** @brief Constructor of the three steps solver with transfer matrix */
	TPZMGSolver(TPZAutoPointer<TPZMatrix<TVar> > trf, const TPZMatrixSolver<TVar> &sol,
				int nvar, TPZAutoPointer<TPZMatrix<TVar> > refmat);
	/** @brief Constructor of the three steps solver */
	TPZMGSolver(TPZAutoPointer<TPZMatrix<TVar> > trf, const TPZMatrixSolver<TVar> &sol,
				int nvar);
	
	/** @brief Copy constructor */
	TPZMGSolver(const TPZMGSolver<TVar> & copy);
	/** @brief Default destructor */
	~TPZMGSolver();
	
	/** @brief Sets the transfer matrix */
	void SetTransferMatrix(TPZAutoPointer<TPZMatrix<TVar> > Refmat);
	/** @brief Clean the transfer matrix */
	void ResetTransferMatrix();
	
	/** @brief Gets the transfer matrix */
	TPZAutoPointer<TPZMatrix<TVar> > TransferMatrix()
	{
		return this->fTransfer;
	}
	
	TPZSolver * Clone() const override;
	
	void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual = 0) override;
	
	public:
int ClassId() const override;

	void Write(TPZStream &buf, int withclassid) const override;
	void Read(TPZStream &buf, void *context) override;
	
	
private:
	TPZMatrixSolver<TVar> * fCoarse{nullptr};
	int fNVar;
	/** @brief Transfer matrix */
	TPZAutoPointer<TPZMatrix<TVar> > fTransfer;
	//    TPZMatrixSolver::TPZContainer *fTransfer;
};

template <class TVar>
int TPZMGSolver<TVar>::ClassId() const{
    return Hash("TPZMGSolver") ^ TPZMatrixSolver<TVar>::ClassId() << 1;
}

#endif //TPZMGSOLVER_H
