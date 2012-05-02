/**
 * @file
 * @brief Contains the TPZMGSolver class which represents a solution process in three steps.
 */

#ifndef TPZMGSOLVER_H
#define TPZMGSOLVER_H
#include "pzsolve.h"
#include "pzstepsolver.h"

/** @brief Id for MG solver */
#define TPZMGSOLVER_ID 28291008

template <class TVar>
class TPZFMatrix;
template<class TVar>
class TPZTransfer;

/**
 * @ingroup solver
 * @brief Represents a solution process in three steps: transfer of the residual, execute a solver on the coarse mesh, extend the solution. \ref solver "Solver"
 */
template <class TVar>
class TPZMGSolver: public TPZMatrixSolver<TVar>
{
public:
	/** @brief Default constructor */
	TPZMGSolver() : TPZMatrixSolver<TVar>() {}
	/** @brief Constructor of the three steps solver with transfer matrix */
	TPZMGSolver(TPZAutoPointer<TPZTransfer<TVar> > trf, const TPZMatrixSolver<TVar> &sol,
				int nvar, TPZAutoPointer<TPZMatrix<TVar> > refmat);
	/** @brief Constructor of the three steps solver */
	TPZMGSolver(TPZAutoPointer<TPZTransfer<TVar> > trf, const TPZMatrixSolver<TVar> &sol,
				int nvar);
	
	/** @brief Copy constructor */
	TPZMGSolver(const TPZMGSolver<TVar> & copy);
	/** @brief Default destructor */
	~TPZMGSolver();
	
	/** @brief Sets the transfer matrix */
	void SetTransferMatrix(TPZAutoPointer<TPZTransfer<TVar
                           > > Refmat);
	/** @brief Clean the transfer matrix */
	void ResetTransferMatrix();
	
	/** @brief Gets the transfer matrix */
	TPZAutoPointer<TPZTransfer<TVar> > TransferMatrix()
	{
		return this->fStep;
	}
	
	TPZSolver<TVar> * Clone() const;
	
	void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual = 0);
	
	virtual int ClassId() const
	{
		return TPZMGSOLVER_ID;
	}
	virtual void Write(TPZStream &buf, int withclassid);
	virtual void Read(TPZStream &buf, void *context);
	
	
private:
	TPZMatrixSolver<TVar> * fCoarse;
	int fNVar;
	/** @brief Transfer matrix */
	TPZAutoPointer<TPZTransfer<TVar> > fStep;
	//    TPZMatrixSolver::TPZContainer *fTransfer;
};

#endif //TPZMGSOLVER_H
