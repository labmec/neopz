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

class TPZFMatrix;
class TPZTransfer;

/**
 * @ingroup solver
 * @brief Represents a solution process in three steps: transfer of the residual, execute a solver on the coarse mesh, extend the solution. \ref solver "Solver"
 */
class TPZMGSolver: public TPZMatrixSolver
{
public:
	TPZMGSolver() : TPZMatrixSolver() {}
	TPZMGSolver(TPZAutoPointer<TPZTransfer> trf, const TPZMatrixSolver &sol,
				int nvar, TPZAutoPointer<TPZMatrix> refmat);
	TPZMGSolver(TPZAutoPointer<TPZTransfer> trf, const TPZMatrixSolver &sol,
				int nvar);
	
	TPZMGSolver(const TPZMGSolver & copy);
	
	~TPZMGSolver();
	
	void SetTransferMatrix(TPZAutoPointer<TPZTransfer> Refmat);
	
	void ResetTransferMatrix();
	
	TPZAutoPointer<TPZTransfer> TransferMatrix()
	{
		return fStep;
	}
	
	TPZSolver * Clone() const;
	
	void Solve(const TPZFMatrix &F, TPZFMatrix &result, TPZFMatrix *residual = 0);
	
	virtual int ClassId() const
	{
		return TPZMGSOLVER_ID;
	}
	virtual void Write(TPZStream &buf, int withclassid);
	virtual void Read(TPZStream &buf, void *context);
	
	
private:
	TPZMatrixSolver * fCoarse;
	int fNVar;
	TPZAutoPointer<TPZTransfer> fStep;
	//    TPZMatrixSolver::TPZContainer *fTransfer;
};

#endif //TPZMGSOLVER_H
