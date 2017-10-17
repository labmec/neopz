/**
 * @file
 * @brief Contains the implementation of the TPZMGSolver methods. 
 */

#include "pzmgsolver.h"
#include "pztransfer.h"
#include "TPZPersistenceManager.h"

using namespace std;

template <class TVar>
TPZMGSolver<TVar>::TPZMGSolver(TPZAutoPointer<TPZTransfer<TVar> > trf, const TPZMatrixSolver<TVar> &sol, int nvar, 
							   TPZAutoPointer<TPZMatrix<TVar> > refmat) : 
TPZRegisterClassId(&TPZMGSolver::ClassId),TPZMatrixSolver<TVar>(refmat), fStep(trf) 
{
	this->fCoarse = (TPZMatrixSolver<TVar> *) sol.Clone();
	this->fNVar = nvar;
}

template <class TVar>
TPZMGSolver<TVar>::TPZMGSolver(TPZAutoPointer<TPZTransfer<TVar> > trf, const TPZMatrixSolver<TVar> &sol, int nvar) : 
TPZRegisterClassId(&TPZMGSolver::ClassId),
TPZMatrixSolver<TVar>(), fStep(trf) 
{
	this->fCoarse = (TPZMatrixSolver<TVar> *) sol.Clone();
	//  fTransfer = new TPZMatrixSolver::TPZContainer(trf);
	this->fNVar = nvar;
}

template <class TVar>
void TPZMGSolver<TVar>::Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual){
	if(!this->Matrix() || !TransferMatrix()) {
		cout << "TPZMGSolver::Solve called without a matrix pointer\n";
		DebugStop();
	}
	TPZAutoPointer<TPZMatrix<TVar> > mat = this->Matrix();
	if(result.Rows() != mat->Rows() || result.Cols() != F.Cols()) {
		result.Redim(mat->Rows(),F.Cols());
	}
	
	TPZFMatrix<TVar> FCoarse,UCoarse;
	TPZAutoPointer<TPZTransfer<TVar> > tr = TransferMatrix();
	tr->TransferResidual(F,FCoarse);
	fCoarse->Solve(FCoarse,UCoarse);
	tr->TransferSolution(UCoarse,result);
	if(residual) this->Matrix()->Residual(F,result,*residual);
}

template <class TVar>
TPZMGSolver<TVar>::TPZMGSolver(const TPZMGSolver<TVar> & copy): TPZRegisterClassId(&TPZMGSolver::ClassId),
TPZMatrixSolver<TVar>(copy), fStep(copy.fStep) {
    fCoarse = (TPZMatrixSolver<TVar> *) copy.fCoarse->Clone();
    fNVar = copy.fNVar;
}

template <class TVar>
TPZSolver<TVar> * TPZMGSolver<TVar>::Clone() const {
    return new TPZMGSolver<TVar>(*this);
}

template <class TVar>
TPZMGSolver<TVar>::~TPZMGSolver(){
    delete fCoarse;
}

template <class TVar>
void TPZMGSolver<TVar>::ResetTransferMatrix(){
	TPZAutoPointer<TPZTransfer<TVar> > reset;
	fStep = reset;
}

template <class TVar>
void TPZMGSolver<TVar>::SetTransferMatrix(TPZAutoPointer<TPZTransfer<TVar> > Refmat){
	fStep = Refmat;
}

template <class TVar>
void TPZMGSolver<TVar>::Write(TPZStream &buf, int withclassid) const
{
	TPZMatrixSolver<TVar>::Write(buf, withclassid);
        TPZPersistenceManager::WritePointer(fCoarse, &buf);
        buf.Write(&fNVar);
        TPZPersistenceManager::WritePointer(fStep.operator ->(), &buf);
}

template <class TVar>
void TPZMGSolver<TVar>::Read(TPZStream &buf, void *context)
{
	TPZMatrixSolver<TVar>::Read(buf, context);
	fCoarse = dynamic_cast<TPZMatrixSolver<TVar> *>(TPZPersistenceManager::GetInstance(&buf));
	buf.Read(&fNVar, 1);
	fStep = TPZAutoPointerDynamicCast<TPZTransfer<TVar>>(TPZPersistenceManager::GetAutoPointer(&buf));
}

template class TPZMGSolver<float>;
template class TPZMGSolver<double>;
template class TPZMGSolver<long double>;


#ifndef BORLAND
template class TPZRestoreClass<TPZMGSolver<REAL>>;
#endif
