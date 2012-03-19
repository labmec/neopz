/**
 * @file
 * @brief Contains the implementation of the TPZMGSolver methods. 
 */

#include "pzmgsolver.h"
#include "pztransfer.h"
using namespace std;

template <class TVar>
TPZMGSolver<TVar>::TPZMGSolver(TPZAutoPointer<TPZTransfer> trf, const TPZMatrixSolver<TVar> &sol, int nvar, 
                         TPZAutoPointer<TPZMatrix<TVar> > refmat) : 
                         TPZMatrixSolver<TVar>(refmat), fStep(trf) 
{
	this->fCoarse = (TPZMatrixSolver<TVar> *) sol.Clone();
	//  fTransfer = new TPZMatrixSolver::TPZContainer(trf);
	this->fNVar = nvar;
}

template <class TVar>
TPZMGSolver<TVar>::TPZMGSolver(TPZAutoPointer<TPZTransfer> trf, const TPZMatrixSolver<TVar> &sol, int nvar) : 
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
	TPZAutoPointer<TPZTransfer> tr = TransferMatrix();
	tr->TransferResidual(F,FCoarse);
	fCoarse->Solve(FCoarse,UCoarse);
	tr->TransferSolution(UCoarse,result);
	if(residual) this->Matrix()->Residual(F,result,*residual);
}

template <class TVar>
TPZMGSolver<TVar>::TPZMGSolver(const TPZMGSolver<TVar> & copy): TPZMatrixSolver<TVar>(copy), fStep(copy.fStep) {
    fCoarse = (TPZMatrixSolver<TVar> *) copy.fCoarse->Clone();
    fNVar = copy.fNVar;
    //fTransfer = copy.fTransfer;
    //    fTransfer->IncreaseRefCount();
}

template <class TVar>
TPZSolver<TVar> * TPZMGSolver<TVar>::Clone() const {
    return new TPZMGSolver<TVar>(*this);
}

template <class TVar>
TPZMGSolver<TVar>::~TPZMGSolver(){
    delete fCoarse;
    //    fTransfer->DecreaseRefCount();
	
}

template <class TVar>
void TPZMGSolver<TVar>::ResetTransferMatrix(){
	//  fTransfer->SetMatrix(0);
	TPZAutoPointer<TPZTransfer> reset;
	fStep = reset;
}

template <class TVar>
void TPZMGSolver<TVar>::SetTransferMatrix(TPZAutoPointer<TPZTransfer> Refmat){
	fStep = Refmat;
	//    if(fTransfer->Matrix() == Refmat || !fTransfer->Matrix()) {
	//        fTransfer->SetMatrix(Refmat);
	//    } else {
	//        fTransfer->DecreaseRefCount();
	//        fTransfer = new TPZContainer(Refmat);
	//    }
}

template <class TVar>
void TPZMGSolver<TVar>::Write(TPZStream &buf, int withclassid)
{
	TPZMatrixSolver<TVar>::Write(buf, withclassid);
	fCoarse->Write(buf, 1);
	buf.Write(&fNVar);
	if(fStep)
	{
		fStep->Write(buf, 1);
	}
	else
	{
		int flag = -1;
		buf.Write(&flag, 1);
	}
}

template <class TVar>
void TPZMGSolver<TVar>::Read(TPZStream &buf, void *context)
{
	TPZMatrixSolver<TVar>::Read(buf, context);
	fCoarse = dynamic_cast<TPZMatrixSolver<TVar> *>(TPZSaveable::Restore(buf, context));
	buf.Read(&fNVar, 1);
	fStep = dynamic_cast<TPZTransfer *>(TPZSaveable::Restore(buf, context));
}

template class TPZRestoreClass<TPZMGSolver<REAL>, TPZMGSOLVER_ID>;


