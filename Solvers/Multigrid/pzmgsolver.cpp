/**
 * @file
 * @brief Contains the implementation of the TPZMGSolver methods. 
 */

#include "pzmgsolver.h"
#include "pztransfer.h"
#include "TPZPersistenceManager.h"

using namespace std;

template <class TVar>
TPZMGSolver<TVar>::TPZMGSolver(TPZAutoPointer<TPZMatrix<TVar> > trf, const TPZMatrixSolver<TVar> &sol, int nvar,
							   TPZAutoPointer<TPZMatrix<TVar> > refmat) : 
TPZRegisterClassId(&TPZMGSolver::ClassId),TPZMatrixSolver<TVar>(refmat), fTransfer(trf)
{
	this->fCoarse = (TPZMatrixSolver<TVar> *) sol.Clone();
	this->fNVar = nvar;
}

template <class TVar>
TPZMGSolver<TVar>::TPZMGSolver(TPZAutoPointer<TPZMatrix<TVar> > trf, const TPZMatrixSolver<TVar> &sol, int nvar) :
TPZRegisterClassId(&TPZMGSolver::ClassId),
TPZMatrixSolver<TVar>(), fTransfer(trf)
{
	this->fCoarse = (TPZMatrixSolver<TVar> *) sol.Clone();
	//  fTransfer = new TPZMatrixSolver::TPZContainer(trf);
	this->fNVar = nvar;
}

template <class TVar>
void TPZMGSolver<TVar>::Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual){
	if((!this->Matrix() && residual != 0) || !TransferMatrix()) {
		cout << "TPZMGSolver::Solve called without a matrix pointer\n";
		DebugStop();
	}
    TPZAutoPointer<TPZMatrix<TVar> > tr = TransferMatrix();
	if(result.Rows() != tr->Cols() || result.Cols() != F.Cols()) {
		result.Redim(tr->Rows(),F.Cols());
	}
	
	TPZFMatrix<TVar> FCoarse,UCoarse;
	tr->Multiply(F,FCoarse,0);
    double norm = Norm(FCoarse);
    if(norm != 0.)
    {
        fCoarse->Solve(FCoarse,UCoarse);
    }
	tr->Multiply(UCoarse,result,1);
	if(residual) this->Matrix()->Residual(F,result,*residual);
}

template <class TVar>
TPZMGSolver<TVar>::TPZMGSolver(const TPZMGSolver<TVar> & copy): TPZRegisterClassId(&TPZMGSolver::ClassId),
TPZMatrixSolver<TVar>(copy), fTransfer(copy.fTransfer) {
    fCoarse = (TPZMatrixSolver<TVar> *) copy.fCoarse->Clone();
    fNVar = copy.fNVar;
}

template <class TVar>
TPZSolver * TPZMGSolver<TVar>::Clone() const {
    return new TPZMGSolver<TVar>(*this);
}

template <class TVar>
TPZMGSolver<TVar>::~TPZMGSolver(){
  if(fCoarse) delete fCoarse;
}

template <class TVar>
void TPZMGSolver<TVar>::ResetTransferMatrix(){
	TPZAutoPointer<TPZMatrix<TVar> > reset;
	fTransfer = reset;
}

template <class TVar>
void TPZMGSolver<TVar>::SetTransferMatrix(TPZAutoPointer<TPZMatrix<TVar> > Refmat){
	fTransfer = Refmat;
}

template <class TVar>
void TPZMGSolver<TVar>::Write(TPZStream &buf, int withclassid) const
{
	TPZMatrixSolver<TVar>::Write(buf, withclassid);
        TPZPersistenceManager::WritePointer(fCoarse, &buf);
        buf.Write(&fNVar);
        TPZPersistenceManager::WritePointer(fTransfer.operator ->(), &buf);
}

template <class TVar>
void TPZMGSolver<TVar>::Read(TPZStream &buf, void *context)
{
	TPZMatrixSolver<TVar>::Read(buf, context);
	fCoarse = dynamic_cast<TPZMatrixSolver<TVar> *>(TPZPersistenceManager::GetInstance(&buf));
	buf.Read(&fNVar, 1);
	fTransfer = TPZAutoPointerDynamicCast<TPZMatrix<TVar>>(TPZPersistenceManager::GetAutoPointer(&buf));
}

template class TPZMGSolver<float>;
template class TPZMGSolver<double>;
template class TPZMGSolver<long double>;


#ifndef BORLAND
template class TPZRestoreClass<TPZMGSolver<REAL>>;
#endif
