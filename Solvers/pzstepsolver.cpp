/**
 * @file
 * @brief Contains the implementation of the TPZStepSolver methods.
 */

#include "pzstepsolver.h"
#include <stdlib.h>
using namespace std;

#include "pzlog.h"

#include "TPZPersistenceManager.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.converge");
#endif

template <class TVar>
TPZStepSolver<TVar>::TPZStepSolver(TPZAutoPointer<TPZMatrix<TVar> > refmat) : TPZRegisterClassId(&TPZStepSolver::ClassId),TPZMatrixSolver<TVar>(refmat), fNumIterations(-1) {
	fPrecond = 0;
	ResetSolver();
}

template <class TVar>
TPZStepSolver<TVar>::TPZStepSolver(const TPZStepSolver<TVar> & copy) : TPZRegisterClassId(&TPZStepSolver::ClassId),
TPZMatrixSolver<TVar>(copy), fNumIterations(copy.fNumIterations) , fSingular(copy.fSingular){
    fSolver = copy.fSolver;
    fDecompose = copy.fDecompose;
    fMaxIterations = copy.fMaxIterations;
    fTol = copy.fTol;
    fOverRelax = copy.fOverRelax;
    fPrecond = 0;
    if(copy.fPrecond) fPrecond = dynamic_cast<TPZMatrixSolver<TVar>*>(copy.fPrecond->Clone());
    fFromCurrent = copy.fFromCurrent;
    fNumVectors = copy.fNumVectors;//Cedric: 24/04/2003 - 12:39
}

template <class TVar>
TPZStepSolver<TVar>::~TPZStepSolver() {
	if(fPrecond) delete fPrecond;
}

/**
 * This method will reset the matrix associated with the solver \n
 * This is useful when the matrix needs to be recomputed in a non linear problem
 */
template <class TVar>
void TPZStepSolver<TVar>::ResetMatrix()
{
	TPZMatrixSolver<TVar>::ResetMatrix();
}

/**
 * @brief Decompose the system of equations if a direct solver is used
 */
template <class TVar>
void TPZStepSolver<TVar>::Decompose()
{
    if (fSolver == this->EDirect) {
        this->Matrix()->Decompose(fDecompose);
    }
}

template <class TVar>
void TPZStepSolver<TVar>::Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual){
	if(!this->Matrix()) {
		cout << "TPZMatrixSolver::Solve called without a matrix pointer\n";
		DebugStop();
	}
	
	TPZAutoPointer<TPZMatrix<TVar> > mat = this->Matrix();
    // update the matrix to which the preconditioner refers
    if(fPrecond)
    {
        
        fPrecond->UpdateFrom(this->Matrix());
    }
    
	if(result.Rows() != mat->Rows() || result.Cols() != F.Cols()) {
		result.Redim(mat->Rows(),F.Cols());
	}
	
	if(this->fScratch.Rows() != result.Rows() || this->fScratch.Cols() != result.Cols()) {
		this->fScratch.Redim(result.Rows(),result.Cols());
	}
	
	REAL tol = fTol;
	int64_t numiterations = fMaxIterations;
	switch(fSolver) {
		case TPZStepSolver::ENoSolver:
		default:
			cout << "TPZMatrixSolver::Solve called without initialized solver, Jacobi used\n";
			SetJacobi(1,0.,0);
		case TPZStepSolver::EJacobi:
			//    cout << "fScratch dimension " << fScratch.Rows() << ' ' << fScratch.Cols() << endl;
			mat->SolveJacobi(numiterations,F,result,residual,this->fScratch,tol,fFromCurrent);
            fNumIterations = numiterations;
			break;
		case TPZStepSolver::ESOR:
			mat->SolveSOR(numiterations,F,result,residual,this->fScratch,fOverRelax,tol,fFromCurrent);
            fNumIterations = numiterations;
			break;
		case TPZStepSolver::ESSOR:
			mat->SolveSSOR(numiterations,F,result,residual,this->fScratch,fOverRelax,tol,fFromCurrent);
            fNumIterations = numiterations;
			break;
		case TPZStepSolver::ECG:
			mat->SolveCG(numiterations,*fPrecond,F,result,residual,tol,fFromCurrent);
			cout << "Number of equations " << mat->Rows() << std::endl;
			cout << "Number of CG iterations " << numiterations << " tol = " << tol << endl;
            fNumIterations = numiterations;
            fTol = tol;
#ifdef PZ_LOG
            if(logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Number of equations " << mat->Rows() << std::endl;
                sout << "Number of CG iterations " << numiterations << " tol = " << tol;
                LOGPZ_DEBUG(logger,sout.str().c_str());
            }
#endif
			break;
		case TPZStepSolver::EGMRES: {
			TPZFMatrix<TVar> H(fNumVectors+1,fNumVectors+1,0.);
			mat->SolveGMRES(numiterations,*fPrecond,H,fNumVectors,F,result,residual,tol,fFromCurrent);
            fNumIterations = numiterations;
            cout << "Number of GMRES iterations " << numiterations << " tol = " << tol;
			if(numiterations == fMaxIterations || tol >= fTol)
			{
				std::cout << "GMRes tolerance was not achieved : numiter " << numiterations <<
				" tol " << tol << endl;
			}
#ifdef PZ_LOG
			{
				std::stringstream sout;
				sout << "Number of GMRES iterations " << numiterations << " tol = " << tol;
				if(logger.isDebugEnabled()) LOGPZ_DEBUG(logger,sout.str().c_str());
			}
#endif
		}
			break;
		case TPZStepSolver::EBICGSTAB: 
			mat->SolveBICGStab(numiterations, *fPrecond, F, result,residual,tol,fFromCurrent);
            fNumIterations = numiterations;
			
			if(numiterations == fMaxIterations || tol >= fTol)
			{
				std::cout << "BiCGStab tolerance was not achieved : numiter " << numiterations <<
				" tol " << tol << endl;
			}
#ifdef PZ_LOG
		{
			std::stringstream sout;
			sout << "Number of BiCGStab iterations " << numiterations << " tol = " << tol;
			LOGPZ_DEBUG(logger,sout.str().c_str());
		}
#endif
			break;
		case TPZStepSolver::EDirect:
			result = F;
			mat->SolveDirect(result,fDecompose);
			if(residual) residual->Redim(F.Rows(),F.Cols());
			break;
		case TPZStepSolver::EMultiply:
			mat->Multiply(F,result);
			if(residual) mat->Residual(result,F,*residual);
			
	}
}
template<class TVar>
void TPZStepSolver<TVar>::ResetSolver() {
	fSolver = this->ENoSolver;
	fDecompose  = ENoDecompose;
	fMaxIterations = 0;
    fNumIterations = -1;
	fTol = 0.;
	fNumVectors = 0;
	fOverRelax = 0.;
	if(fPrecond) delete fPrecond;
	fPrecond = 0;
	fFromCurrent = 0;
}
template <class TVar>
void TPZStepSolver<TVar>::SetDirect (const DecomposeType decomp){
	ResetSolver();
	fSolver = this->EDirect;
	fDecompose = decomp;
}
template <class TVar>
void TPZStepSolver<TVar>::SetCG(const int64_t numiterations, const TPZMatrixSolver<TVar> &pre, const REAL tol, const int64_t FromCurrent){
	ResetSolver();
	fSolver = this->ECG;
	fMaxIterations = numiterations;
    fNumIterations = -1;
	fTol = tol;
	//	fPrecond = &pre;
	if(fPrecond) delete fPrecond;
	fPrecond = dynamic_cast<TPZMatrixSolver<TVar>*>(pre.Clone());
	fFromCurrent = FromCurrent;
}
template<class TVar>
void TPZStepSolver<TVar>::SetGMRES(const int64_t numiterations, const int numvectors, const TPZMatrixSolver<TVar> &pre, const REAL tol, const int64_t FromCurrent){
	ResetSolver();
	fSolver = this->EGMRES;
	fNumVectors = numvectors;
	fMaxIterations = numiterations;
    fNumIterations = -1;
	fTol = tol;
	//	fPrecond = &pre;
	if(fPrecond) delete fPrecond;
	fPrecond = dynamic_cast<TPZMatrixSolver<TVar>*>(pre.Clone());
	fFromCurrent = FromCurrent;
}
template<class TVar>
void TPZStepSolver<TVar>::SetBiCGStab(const int64_t numiterations, const TPZMatrixSolver<TVar>&pre,const REAL tol,const int64_t FromCurrent){
	ResetSolver();
	fSolver = this->EBICGSTAB;
	fMaxIterations = numiterations;
    fNumIterations = -1;
	fTol = tol;
	//	fPrecond = &pre;
	if(fPrecond) delete fPrecond;
	fPrecond = dynamic_cast<TPZMatrixSolver<TVar>*>(pre.Clone());
	fFromCurrent = FromCurrent;
}
template<class TVar>
void TPZStepSolver<TVar>::SetJacobi(const int64_t numiterations, const REAL tol, const int64_t FromCurrent) {
	ResetSolver();
	fSolver = this->EJacobi;
	fMaxIterations = numiterations;
    fNumIterations = -1;
	fTol = tol;
	fFromCurrent = FromCurrent;
}
template <class TVar>void TPZStepSolver<TVar>::SetSSOR(const int64_t numiterations,const REAL overrelax,const REAL tol,const int64_t FromCurrent) {
	ResetSolver();
	fSolver = this->ESSOR;
	fOverRelax = overrelax;
	fMaxIterations = numiterations;
    fNumIterations = -1;
	fTol = tol;
	fFromCurrent = FromCurrent;
}
template <class TVar>
void TPZStepSolver<TVar>::SetSOR(const int64_t numiterations,const REAL overrelax,const REAL tol,const int64_t FromCurrent){
	ResetSolver();
	fSolver = this->ESOR;
	fMaxIterations = numiterations;
    fNumIterations = -1;
	fOverRelax = overrelax;
	fTol = tol;
	fFromCurrent = FromCurrent;
}
template<class TVar>
void TPZStepSolver<TVar>::SetMultiply() {
	ResetSolver();
	fSolver = this->EMultiply;
}


/*!
 \fn TPZStepSolver::SetPreconditioner(TPZSolver &solve);
 */
template <class TVar>
void TPZStepSolver<TVar>::SetPreconditioner(TPZMatrixSolver<TVar> &solve)
{
    if (fSolver == this->EDirect) {
        DebugStop();
    }
	if(fPrecond) delete fPrecond;
	fPrecond = dynamic_cast<TPZMatrixSolver<TVar>*>(solve.Clone());
}

template <class TVar>
void TPZStepSolver<TVar>::Write(TPZStream &buf, int withclassid) const {
    TPZMatrixSolver<TVar>::Write(buf, withclassid);
    TPZPersistenceManager::WritePointer(fPrecond, &buf);
    int lfSolver = fSolver;
    buf.Write(&lfSolver, 1);
    int lfDT = fDecompose;
    buf.Write(&lfDT, 1);
    buf.Write(&fMaxIterations, 1);
    buf.Write(&fNumVectors, 1);
    buf.Write(&fTol, 1);
    buf.Write(&fOverRelax, 1);
    buf.Write(&fFromCurrent, 1);
    int64_t size = fSingular.size();
    buf.Write(&size, 1);
    std::list<int64_t>::const_iterator it = fSingular.begin();
    for (; it != fSingular.end(); it++) {
        buf.Write(&*it, 1);
    }
}

template <class TVar>
void TPZStepSolver<TVar>::Read(TPZStream &buf, void *context)
{
	TPZMatrixSolver<TVar>::Read(buf, context);
	fPrecond = dynamic_cast<TPZMatrixSolver<TVar> *>(TPZPersistenceManager::GetInstance(&buf));
	
	int lfSolver = 0;
	buf.Read(&lfSolver, 1);
	fSolver = (typename TPZMatrixSolver<TVar>::MSolver)lfSolver;
	int lfDT = 0;
	buf.Read(&lfDT, 1);
	fDecompose = (DecomposeType)lfDT;
	buf.Read(&fMaxIterations, 1);
	buf.Read(&fNumVectors, 1);
	buf.Read(&fTol, 1);
	buf.Read(&fOverRelax, 1);
	buf.Read(&fFromCurrent, 1);
	int64_t size = 0;
	buf.Read(&size, 1);
	fSingular.resize(size);
	std::list<int64_t>::iterator it = fSingular.begin();
	for(;it != fSingular.end(); it++)
	{
		buf.Read(&*it, 1);
	}
}

template class TPZStepSolver<float>;
template class TPZStepSolver<double>;
template class TPZStepSolver<long double>;

template class TPZStepSolver<std::complex<float> >;
template class TPZStepSolver<std::complex<double> >;
template class TPZStepSolver<std::complex<long double> >;

#ifndef BORLAND
template class TPZRestoreClass< TPZStepSolver<float>>;
template class TPZRestoreClass< TPZStepSolver<double>>;
template class TPZRestoreClass< TPZStepSolver<long double>>;

template class TPZRestoreClass< TPZStepSolver<std::complex<float>>>;
template class TPZRestoreClass< TPZStepSolver<std::complex<double>>>;
template class TPZRestoreClass< TPZStepSolver<std::complex<long double>>>;
#endif
