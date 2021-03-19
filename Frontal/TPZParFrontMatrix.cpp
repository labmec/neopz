/**
 * @file
 * @brief Contains the implementation of the TPZParFrontMatrix methods.
 */

#include "TPZParFrontMatrix.h"
#include "TPZFrontMatrix.h"
#include "pzsfulmat.h"
#include "TPZFront.h"
#include "pzstack.h"
#include "pzreal.h"
#include <math.h>

#include "tpzeqnarray.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

using namespace std;

#include "pzlog.h"

#ifdef PZ_LOG

static TPZLogger logger("pz.strmatrix.frontstructmatrix");
static TPZLogger loggerfw("pz.frontal.frontmatrix.fw");

#endif

// At the class constructor creates a thread
// this thread will be active while ParFrontMatrix is active
// It will check if a stack contains some equations to be writen to the disk

template<class TVar, class store, class front>
TPZParFrontMatrix<TVar, store, front>::TPZParFrontMatrix():
TPZRegisterClassId(&TPZParFrontMatrix::ClassId), fFinish(0), fwritelock(), fwritecond()
{
	fEqnStack.Resize(0);
}

template<class TVar, class store, class front>
TPZParFrontMatrix<TVar, store, front>::TPZParFrontMatrix(int64_t globalsize) :
TPZRegisterClassId(&TPZParFrontMatrix::ClassId), TPZFrontMatrix<TVar, store, front>(globalsize),
fFinish(0), fwritelock(), fwritecond()
{
	fEqnStack.Resize(0);
}

template<class TVar, class store, class front>
TPZParFrontMatrix<TVar, store, front>::~TPZParFrontMatrix(){
}

template<class TVar, class store, class front>
void TPZParFrontMatrix<TVar, store, front>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec < int64_t > & destinationindex)
{
	
	// message #1.3 to fFront:TPZFront
	this->fFront.AddKel(elmat, destinationindex);
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << "Frondwidth after AddKel "<< this->fFront.FrontSize();
		LOGPZ_INFO(loggerfw,sout.str())
	}
#endif

	int64_t mineq, maxeq;
	this->EquationsToDecompose(destinationindex, mineq, maxeq);
	if(maxeq >= mineq) {

		TPZEqnArray<TVar> *AuxEqn = new TPZEqnArray<TVar>;
		
		this->fFront.DecomposeEquations(mineq,maxeq,*AuxEqn);
		this->CheckCompress();
        {
            std::lock_guard<std::mutex> lock(fwritelock);
            fEqnStack.Push(AuxEqn);
            if(maxeq == this->Rows()-1){
                cout << "Decomposition finished" << endl;
                cout.flush();
                FinishWriting();
                //fStorage.ReOpen();
            }
        }
        fwritecond.notify_all();
	}
	this->fDecomposed = this->fFront.GetDecomposeType();
} 
template<class TVar, class store, class front>
void TPZParFrontMatrix<TVar, store, front>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec < int64_t > & sourceindex, TPZVec < int64_t > & destinationindex)
{
	this->fFront.AddKel(elmat, sourceindex, destinationindex);
#ifdef PZ_LOG
    if (loggerfw.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Frondwidth after AddKel "<< this->fFront.FrontSize();
		LOGPZ_DEBUG(loggerfw,sout.str())
	}
#endif
	int64_t mineq, maxeq;
	this->EquationsToDecompose(destinationindex, mineq, maxeq);

	if(maxeq >= mineq) {
		TPZEqnArray<TVar> *AuxEqn = new TPZEqnArray<TVar>;

		this->fFront.DecomposeEquations(mineq,maxeq,*AuxEqn);
		this->CheckCompress();
        std::lock_guard<std::mutex> lock(fwritelock);
		fEqnStack.Push(AuxEqn);
		if(maxeq == this->Rows()-1){
            //check if writeing is over and closes file
			cout << endl << "Decomposition finished" << endl;
			cout.flush();
			FinishWriting();
			this->fFront.Reset(0);
			//fStorage.ReOpen();
		}
        fwritecond.notify_all();
	}
	this->fDecomposed = this->fFront.GetDecomposeType();
}

template<class TVar, class store, class front>
void TPZParFrontMatrix<TVar, store, front>::FinishWriting(){
	// FinishWriting already has a lock
	cout << endl << "FinishWriting" << endl;
	cout.flush();     
	fFinish = 1;
	// FinishWriting already has a lock
    fwritecond.notify_all();
}

template<class TVar, class store, class front>
void * TPZParFrontMatrix<TVar, store, front>::WriteFile(void *t){
	TPZParFrontMatrix<TVar, store, front> *parfront = (TPZParFrontMatrix<TVar, store, front>*) t;    
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Entering WriteFile thread execution";
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	cout << endl << "Entering Decomposition" << endl;
	cout.flush();
	while(1){
		TPZStack<TPZEqnArray<TVar> *> local;
        {
            std::unique_lock<std::mutex> lock(parfront->fwritelock);
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Acquired writelock";
                LOGPZ_DEBUG(logger,sout.str())
            }
#endif
            if(parfront->fEqnStack.NElements() == 0){
                if(parfront->fFinish == 1) {
#ifdef PZ_LOG
                    if (logger.isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "Terminating WriteFile thread execution";
                        LOGPZ_DEBUG(logger,sout.str())
                    }
#endif
                    cout << "Leaving WHILE" << endl;
                    cout.flush();
                    break;
                }
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Entering cond_wait on fwritecond variable";
                    LOGPZ_DEBUG(logger,sout.str())
                }
#endif
                parfront->fwritecond.wait(lock);
            }
            
            local = parfront->fEqnStack;
            parfront->fEqnStack.Resize(0);
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Copied the equation stack releasing the writelock";
                LOGPZ_DEBUG(logger,sout.str())
            }
#endif
            
        }
		int64_t neqn = local.NElements();

		int64_t eq;
		for(eq=0; eq<neqn; eq++) {
			parfront->fStorage.AddEqnArray(local[eq]);
			delete local[eq];
		}
	}
	parfront->fStorage.FinishWriting();
	parfront->fStorage.ReOpen();
	parfront->fFinish = 0;
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Releasing writelock";
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Falling through on the write thread";
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	std::cout << "Terminating write thread\n";
	return (0);
} 

template<class TVar>
class TPZStackEqnStorage;
template<class TVar>
class TPZFileEqnStorage;
template<class TVar>
class TPZFrontSym;
template<class TVar>
class TPZFrontNonSym;


template class TPZParFrontMatrix<float, TPZStackEqnStorage<float>, TPZFrontSym<float> >;
template class TPZParFrontMatrix<float, TPZFileEqnStorage<float>, TPZFrontSym<float> >;
template class TPZParFrontMatrix<float, TPZStackEqnStorage<float>, TPZFrontNonSym<float> >;
template class TPZParFrontMatrix<float, TPZFileEqnStorage<float>, TPZFrontNonSym<float> >;

template class TPZParFrontMatrix<double, TPZStackEqnStorage<double>, TPZFrontSym<double> >;
template class TPZParFrontMatrix<double, TPZFileEqnStorage<double>, TPZFrontSym<double> >;
template class TPZParFrontMatrix<double, TPZStackEqnStorage<double>, TPZFrontNonSym<double> >;
template class TPZParFrontMatrix<double, TPZFileEqnStorage<double>, TPZFrontNonSym<double> >;

template class TPZParFrontMatrix<long double, TPZStackEqnStorage<long double>, TPZFrontSym<long double> >;
template class TPZParFrontMatrix<long double, TPZFileEqnStorage<long double>, TPZFrontSym<long double> >;
template class TPZParFrontMatrix<long double, TPZStackEqnStorage<long double>, TPZFrontNonSym<long double> >;
template class TPZParFrontMatrix<long double, TPZFileEqnStorage<long double>, TPZFrontNonSym<long double> >;

template class TPZParFrontMatrix<std::complex<float>, TPZStackEqnStorage<std::complex<float> >, TPZFrontSym<std::complex<float> > >;
template class TPZParFrontMatrix<std::complex<float>, TPZFileEqnStorage<std::complex<float> >, TPZFrontSym<std::complex<float> > >;
template class TPZParFrontMatrix<std::complex<float>, TPZStackEqnStorage<std::complex<float> >, TPZFrontNonSym<std::complex<float> > >;
template class TPZParFrontMatrix<std::complex<float>, TPZFileEqnStorage<std::complex<float> >, TPZFrontNonSym<std::complex<float> > >;

template class TPZParFrontMatrix<std::complex<double>, TPZStackEqnStorage<std::complex<double> >, TPZFrontSym<std::complex<double> > >;
template class TPZParFrontMatrix<std::complex<double>, TPZFileEqnStorage<std::complex<double> >, TPZFrontSym<std::complex<double> > >;
template class TPZParFrontMatrix<std::complex<double>, TPZStackEqnStorage<std::complex<double> >, TPZFrontNonSym<std::complex<double> > >;
template class TPZParFrontMatrix<std::complex<double>, TPZFileEqnStorage<std::complex<double> >, TPZFrontNonSym<std::complex<double> > >;

template class TPZParFrontMatrix<std::complex<long double>, TPZStackEqnStorage<std::complex<long double> >, TPZFrontSym<std::complex<long double> > >;
template class TPZParFrontMatrix<std::complex<long double>, TPZFileEqnStorage<std::complex<long double> >, TPZFrontSym<std::complex<long double> > >;
template class TPZParFrontMatrix<std::complex<long double>, TPZStackEqnStorage<std::complex<long double> >, TPZFrontNonSym<std::complex<long double> > >;
template class TPZParFrontMatrix<std::complex<long double>, TPZFileEqnStorage<std::complex<long double> >, TPZFrontNonSym<std::complex<long double> > >;

