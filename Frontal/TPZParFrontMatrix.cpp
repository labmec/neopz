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
#include <pthread.h>

#include "tpzeqnarray.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

using namespace std;

#include "pzlog.h"

#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("pz.strmatrix.frontstructmatrix"));
static LoggerPtr loggerfw(Logger::getLogger("pz.frontal.frontmatrix.fw"));

#endif

/** @brief Initializing semaphore */
pthread_mutex_t mutex_write = PTHREAD_MUTEX_INITIALIZER;
/** @brief Initializing condition */
pthread_cond_t conda_write = PTHREAD_COND_INITIALIZER;

// At the class constructor creates a thread
// this thread will be active while ParFrontMatrix is active
// It will check if a stack contains some equations to be writen to the disk

template<class store, class front>
TPZParFrontMatrix<store, front>::TPZParFrontMatrix():
fFinish(0)
{
	fEqnStack.Resize(0);
	pthread_mutex_t mlocal = PTHREAD_MUTEX_INITIALIZER;
	fwritelock = mlocal;
	pthread_cond_t clocal = PTHREAD_COND_INITIALIZER;
	fwritecond = clocal;
	/*	fFront.Reset();
	 fStorage.Reset();
	 fNumElConnected.Resize(0);
	 fLastDecomposed = -1;
	 fNumEq=0;
	 */
}

template<class store, class front>
TPZParFrontMatrix<store, front>::TPZParFrontMatrix(int globalsize) :
TPZFrontMatrix<store, front>(globalsize),
fFinish(0)
{
	fEqnStack.Resize(0);
	pthread_mutex_t mlocal = PTHREAD_MUTEX_INITIALIZER;
	fwritelock = mlocal;
	pthread_cond_t clocal = PTHREAD_COND_INITIALIZER;
	fwritecond = clocal;
	/*	fFront.Reset(globalsize);
	 fStorage.Reset();
	 fNumElConnected.Resize(0);
	 fLastDecomposed = -1;
	 fNumEq=globalsize;*/
}
template<class store, class front>
TPZParFrontMatrix<store, front>::~TPZParFrontMatrix(){
}

template<class store, class front>
void TPZParFrontMatrix<store, front>::AddKel(TPZFMatrix & elmat, TPZVec < int > & destinationindex)
{
	
	// message #1.3 to fFront:TPZFront
	this->fFront.AddKel(elmat, destinationindex);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Frondwidth after AddKel "<< this->fFront.FrontSize();
		LOGPZ_INFO(loggerfw,sout.str())
	}
#endif
    /*      cout << "destination index" << endl;
	 int i;
	 for(i=0;i<destinationindex.NElements();i++) cout << destinationindex[i] << " ";
	 cout << endl;
	 cout.flush();
	 elmat.Print("Element Matrix");
	 */
	int mineq, maxeq;
	this->EquationsToDecompose(destinationindex, mineq, maxeq);
	TPZEqnArray *AuxEqn = new TPZEqnArray;
	if(maxeq >= mineq) {
		//               if(!(maxeq%10)){
		//	               cout << (100*maxeq/fNumEq) << " % Decomposed" << endl;
		//	               cout.flush();
		//              }
		
		this->fFront.DecomposeEquations(mineq,maxeq,*AuxEqn);
		this->CheckCompress();
		pthread_mutex_lock(&fwritelock);
		fEqnStack.Push(AuxEqn);
		if(maxeq == this->Rows()-1){
			cout << "Decomposition finished" << endl;
			cout.flush();
			FinishWriting();
			//fStorage.ReOpen();
		}
		pthread_mutex_unlock(&fwritelock);
		pthread_cond_signal(&fwritecond);
	}
	this->fDecomposed = this->fFront.GetDecomposeType();
} 
template<class store, class front>
void TPZParFrontMatrix<store, front>::AddKel(TPZFMatrix & elmat, TPZVec < int > & sourceindex, TPZVec < int > & destinationindex)
{
	this->fFront.AddKel(elmat, sourceindex, destinationindex);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Frondwidth after AddKel "<< this->fFront.FrontSize();
		LOGPZ_INFO(loggerfw,sout.str())
	}
#endif
	//	EquationsToDecompose(destinationindex);
	//         cout << "AddKel::destination index 2" << endl;
	//          for(i=0;i<destinationindex.NElements();i++) cout << destinationindex[i] << " ";
	//         cout << endl;
	//         cout.flush();
	//          elmat.Print("AddKel: Element Matrix 2");
	int mineq, maxeq;
	this->EquationsToDecompose(destinationindex, mineq, maxeq);
	TPZEqnArray *AuxEqn = new TPZEqnArray;
	if(maxeq >= mineq) {
		//	     if(!(maxeq%10)){
		//	          cout << (100*maxeq/fNumEq) << " % Decomposed" << endl;
		//	          cout.flush();
		//          }
		
		this->fFront.DecomposeEquations(mineq,maxeq,*AuxEqn);
		this->CheckCompress();
		//fStorage.AddEqnArray(&AuxEqn);
		//adds an equation to a stack!!!
		//some sort of lock here
		//          fEqnStack->Push(&AuxEqn);
		pthread_mutex_lock(&fwritelock);
		fEqnStack.Push(AuxEqn);
		if(maxeq == this->Rows()-1){
            //check if writeing is over and closes file
			cout << endl << "Decomposition finished" << endl;
			cout.flush();
			FinishWriting();
			this->fFront.Reset(0);
			//fStorage.ReOpen();
		}
		pthread_mutex_unlock(&fwritelock);
     	pthread_cond_signal(&fwritecond);
	}
	this->fDecomposed = this->fFront.GetDecomposeType();
}

template<class store, class front>
void TPZParFrontMatrix<store, front>::FinishWriting(){
	// FinishWriting already has a lock
	//pthread_mutex_lock(&fwritelock);
	cout << endl << "FinishWriting" << endl;
	cout.flush();     
	fFinish = 1;
	// FinishWriting already has a lock
	//pthread_mutex_unlock(&fwritelock);
	pthread_cond_signal(&fwritecond);
} 

template<class store, class front>
void * TPZParFrontMatrix<store, front>::WriteFile(void *t){
	TPZParFrontMatrix<store, front> *parfront = (TPZParFrontMatrix<store, front>*) t;    
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Entering WriteFile thread execution";
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	cout << endl << "Entering Decomposition" << endl;
	cout.flush();
	while(1){
		TPZStack<TPZEqnArray *> local;
		pthread_mutex_lock(&parfront->fwritelock);
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Acquired writelock";
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		if(parfront->fEqnStack.NElements() == 0){
			if(parfront->fFinish == 1) {
#ifdef LOG4CXX
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
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "Entering cond_wait on fwritecond variable";
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			pthread_cond_wait(&parfront->fwritecond, &parfront->fwritelock);
		}
		
		local = parfront->fEqnStack;
		parfront->fEqnStack.Resize(0);
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Copied the equation stack releasing the writelock";
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		
		pthread_mutex_unlock(&parfront->fwritelock);
		int neqn = local.NElements();
		
		/*          nlocal++;
		 if(!(nlocal%200)) cout << endl << "         Decomposing  " << neqn << " " << nlocal << " on thread " << pthread_self() << endl;
		 if(!(nlocal%20)) cout << nlocal << endl;
		 cout << '#';
		 cout.flush();
		 */          
		int eq;
		for(eq=0; eq<neqn; eq++) {
			parfront->fStorage.AddEqnArray(local[eq]);
			delete local[eq];
		}
		
		
	}
	parfront->fStorage.FinishWriting();
	parfront->fStorage.ReOpen();
	parfront->fFinish = 0;
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Releasing writelock";
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	pthread_mutex_unlock(&parfront->fwritelock);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Falling through on the write thread";
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	std::cout << "Terminating write thread\n";
	return (0);
} 

class TPZStackEqnStorage;
class TPZFileEqnStorage;
class TPZFrontSym;
class TPZFrontNonSym;

template class TPZParFrontMatrix<TPZStackEqnStorage, TPZFrontSym>;
template class TPZParFrontMatrix<TPZFileEqnStorage, TPZFrontSym>;
template class TPZParFrontMatrix<TPZStackEqnStorage, TPZFrontNonSym>;
template class TPZParFrontMatrix<TPZFileEqnStorage, TPZFrontNonSym>;



