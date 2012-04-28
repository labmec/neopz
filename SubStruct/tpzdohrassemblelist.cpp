/**
 * @file
 * @brief Contains the implementation of the TPZDohrAssembleList methods. 
 */
/*
 *  tpzdohrassemblelist.cpp
 *  SubStruct
 *
 *  Created by Philippe Devloo on 20/07/10.
 *  Copyright 2010 UNICAMP. All rights reserved.
 *
 */

#include "tpzdohrassemblelist.h"

template<class TVar>
TPZDohrAssembleList<TVar>::TPZDohrAssembleList(int numitems, TPZFMatrix<TVar> &output, TPZAutoPointer<TPZDohrAssembly<TVar> > assembly) : fNumItems(numitems),
fAssembleIndexes(assembly), fOutput(&output)
{
	/*
	 #ifdef MACOSX
	 std::stringstream sout;
	 static int counter = 0;
	 sout << "DohrAssemblySem" << counter++;
	 fSemaphore = sem_open(sout.str().c_str(), O_CREAT,777,1);
	 if(fSemaphore == SEM_FAILED)
	 {
	 std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
	 }
	 #else
	 int sem_result = sem_init(&fSemaphore,0,0);
	 if(sem_result != 0)
	 {
	 std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
	 }
	 #endif
	 */	
	pthread_mutex_init(&fAssemblyLock, 0);
	pthread_mutex_init(&fListAccessLock, 0);
}

template<class TVar>
TPZDohrAssembleList<TVar>::~TPZDohrAssembleList()
{
	pthread_mutex_destroy(&fAssemblyLock);
	pthread_mutex_destroy(&fListAccessLock);
	/*
	 #ifdef MACOSX
	 sem_close(fSemaphore);
	 #else
	 sem_destroy(&fSemaphore);
	 #endif
	 */
}

// Add an item to the list in a thread safe way
template<class TVar>
void TPZDohrAssembleList<TVar>::AddItem(TPZAutoPointer<TPZDohrAssembleItem<TVar> > assembleItem)
{
	pthread_mutex_lock(&fListAccessLock);
	fWork.push_back(assembleItem);
	fSemaphore.Post();
	/*
	 #ifdef MACOSX
	 sem_post(fSemaphore);
	 #else
	 sem_post(&fSemaphore);
	 #endif
	 */
	pthread_mutex_unlock(&fListAccessLock);
}
// remove an item from the list
template<class TVar>
TPZAutoPointer<TPZDohrAssembleItem<TVar> > TPZDohrAssembleList<TVar>::PopItem()
{
	TPZAutoPointer<TPZDohrAssembleItem<TVar> > result;
	pthread_mutex_lock(&fListAccessLock);
	if (fWork.begin() != fWork.end()) {
		fNumItems--;
		result = *fWork.begin();
		fWork.pop_front();
	}
	pthread_mutex_unlock(&fListAccessLock);
	return result;
}

template<class TVar>
void *TPZDohrAssembleList<TVar>::Assemble(void *voidptr)
{
	TPZDohrAssembleList *myptr = (TPZDohrAssembleList *) voidptr;
	while (myptr->fNumItems > 0) {
		TPZAutoPointer<TPZDohrAssembleItem<TVar> > work = myptr->PopItem();
		if (work) {
			pthread_mutex_lock(&myptr->fAssemblyLock);
			myptr->fAssembleIndexes->Assemble(work->fSubIndex,work->fAssembleData,*(myptr->fOutput));
			pthread_mutex_unlock(&myptr->fAssemblyLock);
		}
		else {
			// wait for a signal
			myptr->fSemaphore.Wait();
			/*
			 #ifdef MACOSX
			 sem_wait(myptr->fSemaphore);
			 #else
			 sem_wait(&myptr->fSemaphore);
			 #endif
			 */
		}
		
	}
	return voidptr;
}

template struct TPZDohrAssembleList<double>;
template struct TPZDohrAssembleList<std::complex<double> >;