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

#include "pzp_thread.h"

TPZDohrAssembleList::TPZDohrAssembleList(int numitems, TPZFMatrix &output, TPZAutoPointer<TPZDohrAssembly> assembly) : fNumItems(numitems),
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
  PZP_THREAD_MUTEX_INIT(&fAssemblyLock, 0, "TPZDohrAssembleList::TPZDohrAssembleList");
  PZP_THREAD_MUTEX_INIT(&fListAccessLock, 0, "TPZDohrAssembleList::TPZDohrAssembleList");
}

TPZDohrAssembleList::~TPZDohrAssembleList()
{
  PZP_THREAD_MUTEX_DESTROY(&fAssemblyLock, "TPZDohrAssembleList::TPZDohrAssembleList");
  PZP_THREAD_MUTEX_DESTROY(&fListAccessLock, "TPZDohrAssembleList::TPZDohrAssembleList");
	/*
	 #ifdef MACOSX
	 sem_close(fSemaphore);
	 #else
	 sem_destroy(&fSemaphore);
	 #endif
	 */
}

// Add an item to the list in a thread safe way
void TPZDohrAssembleList::AddItem(TPZAutoPointer<TPZDohrAssembleItem> assembleItem)
{
  PZP_THREAD_MUTEX_LOCK(&fListAccessLock,"TPZDohrAssembleList::AddItem()");
	fWork.push_back(assembleItem);
	fSemaphore.Post();
	/*
	 #ifdef MACOSX
	 sem_post(fSemaphore);
	 #else
	 sem_post(&fSemaphore);
	 #endif
	 */
  PZP_THREAD_MUTEX_UNLOCK(&fListAccessLock,"TPZDohrAssembleList::AddItem()");
}
// remove an item from the list
TPZAutoPointer<TPZDohrAssembleItem> TPZDohrAssembleList::PopItem()
{
	TPZAutoPointer<TPZDohrAssembleItem> result;
	PZP_THREAD_MUTEX_LOCK(&fListAccessLock,"TPZDohrAssembleList::PopItem()");
	if (fWork.begin() != fWork.end()) {
		fNumItems--;
		result = *fWork.begin();
		fWork.pop_front();
	}
	PZP_THREAD_MUTEX_UNLOCK(&fListAccessLock,"TPZDohrAssembleList::PopItem()");
	return result;
}

void *TPZDohrAssembleList::Assemble(void *voidptr)
{
	TPZDohrAssembleList *myptr = (TPZDohrAssembleList *) voidptr;
	while (myptr->fNumItems > 0) {
		TPZAutoPointer<TPZDohrAssembleItem> work = myptr->PopItem();
		if (work) {
		        PZP_THREAD_MUTEX_LOCK(&myptr->fAssemblyLock,"TPZDohrAssembleList::Assemble()");
			myptr->fAssembleIndexes->Assemble(work->fSubIndex,work->fAssembleData,*(myptr->fOutput));
		        PZP_THREAD_MUTEX_UNLOCK(&myptr->fAssemblyLock,"TPZDohrAssembleList::Assemble()");
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
