/*
 *  tpzdohrassemblelist.h
 *  SubStruct
 *
 *  Created by Philippe Devloo on 20/07/10.
 *  Copyright 2010 UNICAMP. All rights reserved.
 *
 */
#ifndef TPZDOHRASSEMBLELIST
#define TPZDOHRASSEMBLELIST

#include "pzfmatrix.h"
#include "tpzautopointer.h"
#include "tpzdohrassembly.h"
#include "TPZSemaphore.h"

#include <pthread.h>
#include <list>
#include <semaphore.h>


struct TPZDohrAssembleItem {
	/// Initialize the assemble item based on the submesh index and size of the local contribution
	TPZDohrAssembleItem(int subindex, int size) : fSubIndex(subindex), fAssembleData(size,1,0.)
	{
	}
	/// Initialize the assemble item based on the submesh index and size of the local contribution
	TPZDohrAssembleItem(int subindex, int nrow, int ncol) : fSubIndex(subindex), fAssembleData(nrow,ncol,0.)
	{
	}
	/// substructure index
	int fSubIndex;
	/// the data which should be assembled
	TPZFMatrix fAssembleData;
};

struct TPZDohrAssembleList {
	
	/// constructor indicating the number of items that will be assembled and the target matrix
	TPZDohrAssembleList(int numitems, TPZFMatrix &output, TPZAutoPointer<TPZDohrAssembly> assembly);
	/// destructor
	~TPZDohrAssembleList();
	/// the number of items that will be assembled before returning
	int fNumItems;
	// semaphore (to wake up assembly thread)
	/*
#ifdef MACOSX
	sem_t *fSemaphore;
#else
	sem_t fSemaphore;
#endif
	 */
	TPZSemaphore fSemaphore;
	/// this is the mutex which controls the access to the list
	pthread_mutex_t fListAccessLock;
	/// this is the mutex which controls the assembly
	pthread_mutex_t fAssemblyLock;
	/// list of objects needed to be assembled
	std::list<TPZAutoPointer<TPZDohrAssembleItem> > fWork;
	/// Add an item to the list in a thread safe way
	void AddItem(TPZAutoPointer<TPZDohrAssembleItem> assembleItem);
	/// remove an item from the list
	TPZAutoPointer<TPZDohrAssembleItem> PopItem();
	/// Assembly indexes
	TPZAutoPointer<TPZDohrAssembly> fAssembleIndexes;
	/// Target Matrix
	TPZFMatrix *fOutput;
	/// Procedure which performs the assembly process
	static void *Assemble(void *voidptr);
};

#endif