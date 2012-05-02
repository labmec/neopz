/**
 * @file
 * @brief Contains the implementation of the TPZDohrAssembleList methods. 
 */

#include "tpzdohrassemblelist.h"

template<class TVar>
TPZDohrAssembleList<TVar>::TPZDohrAssembleList(int numitems, TPZFMatrix<TVar> &output, TPZAutoPointer<TPZDohrAssembly<TVar> > assembly) : fNumItems(numitems),
fAssembleIndexes(assembly), fOutput(&output)
{
	pthread_mutex_init(&fAssemblyLock, 0);
	pthread_mutex_init(&fListAccessLock, 0);
}

template<class TVar>
TPZDohrAssembleList<TVar>::~TPZDohrAssembleList()
{
	pthread_mutex_destroy(&fAssemblyLock);
	pthread_mutex_destroy(&fListAccessLock);
}

// Add an item to the list in a thread safe way
template<class TVar>
void TPZDohrAssembleList<TVar>::AddItem(TPZAutoPointer<TPZDohrAssembleItem<TVar> > assembleItem)
{
	pthread_mutex_lock(&fListAccessLock);
	fWork.push_back(assembleItem);
	fSemaphore.Post();
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
		}
	}
	return voidptr;
}

template struct TPZDohrAssembleList<float>;
template struct TPZDohrAssembleList<double>;
template struct TPZDohrAssembleList<long double>;

template struct TPZDohrAssembleList<std::complex<float> >;
template struct TPZDohrAssembleList<std::complex<double> >;
template struct TPZDohrAssembleList<std::complex<long double> >;