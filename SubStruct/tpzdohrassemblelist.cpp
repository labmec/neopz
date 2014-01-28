/**
 * @file
 * @brief Contains the implementation of the TPZDohrAssembleList methods.
 */

#include "tpzdohrassemblelist.h"

#include "pz_pthread.h"

template<class TVar>
TPZDohrAssembleList<TVar>::TPZDohrAssembleList(int numitems, TPZFMatrix<TVar> &output, TPZAutoPointer<TPZDohrAssembly<TVar> > assembly) : fNumItems(numitems),
fAssembleIndexes(assembly), fOutput(&output)
{
    PZ_PTHREAD_MUTEX_INIT(&fAssemblyLock, 0, "TPZDohrAssembleList::TPZDohrAssembleList");
    PZ_PTHREAD_MUTEX_INIT(&fListAccessLock, 0, "TPZDohrAssembleList::TPZDohrAssembleList");
}

template<class TVar>
TPZDohrAssembleList<TVar>::~TPZDohrAssembleList()
{
    PZ_PTHREAD_MUTEX_DESTROY(&fAssemblyLock, "TPZDohrAssembleList::TPZDohrAssembleList");
    PZ_PTHREAD_MUTEX_DESTROY(&fListAccessLock, "TPZDohrAssembleList::TPZDohrAssembleList");
}

// Add an item to the list in a thread safe way
template<class TVar>
void TPZDohrAssembleList<TVar>::AddItem(TPZAutoPointer<TPZDohrAssembleItem<TVar> > assembleItem)
{
    PZ_PTHREAD_MUTEX_LOCK(&fListAccessLock,"TPZDohrAssembleList::AddItem()");
	fWork.push_back(assembleItem);
	fSemaphore.Post();
    PZ_PTHREAD_MUTEX_UNLOCK(&fListAccessLock,"TPZDohrAssembleList::AddItem()");
}
// remove an item from the list
template<class TVar>
TPZAutoPointer<TPZDohrAssembleItem<TVar> > TPZDohrAssembleList<TVar>::PopItem()
{
	TPZAutoPointer<TPZDohrAssembleItem<TVar> > result;
	PZ_PTHREAD_MUTEX_LOCK(&fListAccessLock,"TPZDohrAssembleList::PopItem()");
	if (fWork.begin() != fWork.end()) {
		fNumItems--;
		result = *fWork.begin();
		fWork.pop_front();
	}
	PZ_PTHREAD_MUTEX_UNLOCK(&fListAccessLock,"TPZDohrAssembleList::PopItem()");
	return result;
}

template<class TVar>
void *TPZDohrAssembleList<TVar>::Assemble(void *voidptr)
{
	TPZDohrAssembleList *myptr = (TPZDohrAssembleList *) voidptr;
	while (myptr->fNumItems > 0) {
		TPZAutoPointer<TPZDohrAssembleItem<TVar> > work = myptr->PopItem();
		if (work) {
            PZ_PTHREAD_MUTEX_LOCK(&myptr->fAssemblyLock,"TPZDohrAssembleList::Assemble()");
			myptr->fAssembleIndexes->Assemble(work->fSubIndex,work->fAssembleData,*(myptr->fOutput));
            PZ_PTHREAD_MUTEX_UNLOCK(&myptr->fAssemblyLock,"TPZDohrAssembleList::Assemble()");
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
