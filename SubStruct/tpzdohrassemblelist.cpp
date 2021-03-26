/**
 * @file
 * @brief Contains the implementation of the TPZDohrAssembleList methods.
 */

#include "tpzdohrassemblelist.h"


template<class TVar>
TPZDohrAssembleList<TVar>::TPZDohrAssembleList(int numitems, TPZFMatrix<TVar> &output, TPZAutoPointer<TPZDohrAssembly<TVar> > assembly) : fNumItems(numitems),
fAssembleIndexes(assembly), fOutput(&output), fAssemblyLock(), fListAccessLock()
{
}

template<class TVar>
TPZDohrAssembleList<TVar>::~TPZDohrAssembleList()
{
}

// Add an item to the list in a thread safe way
template<class TVar>
void TPZDohrAssembleList<TVar>::AddItem(TPZAutoPointer<TPZDohrAssembleItem<TVar> > assembleItem)
{
    std::lock_guard<std::mutex> lock(fListAccessLock);
	fWork.push_back(assembleItem);
	fSemaphore.Post();
}
// remove an item from the list
template<class TVar>
TPZAutoPointer<TPZDohrAssembleItem<TVar> > TPZDohrAssembleList<TVar>::PopItem()
{
	TPZAutoPointer<TPZDohrAssembleItem<TVar> > result;
    std::lock_guard<std::mutex> lock(fListAccessLock);
	if (fWork.begin() != fWork.end()) {
		fNumItems--;
		result = *fWork.begin();
		fWork.pop_front();
	}
	return result;
}

template<class TVar>
void *TPZDohrAssembleList<TVar>::Assemble(void *voidptr)
{
	TPZDohrAssembleList *myptr = (TPZDohrAssembleList *) voidptr;
	while (myptr->fNumItems > 0) {
		TPZAutoPointer<TPZDohrAssembleItem<TVar> > work = myptr->PopItem();
		if (work) {
            std::lock_guard<std::mutex> lock(myptr->fAssemblyLock);
			myptr->fAssembleIndexes->Assemble(work->fSubIndex,work->fAssembleData,*(myptr->fOutput));
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
