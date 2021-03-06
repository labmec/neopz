/**
 * @file
 * @brief Contains the TPZDohrAssembleItem and TPZDohrAssembleList structs to assembling using Dohrman algorithm.
 * @since 20/07/10.
 */

#ifndef TPZDOHRASSEMBLELIST
#define TPZDOHRASSEMBLELIST

#include "pzfmatrix.h"
#include "tpzautopointer.h"
#include "tpzdohrassembly.h"
#include "TPZSemaphore.h"

#include <list>

/**
 * @ingroup substructure
 * @brief To assembling one item using Dohrmann algorithm. \ref substructure "Sub structure"
 */
template<class TVar>
struct TPZDohrAssembleItem {
	/** @brief Initialize the assemble item based on the submesh index and size of the local contribution */
	TPZDohrAssembleItem(int subindex, int size) : fSubIndex(subindex), fAssembleData(size,1,0.)
	{
	}
	/** @brief Initialize the assemble item based on the submesh index and size of the local contribution */
	TPZDohrAssembleItem(int subindex, int nrow, int ncol) : fSubIndex(subindex), fAssembleData(nrow,ncol,0.)
	{
	}
	/** @brief Substructure index */
	int fSubIndex;
	/** @brief The data which should be assembled */
	TPZFMatrix<TVar> fAssembleData;
};

/**
 * @ingroup substruture
 * @brief List of items to assembling using Dohrmann algorithm
 */ 
template<class TVar>
struct TPZDohrAssembleList {
	/** @brief Constructor indicating the number of items that will be assembled and the target matrix */
	TPZDohrAssembleList(int numitems, TPZFMatrix<TVar> &output, TPZAutoPointer<TPZDohrAssembly<TVar> > assembly);
	/// destructor
	~TPZDohrAssembleList();
	/** @brief The number of items that will be assembled before returning */
	int fNumItems;
	/** @brief Semaphore (to wake up assembly thread) */
	TPZSemaphore fSemaphore;
	/** @brief This is the mutex which controls the access to the list */
    std::mutex fListAccessLock;
	/** @brief This is the mutex which controls the assembly */
    std::mutex fAssemblyLock;
	/** @brief List of objects needed to be assembled */
	std::list<TPZAutoPointer<TPZDohrAssembleItem<TVar> > > fWork;
	/** @brief Add an item to the list in a thread safe way */
	void AddItem(TPZAutoPointer<TPZDohrAssembleItem<TVar> > assembleItem);
	/** @brief Remove an item from the list */
	TPZAutoPointer<TPZDohrAssembleItem<TVar> > PopItem();
	/** @brief Assembly indexes */
	TPZAutoPointer<TPZDohrAssembly<TVar> > fAssembleIndexes;
	/** @brief Target Matrix */
	TPZFMatrix<TVar> *fOutput;
	/** @brief Procedure which performs the assembly process */
	static void *Assemble(void *voidptr);
};

#endif
