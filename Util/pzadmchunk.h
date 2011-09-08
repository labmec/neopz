/**
 * @file pzadmchunk.h
 * @brief Free store vector implementation.
 */
// $Id: pzadmchunk.h,v 1.6 2005-04-25 02:55:51 phil Exp $

#ifndef PZADMCHUNK_H
#define PZADMCHUNK_H

#include "pzchunk.h"
#include "pzstack.h"
#include "pzerror.h"

class TPZSaveable;

/**
 * @ingroup util
 * @brief Implements a chunk vector with free store administration. \ref util "Utility"
 */
/** 
 * An object of this class allows the user to request a new object of
 * the type administered and allows the user to flag given elements as
 * unused.
 */
template <class T,int EXP=10>
class TPZAdmChunkVector : public TPZChunkVector<T,EXP>
{
	public :
	/** 
	 * @brief Assignment operator. 
	 * @param TPZAdmCh Vector which will be duplicated.
	 */
	/**
	 * It will copy the objects from TPZAdmCh will call the
	 * assignment operator on all objects (also the freed objects).
	 */
	TPZAdmChunkVector<T,EXP> & operator=(const TPZAdmChunkVector<T,EXP> &TPZAdmCh);
	
	/**
	 * @brief Copy constructor.
	 * @param AdmCh Object whose elements will be copied.
	 */
	TPZAdmChunkVector(const TPZAdmChunkVector<T,EXP> &AdmCh);
	
	/**
	 * Constructor with indication of the initial size of the chunk
	 * allocation vector and the size of the chunks these sizes
	 * cannot be modified during the lifecycle of the object.
	 */
	/**
	 * @brief Constructor.
	 * @param numberofchunks Indicates how large the initial chunk
	 * vector will be.
	 */
	TPZAdmChunkVector(int numberofchunks = DEFAULTNUMBEROFCHUNKS);
	
	/** @brief Destructor */
	virtual ~TPZAdmChunkVector();
	
	/**
	 * @brief Makes more room for new elements.
	 * @return The index of a free element.
	 */
	/** 
	 * This method will search the list of free locations to return
	 * the next free index \n in case there are no free indexes, this
	 * method will increase the size of the chunk vector \n and returns
	 * the allocated element.
	 */
	int AllocateNewElement();
	
	/** 
	 * @brief Indicate an element as free.
	 * @note The object does not verify whether an element has been freed several times.
	 * @param index The index of the element being put on the free stack.
	 */
	void SetFree(int index);
	
	/**
	 * @brief Access method to return the number of free elements.
	 * @return Number of free elements.
	 */
	inline int NFreeElements() { return fFree.NElements(); }
	
	/**
	 * @brief Sets the method to compact the data structure based on the
	 * @param type Type of compacting scheme to be used.
	 */
	/**
	 * parameter type:
	 * <ul>
	 * <li> when = 0 : never compact the data structure;
	 * <li> when = 1 : compact the data structure now;
	 * <li> when = 2 : compact the data structure always (default).
	 * </ul>
	 */
	void CompactDataStructure(int type=2);
	
	/** Print index i into the fFree vector. */
	inline int PrintFree(int i)
	{
		// Jorge 12/01/2000
		return fFree[i];
	}
	
	/**
	 * @brief  Increase the size of the chunk vector.
	 * @param newsize Requested new size of the vector.
	 */
	void Resize(const int newsize);
	
private:
	
	friend class TPZSaveable;
	
	/**
	 * @brief Internal variable indicating the type of compacting scheme.
	 * @see CompactDataStructure.
	 */
	int fCompactScheme;
	
	/** @brief Number of free elements within each chunk. */
	TPZManVector<int> fNFree;
	
	/** @brief List of indexes of freed elements. */
	TPZStack<int> fFree;
};

//--| IMPLEMENTATION |----------------------------------------------------------

template< class T , int EXP>
TPZAdmChunkVector<T,EXP>::TPZAdmChunkVector(int numberofchunks)
: TPZChunkVector<T,EXP>(numberofchunks),
fCompactScheme( 0 ),    // never compact the data structure
fNFree( numberofchunks ),
fFree()
{
	for( int i = 0; i < numberofchunks; i++ )
		fNFree[i] = 0;
	
	fNFree.Resize(0);
}

template< class T, int EXP >
TPZAdmChunkVector<T,EXP>::~TPZAdmChunkVector() { }

// Return the index of a free element
template< class T,int EXP >
int TPZAdmChunkVector<T,EXP>::AllocateNewElement() {
	if(fFree.NElements() >0)
	{
		int index = fFree.Pop(),chunk;
		chunk = index >> EXP;
		fNFree[chunk]--;
		return index;
	}
	
	Resize(this->NElements()+1);
	
	return this->NElements()-1;
}

// Indicate an element as free
template< class T,int EXP >
void TPZAdmChunkVector<T,EXP>::SetFree(int index) {
#ifndef NODEBUG
	if(index<0) {
		PZError << "TPZAdmChunkVector::SetFree. Bad parameter index." << std::endl;
		PZError.flush();
		return;
	}
#endif
	
	int chunk = index >> EXP;
	
	fNFree[chunk]++;
	fFree.Push(index);
	
	if(fCompactScheme == 2)
		CompactDataStructure(2);
}

// Let to compact the data structure.
// If type=0 never compact, (type=1 let to compact now),(type=2  let to compact always)
template< class T, int EXP >
void TPZAdmChunkVector<T,EXP>::CompactDataStructure(int type) {
#ifndef NODEBUG
	if(type<0) {
		PZError << "TPZAdmChunkVector::CompactDataStructure. Bad parameter type."
		<< std::endl;
		PZError.flush();
		
		return;
	}
#endif
	
	if(type == 0) {
		fCompactScheme = 0;
		return;
	}
	if(type == 2) fCompactScheme = 2;
	int i,chunksize = 1<<EXP;
	int nchunksused = 0;
	if(this->NElements()) nchunksused = ((this->NElements()-1) >> EXP)+1;
	i = nchunksused-1;
	int maxfree = this->NElements()-((nchunksused-1) << EXP);
	
	if(i>=0 && this->fVec[i] && fNFree[i] == maxfree) {
		Resize(chunksize*i);
		i--;
		while(i >=0 && this->fVec[i] && fNFree[i] == chunksize) {
			Resize(chunksize*i);
			i--;
		}
	}
	this->fVec.Shrink();
	fNFree.Shrink();
	fFree.Shrink();
}

template < class T,int EXP >
TPZAdmChunkVector<T,EXP>::TPZAdmChunkVector(const TPZAdmChunkVector<T,EXP> &AdmCh) :
TPZChunkVector<T,EXP>(AdmCh), fCompactScheme( AdmCh.fCompactScheme ),
fNFree(AdmCh.fNFree), fFree(AdmCh.fFree)
{
	// NOTHING TO DO HERE!
}

template < class T , int EXP>
TPZAdmChunkVector<T,EXP> & TPZAdmChunkVector<T,EXP>::operator=(
															   const TPZAdmChunkVector<T,EXP> &AdmCh)
{
	if(this == &AdmCh)
		return *this;
	
	TPZChunkVector<T,EXP>::operator=(AdmCh);
	
	fFree = AdmCh.fFree;
	fNFree = AdmCh.fNFree;
	fCompactScheme = AdmCh.fCompactScheme;
	
	return *this;
}

template< class T, int EXP >
void TPZAdmChunkVector<T,EXP>::Resize(const int newsize) {
#ifndef NODEBUG
	if(newsize<0) {
		PZError << "TPZAdmChunkVector::Resize. Bad parameter newsize." << std::endl;
		PZError.flush();
		return;
	}
#endif
	
	TPZChunkVector<T,EXP>::Resize(newsize);
	
	//   int sizechunk = 1 << EXP;
	int nchunks = fNFree.NElements();
	int chunksneeded = this->fVec.NElements(); // equivalent to newsize>>fExponent??
	
	fNFree.Resize(chunksneeded);
	
	for( int i = nchunks; i < chunksneeded; i++ )
		fNFree[i] = 0;
	
	if(chunksneeded > nchunks) return;
	
	// delete all free indexes which are above the new size
	// update the number of free elements of the last chunk
	TPZStack<int> temp(fFree);
	temp.Resize(0);
	
	while(fFree.NElements() > 0) {
		int index = fFree.Pop();
		if(index < newsize)
			temp.Push(index);
		else {
			int chunk = index >> EXP;
			
			if(chunk == chunksneeded-1)
				fNFree[chunksneeded-1]--;
		}
	}
	
	fFree = temp;
}

#endif // PZADMCHUNK_H

//--| PZ |----------------------------------------------------------------------
