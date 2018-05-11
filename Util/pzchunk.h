/**
 * @file
 * @brief Free store vector implementation in chunks.
 */

#ifndef PZCHUNK_H
#define PZCHUNK_H

#include "pzmanvector.h"
#include "pzerror.h"

#include <stdlib.h>

/**
 * @ingroup util
 * @brief Default number of elements which will be allocated in the chunk vector.
 */
const int DEFAULTNUMBEROFCHUNKS = 100;

/**
 * @ingroup util
 * @brief Default number of elements in each chunk is \f$ pow(2,--) \f$. 
 */
const int DEFAULTCHUNKEXPONENT = 10;

/**
 * @ingroup util
 * @brief An object of this class implements a vector which allocates objects by chunks. \ref util "Utility"
 */
/**
 * The expansion of a TChunkVector object does not
 * involve the copying of the already allocated objects.
 */
template<class T,int EXP=10>
class TPZChunkVector
{
	public :
	/**
	 * @brief Assignment operator, copies all elements from the object TCh.
	 * @param TCh Chunk vector from which the elements will be copied.
	 */
	TPZChunkVector<T,EXP> & operator=(const TPZChunkVector<T,EXP> &TCh);
	
	/**
	 * @brief Copy constructor.
	 * @param TCh The elements of the current object will be copied from TCh.
	 */
	TPZChunkVector(const TPZChunkVector<T,EXP> &TCh);
	
	/**
	 * @brief Constructor with numberofchuncks chuncks.
	 * @param numberofchunks Indicates how large the initial chunk vector will be.
	 */
	TPZChunkVector(int64_t numberofchunks = DEFAULTNUMBEROFCHUNKS);
	/** @brief Destructor. */
	virtual ~TPZChunkVector();
	
	/**
	 * @brief Access method to query the number of elements of the vector.
	 * @return Returns the number of elements of the vector.
	 */
	inline int64_t NElements() const {return fNElements;}
	
	/**
	 * @brief Increase the size of the chunk vector.
	 * @param newsize New size of the vector. Does not indicate how
	 * much memory will be allocated!
	 */
	void Resize(const int64_t newsize);
	
	/** @brief Returns a reference to the ith element of the vector. */
	/**
	 * If NODEBUG is defined, the element is returned without any
	 * bounds checking, else a bounds check is performed and the \n
	 * code exits if the index is out of bounds.  Note that NODEBUG
	 * may modify the result of the code.
	 */
	T &operator[](const int64_t nelem) const;
	
	/** @brief Finds the index of an object by its pointer */
	int64_t FindObject(T *object);
	
protected:
	/** @brief Number of elements of the chunk vector. */
	int64_t fNElements;
	
	/** @brief Vector which points to each chunk of objects. */
	TPZManVector<T*> fVec;
};

//--| IMPLEMENTATION |----------------------------------------------------------

template< class T, int EXP >
TPZChunkVector<T,EXP>::TPZChunkVector(int64_t numberofchunks) :
fVec(numberofchunks)
{
	fNElements = 0;
	
	for(int64_t i=0; i<numberofchunks; i++)
		fVec[i] = 0;
	
	fVec.Resize(0);
}

template< class T,int EXP >
TPZChunkVector<T,EXP>::~TPZChunkVector()
{
	int64_t nchunks=fVec.NElements();
	
	for(int64_t i=0;i<nchunks;i++)
		if(fVec[i]) delete[] fVec[i];
}

// Increase the size of the chunk vector

template< class T ,int EXP>
void TPZChunkVector<T,EXP>::Resize(const int64_t newsize) {
#ifndef NODEBUG
	if(newsize<0) {
		PZError << "TPZChunkVector::Resize. Bad parameter newsize." << std::endl;
		PZError.flush();
		return;
	}
#endif
	fNElements = newsize;
	
	int64_t i;
	int64_t nchunks = fVec.NElements();
	
	int64_t chunksneeded = ((newsize-1) >> EXP)+1;
	if (chunksneeded == nchunks) return;
	
	T *NullPointer = 0;
	
	if(chunksneeded > nchunks)
		fVec.Resize(chunksneeded,NullPointer);
	
	const int64_t sizechunk = 1 << EXP;
	for(i = 0; i<chunksneeded; i++ )
		if(!fVec[i])
			fVec[i] = new T[sizechunk];
	
	for( ; i<nchunks; i++) if(fVec[i]) {
		delete [] fVec[i];
		fVec[i] = 0;
	}
	
	fVec.Resize(chunksneeded);
}

template<class T, int EXP>
int64_t TPZChunkVector<T,EXP>::FindObject(T *obj) {
	int64_t nch = fVec.NElements();
	int64_t ich;
	int64_t index = 0;
	// number of elements in a chunk
	int64_t nelch = 1<<EXP;
	for(ich=0; ich<nch; ich++) {
		if(fVec[ich] == 0) continue;
		if(obj >= fVec[ich] && obj < fVec[ich]+nelch) {
			return index + (obj-fVec[ich]);
		} else {
			index += nelch;
		}
	}
	if(ich == nch) return -1;
	return index;
}

// Return a reference to the ith element of the vector
template< class T, int EXP >
inline T &TPZChunkVector<T,EXP>::operator[](const int64_t nelem) const {
#ifndef NODEBUG
	if(nelem<0 || nelem >= NElements()) {
		PZError << "TPZChunkVector::operator[]. "
		<< "Bad parameter nelem." << nelem << " NElements "
		<< NElements() << std::endl;
		PZError.flush();
		DebugStop();
		exit (-1);
		return fVec[0][0];
	}
#endif
	
	const int64_t mask = (1 << EXP)-1;
	
	// nelem & mask deve ser nelem ??
	return (fVec[nelem >> EXP])[nelem & mask];
}

template< class T, int EXP >
TPZChunkVector<T,EXP>::TPZChunkVector(const TPZChunkVector<T,EXP> &TCh) :
fVec(TCh.fVec.NElements()) {
	
	fNElements = TCh.NElements();
	int64_t nchunks = TCh.fVec.NElements();
	
	for(int64_t i=0; i<nchunks; i++) {
		if(!TCh.fVec[i]) {
			fVec[i] = 0;
		} else {
			int64_t j,k=1<<EXP;
			T* ptr = new  T[1 << EXP];
			fVec[i] = ptr;
			T* ptrcp = TCh.fVec[i];
			for(j=0;j<k;j++) ptr[j] = ptrcp[j];
		}
	}
}

template < class T, int EXP >
TPZChunkVector<T,EXP> & TPZChunkVector<T,EXP>::operator=(const TPZChunkVector<T,EXP> &TCh){
	if(this == &TCh) return *this;
	int64_t prvsz = fVec.NElements();
	int64_t sz = TCh.fVec.NElements();
	int64_t i;
	for(i=sz;i<prvsz;i++)
		if(fVec[i])
			delete[] fVec[i];
	
	fVec.Resize(sz,0);
	fVec.Shrink();
	fNElements = TCh.NElements();
	
	int64_t nchunks = sz;
	
	for(i=0; i<nchunks; i++) {
		if(!TCh.fVec[i]) {
			fVec[i] = 0;
		} else {
			int64_t j,k=1<<EXP;
			fVec[i] = new T[1 << EXP];
			
			for(j=0;j<k;j++)
				fVec[i][j] = TCh.fVec[i][j];
		}
	}
	
	return *this;
}

#endif // PZCHUNK_H

//--| PZ |----------------------------------------------------------------------
