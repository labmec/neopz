// -*- c++ -*-
/**
 * @file pzadmchunk.h
 * @brief Free store vector implementation.
 */
// $Id: pzchunk.h,v 1.2 2003-10-05 00:31:56 phil Exp $

#ifndef PZCHUNK_H
#define PZCHUNK_H

#include "pzmanvector.h"
#include "pzerror.h"

using namespace std;

/** Default number of elements which will be allocated in the chunk vector. */
const int DEFAULTNUMBEROFCHUNKS = 100;

/** Default number of elements in each chunk is pow(2,--). */
const int DEFAULTCHUNKEXPONENT = 10;

/**
 * @ingroup util
 *
 * @brief An object of this class implements a vector which allocates
 * objects by chunks. The expansion of a TChunkVector object does not
 * involve the copying of the already allocated objects.
 */
template<class T>
class TPZChunkVector
{
   public :
      /**
       * Assignment operator, copies all elements from the object TCh.
       *
       * @param TCh Chunk vector from which the elements will be
       * copied.
       */
      TPZChunkVector<T> & operator=(const TPZChunkVector<T> &TCh);

      /**
       * Copy constructor.
       *
       * @param TCh The elements of the current object will be copied
       * from TCh.
       */
      TPZChunkVector(const TPZChunkVector<T> &TCh);

      /**
       * Access method to query the number of elements of the vector.
       *
       * @return Returns the number of elements of the vector.
       */
      inline int NElements() const {return fNElements;}

      /**
       * @param numberofchunks Indicates how large the initial chunk
       * vector will be.
       *
       * @param chunkexponent Indicates the size of the chunks as an
       * exponent of 2.
       */
      TPZChunkVector(int numberofchunks = DEFAULTNUMBEROFCHUNKS,
		     int chunkexponent  = DEFAULTCHUNKEXPONENT);
      /**
       * Destructor.
       */
      virtual ~TPZChunkVector();

      /**
       * Increase the size of the chunk vector.
       *
       * @param newsize New size of the vector. Does not indicate how
       * much memory will be allocated!
       */
      virtual void Resize(const int newsize);

      /**
       * Returns a reference to the ith element of the vector.
       *
       * If NODEBUG is defined, the element is returned without any
       * bounds checking, else a bounds check is performed and the
       * code exits if the index is out of bounds.  Note that NODEBUG
       * may modify the result of the code.
       */
      T &operator[](const int nelem) const;

      /**
       * Finds the index of an object by its pointer
       */
      int FindObject(T *object);

   protected:
      /** Number of elements of the chunk vector. */
      int fNElements;


      /** Exponent which determines the size of each chunk. */
      int fExponent;

      /** Vector which points to each chunk of objects. */
      TPZManVector<T*> fVec;
};

//--| IMPLEMENTATION |----------------------------------------------------------

template< class T >
TPZChunkVector<T>::TPZChunkVector(int numberofchunks, int chunkexponent) :
   fVec(numberofchunks)
{
#ifndef NOTDEBUG
   if(chunkexponent<0) {
      PZError << "TPZChunkVector constructor. "
	      << "Bad parameter chunkexponent, then chunkexponent = 0."
	      << endl;

      PZError.flush();
      chunkexponent = 0;
   }
#endif

   fExponent = chunkexponent;
   fNElements = 0;

   for(int i=0; i<numberofchunks; i++)
      fVec[i] = 0;

   fVec.Resize(0);
}

template< class T >
TPZChunkVector<T>::~TPZChunkVector()
{
   int nchunks=fVec.NElements();

   for(int i=0;i<nchunks;i++)
      if(fVec[i]) delete[] fVec[i];
}

// Increase the size of the chunk vector

template< class T >
void TPZChunkVector<T>::Resize(const int newsize) {
#ifndef NOTDEBUG
   if(newsize<0) {
      PZError << "TPZChunkVector::Resize. Bad parameter newsize." << endl;
      PZError.flush();
      return;
   }
#endif

   int i,sizechunk = 1 << fExponent;
   int nchunks = fVec.NElements();

   // equivalent to newsize>>fExponent??
   int chunksneeded = ((newsize-1)/sizechunk)+1;

   if(newsize == 0)
      chunksneeded = 0;

   T *NullPointer = 0;

   if(chunksneeded > nchunks)
      fVec.Resize(chunksneeded,NullPointer);

   for(i = 0; i<chunksneeded; i++ )
      if(!fVec[i])
	 fVec[i] = new T[sizechunk];

   for( ; i<nchunks; i++) if(fVec[i]) {
      delete [] fVec[i];
      fVec[i] = 0;
   }
   
   fVec.Resize(chunksneeded);
   fNElements = newsize;
}

template<class T>
int TPZChunkVector<T>::FindObject(T *obj) {
  int nch = fVec.NElements();
  int ich;
  int index = 0;
  // number of elements in a chunk
  int nelch = 1<<fExponent;
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
template< class T >
T &TPZChunkVector<T>::operator[](const int nelem) const {
#ifndef NODEBUG
   if(nelem<0 || nelem >= NElements()) {
      PZError << "TPZChunkVector::operator[]. "
	      << "Bad parameter nelem." << nelem << " NElements "
	      << NElements() << endl;
      PZError.flush();
      exit (-1);
      return fVec[0][0];
   }
#endif

   int mask = (1 << fExponent)-1;

   // nelem & mask deve ser nelem ??
   return (fVec[nelem >> fExponent])[nelem & mask];
}

template< class T >
TPZChunkVector<T>::TPZChunkVector(const TPZChunkVector<T> &TCh) :
   fVec(TCh.fVec.NElements()) {
   fExponent = TCh.fExponent;
   fNElements = TCh.NElements();
   int nchunks = TCh.fVec.NElements();

   for(int i=0; i<nchunks; i++) {
      if(!TCh.fVec[i]) {
	 fVec[i] = 0;
      } else {
	 int j,k=1<<fExponent;
	 fVec[i] = new T[1 << TCh.fExponent];
	 for(j=0;j<k;j++) fVec[i][j] = TCh.fVec[i][j];
      }
   }
}


template < class T >
TPZChunkVector<T> & TPZChunkVector<T>::operator=(const TPZChunkVector<T> &TCh){
   if(this == &TCh) return *this;
   int n=fVec.NElements();

   for(int i=0;i<n;i++)
      if(fVec[i])
	 delete[] fVec[i];
   fVec.Resize(TCh.fVec.NElements());
   fVec.Shrink();
   fExponent = TCh.fExponent;
   fNElements = TCh.NElements();

   int nchunks = TCh.fVec.NElements();

   for(int i=0; i<nchunks; i++) {
      if(!TCh.fVec[i]) {
	 fVec[i] = 0;
      } else {
	 int j,k=1<<fExponent;
	 fVec[i] = new T[1 << TCh.fExponent];

	 for(j=0;j<k;j++)
	    fVec[i][j] = TCh.fVec[i][j];
      }
   }

   return *this;
}

#endif // PZCHUNK_H

//--| PZ |----------------------------------------------------------------------
