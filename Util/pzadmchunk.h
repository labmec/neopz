/**
 * @file pzadmchunk.h
 * @brief Free store vector implementation.
 */
// $Id: pzadmchunk.h,v 1.1.1.1 2003-02-04 16:45:27 cantao Exp $

#ifndef PZADMCHUNK_H
#define PZADMCHUNK_H

#include "pzchunk.h"
#include "pzstack.h"
#include "pzerror.h"

/**
 * @ingroup util
 * @brief Implements a chunk vector with free store administration.
 *
 * An object of this class allows the user to request a new object of
 * the type administered and allows the user to flag given elements as
 * unused.
 */
template <class T>
class TPZAdmChunkVector : public TPZChunkVector<T>
{
   public :
      /**
       * Assignment operator.
       *
       * It will copy the objects from TPZAdmCh will call the
       * assignment operator on all objects (also the freed objects).
       *
       * @param TPZAdmCh Vector which will be duplicated.
       */
      TPZAdmChunkVector<T> & operator=(const TPZAdmChunkVector<T> &TPZAdmCh);

      /**
       * Copy constructor.
       *
       * @param AdmCh Object whose elements will be copied.
       */
      TPZAdmChunkVector(const TPZAdmChunkVector<T> &AdmCh);

      /**
       * Constructor.
       *
       * Constructor with indication of the initial size of the chunk
       * allocation vector and the size of the chunks these sizes
       * cannot be modified during the lifecycle of the object.
       *
       * @param numberofchunks Indicates how large the initial chunk
       * vector will be.
       *
       * @param chunkexponent Indicates the size of the chunks as an
       * exponent of 2.
       */
      TPZAdmChunkVector(int numberofchunks = DEFAULTNUMBEROFCHUNKS,
			int chunkexponent  = DEFAULTCHUNKEXPONENT);

      /**
       * Destructor
       */
      virtual ~TPZAdmChunkVector();

      /**
       * Makes more room for new elements.
       *
       * This method will search the list of free locations to return
       * the next free index in case there are no free indexes, this
       * method will increase the size of the chunk vector and returns
       * the allocated element.
       *
       * @return The index of a free element.
       */
      int AllocateNewElement();
   
      /**
       * Indicate an element as free.
       *
       * @note The object does not verify whether an element has been
       * freed several times.
       *
       * @param index The index of the element being put on the free stack.
       */
      void SetFree(int index);

      /**
       * Access method to return the number of free elements.
       *
       * @return Number of free elements.
       */
      inline int NFreeElements() { return fFree.NElements(); }

      /**
       * Sets the method to compact the data structure based on the
       * parameter type:
       *
       * <ul>
       * <li> when = 0 : never compact the data structure;
       * <li> when = 1 : compact the data structure now;
       * <li> when = 2 : compact the data structure always (default).
       * </ul>
       *
       * @param type Type of compacting scheme to be used.
       */
      void CompactDataStructure(int type=2);

      /**
       * Print index i into the fFree vector.
       */
      inline int PrintFree(int i)
      {
	 // Jorge 12/01/2000
	 return fFree[i];
      }

      /**
       * Increase the size of the chunk vector.
       *
       * @param newsize Requested new size of the vector.
       */
      void Resize(const int newsize);

   private:
      /** @shapeType DependencyLink */
      /*# TPZChunkVector lnkUnnamed */
      
      /**
       * Internal variable indicating the type of compacting scheme.
       *
       * @see CompactDataStructure.
       */
      int fCompactScheme;

      /** Number of free elements within each chunk. */
      TPZManVector<int> fNFree;

      /** List of indexes of freed elements. */
      TPZStack<int> fFree;
};

//--| IMPLEMENTATION |----------------------------------------------------------

template< class T >
TPZAdmChunkVector<T>::TPZAdmChunkVector(int numberofchunks,int chunkexponent)
   : TPZChunkVector<T>(numberofchunks,chunkexponent),
     fCompactScheme( 0 ),    // never compact the data structure
     fNFree( numberofchunks ),
     fFree()
{
   for( int i = 0; i < numberofchunks; i++ )
      fNFree[i] = 0;

   fNFree.Resize(0);
}

template< class T >
TPZAdmChunkVector<T>::~TPZAdmChunkVector() { }

// Return the index of a free element
template< class T >
int TPZAdmChunkVector<T>::AllocateNewElement() {
   if(fFree.NElements() >0)
   {
      int index = fFree.Pop(),chunk;
      chunk = index >> fExponent;
      fNFree[chunk]--;
      return index;
   }

   Resize(NElements()+1);

   return NElements()-1;
}

// Indicate an element as free
template< class T >
void TPZAdmChunkVector<T>::SetFree(int index) {
#ifndef NODEBUG
   if(index<0) {
      PZError << "TPZAdmChunkVector::SetFree. Bad parameter index." << endl;
      PZError.flush();
      return;
   }
#endif

   int chunk = index >> fExponent;

   fNFree[chunk]++;
   fFree.Push(index);

   if(fCompactScheme == 2)
      CompactDataStructure(2);
}

// Let to compact the data structure.
// If type=0 never compact, (type=1 let to compact now),(type=2  let to compact always)
template< class T >
void TPZAdmChunkVector<T>::CompactDataStructure(int type) {
#ifndef NOTDEBUG
   if(type<0) {
      PZError << "TPZAdmChunkVector::CompactDataStructure. Bad parameter type."
	      << endl;
      PZError.flush();

      return;
   }
#endif

   if(type == 0) {
      fCompactScheme = 0;
      return;
   }
   if(type == 2) fCompactScheme = 2;
   int i,chunksize = 1<<fExponent;
   int nchunksused = 0;
   if(NElements()) nchunksused = ((NElements()-1) >> fExponent)+1;
   i = nchunksused-1;
   int maxfree = NElements()-((nchunksused-1) << fExponent);

   if(i>=0 && fVec[i] && fNFree[i] == maxfree) {
      Resize(chunksize*i);
      i--;
      while(i >=0 && fVec[i] && fNFree[i] == chunksize) {
	 Resize(chunksize*i);
	 i--;
      }
   }
   fVec.Shrink();
   fNFree.Shrink();
   fFree.Shrink();
}

template < class T >
TPZAdmChunkVector<T>::TPZAdmChunkVector(const TPZAdmChunkVector<T> &AdmCh) :
   TPZChunkVector<T>(AdmCh), fCompactScheme( AdmCh.fCompactScheme ),
   fNFree(AdmCh.fNFree), fFree(AdmCh.fFree)
{
   // NOTHING TO DO HERE!
}

template < class T >
TPZAdmChunkVector<T> & TPZAdmChunkVector<T>::operator=(
   const TPZAdmChunkVector<T> &AdmCh)
{
   if(this == &AdmCh)
      return *this;

   TPZChunkVector<T>::operator=(AdmCh);

   fFree = AdmCh.fFree;
   fNFree = AdmCh.fNFree;
   fCompactScheme = AdmCh.fCompactScheme;

   return *this;
}

template< class T >
void TPZAdmChunkVector<T>::Resize(const int newsize) {
#ifndef NOTDEBUG
   if(newsize<0) {
      PZError << "TPZAdmChunkVector::Resize. Bad parameter newsize." << endl;
      PZError.flush();
      return;
   }
#endif

   TPZChunkVector<T>::Resize(newsize);

   int sizechunk = 1 << fExponent;
   int nchunks = fNFree.NElements();
   int chunksneeded = fVec.NElements(); // equivalent to newsize>>fExponent??

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
	 int chunk = index/sizechunk;

	 if(chunk == chunksneeded-1)
	    fNFree[chunksneeded-1]--;
      }
   }

   fFree = temp;
}

#endif // PZADMCHUNK_H

//--| PZ |----------------------------------------------------------------------
