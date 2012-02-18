/**
 * @file
 * @brief Free store vector implementation.
 */

#ifndef PZMANVECTOR_H
#define PZMANVECTOR_H

#include "pzvec.h"
#include "pzerror.h"

/// To allocate vector by default
const int DEFAULTVEC_ALLOC = 200;

/**
 * @ingroup util
 * @see Shrink
 * @brief Implements a vector class which allows to use external storage provided by the user. \ref util "Utility"
 */
/** 
 * The external storage will be used if the number of elements of the external storage is \n
 * greater than or equal to the number of elements the object needs. When changing the size \n
 * of the object, this class will only allocate new storage area if the storage currently allocated \n
 * is insufficient to the hold the object. \n
 * This makes the resize method more efficient in terms of dynamic memory allocation the Shrink method \n
 * will reallocate the storage to fit the number of elements exactly.
 */
template < class T, int NumExtAlloc = DEFAULTVEC_ALLOC >
class TPZManVector : public TPZVec< T >
{
public:
	/**
	 * @brief Creates a vector of a given size.
	 * @param size Size of the new vector.
	 */
	/**
	 * It will call the empty constructor on all objects of type T.
	 */
	TPZManVector( const int size = 0 );
	
	/**
	 * @brief Creates a vector of a given size, filling it.
	 * @param size Size of the new vector.
	 * @param copy Model object to initialize the other objects.
	 */
	/**
	 * It will call the empty constructor on all objects of type T
	 * created. Copies the object copy to all elements.
	 */
	TPZManVector( const int size, const T& copy );
	
	/**
	 * @brief Copy constructor.
	 * @param copy Original vector.
	 */
	/** It will call the empty constructor on all objects of type T created. */
	TPZManVector( const TPZManVector< T, NumExtAlloc >& copy );
	
	TPZManVector(const TPZVec<T> & copy );
	
	/**
	 * @brief Assignment operator.
	 * @param copy Vector which will be copied.
	 * @return Reference to the current object.
	 */
	/**
	 * It first deletes the allocated storage before allocating
	 * storage for the copy only if necessary (when there is no
	 * preallocated storage or when the current storage cannot hold
	 * the copied vector
	 */
	TPZManVector< T, NumExtAlloc >& operator=(
											  const TPZManVector< T, NumExtAlloc >& copy);
	
	/** @brief Destructor. */
	/** Deletes the storage allocated. */
	virtual ~TPZManVector();
	
	/** @brief Returns number of elements allocated for this object. */
	inline int NAlloc() const { return fNAlloc; }
	
	/**
	 * @brief Expands the allocated storage to fit the newsize parameter.
	 * @param newsize Storage size which is requested.
	 */ 
	/** Does nothing if the externally provided storage is larger than newsize. Does not change the size of the vector. */
	void Expand( const int newsize );
	
	/** @brief It reallocates storage to fit the necessary storage exactly. */
	void Shrink();
	
	/**
	 * @brief Resizes the vector object.
	 * @see TPZVec<T>::Resize(const int, const T&)
	 * @param newsize Size of the vector.
	 * @param object Object used to initialize the new members.
	 */
	/**
	 * It reallocates storage if necessary, and copies the existing
	 * objects onto the new storage.
	 */
	virtual void Resize( const int newsize, const T& object );
	
	/**
	 * @brief Resizes the vector object reallocating the storage if necessary.
	 * @see TPZVec<T>::Resize(const int);
	 * @param newsize Size of the vector.
	 */
	/**
	 * It copies the existing objects to the new storage. The new
	 * members are not initialized.
	 */
	virtual void Resize( const int newsize );
	
private:
	
	/**
	 * @brief Returns a suggested size for expanding the storage to fit the required storage.
	 * @param proposed Storage needed for the new vector.
	 * @return Expansion size which is suggested.
	 */
	/**
	 * It suggests to expand the storage size at least 20%. This
	 * method will not expand the allocated storage size.
	 */
	int ExpandSize( const int proposed ) const;
	
	/** @brief Number of elements allocated for this object. */
	int fNAlloc;
	
	/** @brief Pointer to the externally allocated space. */
	T fExtAlloc[ NumExtAlloc ];
};

//--| IMPLEMENTATION |----------------------------------------------------------

template< class T, int NumExtAlloc >
TPZManVector< T, NumExtAlloc >::TPZManVector( const int size ) :
TPZVec<T>( 0 ) // There is always some static allocation.
{
	/* If the size requested fits inside the size already provided
	 * statically, provides access to that space, by setting some
	 * TPZVec data.
	 */
	if( size <= NumExtAlloc )
	{	
		// Needed to make TPZVec::operator[] work properly.
		this->fStore     = fExtAlloc;
		this->fNElements = size;
		// No memory was allocated by the constructor.
		fNAlloc    = 0;
	}
	else // The size requested is bigger than the size already provided.
	{
		// Executes the allocation that would be done by TPZVec<T>(size).
		this->fStore     = new T[ size ];
		this->fNElements = size;
		fNAlloc    = size;
	}
}

template< class T, int NumExtAlloc >
TPZManVector< T, NumExtAlloc >::TPZManVector( const int size, const T& copy ) :
TPZVec<T>( 0 )
//fNExtAlloc( NumExtAlloc ) // There is always some static allocation.
{
	/* If the size requested fits inside the size already provided
	 * statically, provides access to that space, by setting some
	 * TPZVec data.
	 */
	if( size <= NumExtAlloc )
	{	
		// Needed to make TPZVec::operator[] work properly.
		this->fStore     = fExtAlloc;
		this->fNElements = size;
		// No memory was allocated by the constructor.
		fNAlloc    = 0;
	}
	else // The size requested is bigger than the size already provided.
	{
		// Executes the allocation that would be done by TPZVec<T>(size).
		this->fStore     = new T[ size ];
		this->fNElements = size;
		fNAlloc    = size;
	}
	
	for( int i = 0; i < size; i++ )
	{
		this->fStore[i] = copy;
	}
}

template< class T, int NumExtAlloc>
inline TPZManVector< T, NumExtAlloc >::TPZManVector(
													const TPZManVector< T, NumExtAlloc >& copy )
{
	const int size = copy.NElements();
	
	/* If the size requested fits inside the size already provided
	 * statically, provides access to that space, by setting some
	 * TPZVec data.
	 */
	if( size <= (int) (sizeof(fExtAlloc)/sizeof(T)))
	{
		// Needed to make TPZVec::operator[] work properly.
		this->fStore     = fExtAlloc;
		this->fNElements = size;
		// No memory was allocated by the constructor.
		fNAlloc    = 0;
	}
	else // The size requested is bigger than the size already provided.
	{
		// Executes the allocation that would be done by TPZVec<T>(size).
		this->fStore     = new T[ size ];
		this->fNElements = size;
		fNAlloc    = size;
	}
	
	for( int i = 0; i < size; i++ )
	{
		this->fStore[i] = copy.fStore[i];
	}
}

template< class T, int NumExtAlloc>
inline TPZManVector< T, NumExtAlloc >::TPZManVector(
													const TPZVec<T> & copy )
{
	const int size = copy.NElements();
	
	/* If the size requested fits inside the size already provided
	 * statically, provides access to that space, by setting some
	 * TPZVec data.
	 */
	if( size <= (int) (sizeof(fExtAlloc)/sizeof(T)))
	{
		// Needed to make TPZVec::operator[] work properly.
		this->fStore     = fExtAlloc;
		this->fNElements = size;
		// No memory was allocated by the constructor.
		fNAlloc    = 0;
	}
	else // The size requested is bigger than the size already provided.
	{
		// Executes the allocation that would be done by TPZVec<T>(size).
		this->fStore     = new T[ size ];
		this->fNElements = size;
		fNAlloc    = size;
	}
	
	for( int i = 0; i < size; i++ )
	{
		this->fStore[i] = copy[i];
	}
}

template< class T, int NumExtAlloc >
TPZManVector< T, NumExtAlloc >& TPZManVector< T, NumExtAlloc >::operator=(
																		  const TPZManVector< T, NumExtAlloc >& copy )
{
	// Checking auto assignment.
	if( this == &copy )
	{
		return *this;
	}
	
	const int nel = copy.NElements();
	
	if( nel > fNAlloc && this->fStore && this->fStore != fExtAlloc )
	{
		delete [] this->fStore;		
		this->fStore  = 0;
		fNAlloc = 0;
	}
	
	if(nel <= NumExtAlloc)
	{
		if(this->fStore != fExtAlloc)
		{
			delete []this->fStore;
		}
		this->fNAlloc = 0;
		this->fStore = fExtAlloc;
		this->fNElements = nel;		
	}
    else if( fNAlloc >= nel )
	{
		this->fNElements = nel;
	}
	else
	{
		this->fStore     = new T[ nel ];
		fNAlloc    = nel;
		this->fNElements = nel;
	}
	
	for( int i = 0; i < nel; i++ )
	{
		this->fStore[i] = copy.fStore[i];
	}
	
	return *this;
}

template< class T, int NumExtAlloc >
TPZManVector< T, NumExtAlloc >::~TPZManVector()
{
	if( this->fStore == fExtAlloc )
	{
		this->fStore = 0;
	}
	
	//   fNExtAlloc = 0;
	fNAlloc    = 0;
}

template< class T, int NumExtAlloc >
void TPZManVector< T, NumExtAlloc >::Expand( const int newsize )
{
	// If newsize is negative then return.
	if( newsize <= fNAlloc || newsize <= NumExtAlloc )
	{
		return;
	}
	else
	{// the size is larger than the allocated memory
		T* newstore = new T[ newsize ];
		
		for ( int i = 0; i < this->fNElements; i++)
		{
			newstore[i] = this->fStore[i];
		}
		
		if( this->fStore != fExtAlloc )
		{
			delete [] this->fStore;
		}
		
		this->fStore  = newstore;
		fNAlloc = newsize;
	}
}

template< class T, int NumExtAlloc >
void TPZManVector< T, NumExtAlloc >::Shrink()
{
	// Philippe : Porque NumExtAlloc <= fNAlloc????
	//    if(this->fNElements <= NumExtAlloc && NumExtAlloc <= fNAlloc) {
	if (this->fNElements <= NumExtAlloc)
	{
		if (this->fStore != fExtAlloc)
		{
			for (int i = 0; i < this->fNElements; i++)
				fExtAlloc[i] = this->fStore[i];
			
			if (this->fStore)
				delete [] this->fStore;
			
			this->fStore  = fExtAlloc;
			fNAlloc = NumExtAlloc;
		}
	}
	else if (fNAlloc != this->fNElements)
	{  // then  fExtAlloc != this->fStore  because  NumExtAlloc != fNAlloc
		// Philippe : Memoria alocada externamente nao pode ser deletada
		//          if(fExtAlloc) delete[] fExtAlloc;
		T *newstore = 0;
		
		if (this->fNElements)
			newstore = new T[this->fNElements];			
		
		for (int i = 0; i < this->fNElements; i++)
			newstore[i] = this->fStore[i];			
		
		if (this->fStore)
			delete[]this->fStore;
		
		this->fStore = newstore;
		fNAlloc = this->fNElements;
		
		// Philippe Isto eh um absurdo
		//          fExtAlloc = this->fStore;
		//          NumExtAlloc = fNAlloc;
	}
}

template< class T, int NumExtAlloc >
void TPZManVector< T, NumExtAlloc >::Resize(const int newsize, const T& object)
{
#ifndef NODEBUG
	if (newsize < 0)
	{
		PZError << "TManVec::Resize. Bad parameter newsize." << std::endl;
		PZError.flush ();
		return;
	}
	if (newsize == this->fNElements)
		return;
#endif
	
	if (newsize <= fNAlloc)
	{
		for (int i = this->fNElements; i < newsize; i++)
			this->fStore[i] = object;
		
		this->fNElements = newsize;
	}
	else if (newsize <= NumExtAlloc)
	{  // that is, fExtAlloc != this->fStore. Moreover : this->fNElements <= fNAlloc
		// <= NumExtAlloc
		
		int i;
		
		for (i = 0; i < this->fNElements; i++)
		{
			fExtAlloc[i] = this->fStore[i];
		}
		
		for (; i < newsize; i++)
			fExtAlloc[i] = object;
		
		if(this->fStore != fExtAlloc) delete [] this->fStore;
		
		this->fStore = fExtAlloc;
		this->fNElements = newsize;
		fNAlloc = NumExtAlloc;
	}
	else
	{  // the size is larger than the allocated memory, then this->fNElements
		// is always lower than newsize because fNElemets <=fNAllocs
		int i, realsize = ExpandSize (newsize);
		
		T* newstore = new T[realsize];
		
		for (i = 0; i < this->fNElements; i++)
		{
			newstore[i] = this->fStore[i];
		}
		
		for (; i < newsize; i++)
			newstore[i] = object;
		
		if (this->fStore != fExtAlloc)
			delete[]this->fStore;
		
		this->fStore = newstore;
		this->fNElements = newsize;
		fNAlloc = realsize;
	}
}

//
template< class T, int NumExtAlloc >
void TPZManVector< T, NumExtAlloc >::Resize(const int newsize)
{
#ifndef NODEBUG
	if (newsize < 0)
	{
		PZError << "TManVec::Resize. Bad parameter newsize." << std::endl;
		PZError.flush ();
		return;
	}
#endif
	
	if (newsize == this->fNElements)
		return;
	
	if (newsize <= fNAlloc)
	{
		this->fNElements = newsize;
	}
	else if (newsize <= NumExtAlloc)
	{  // that is, fExtAlloc != this->fStore
		if(this->fStore != fExtAlloc) {
			for (int i = 0; i < this->fNElements; i++)
			{
				fExtAlloc[i] = this->fStore[i];
			}
			delete [] this->fStore;
			this->fStore     = fExtAlloc;
		}
		
		this->fNElements = newsize;
		fNAlloc    = NumExtAlloc;
	}
	else
	{  // the size is larger than the allocated memory, then this->fNElements
		// is always lower than newsize because fNElemets <=fNAllocs
		
		int realsize = ExpandSize (newsize);
		T *newstore = new T[realsize];
		
		for (int i = 0; i < this->fNElements; i++)
		{
			newstore[i] = this->fStore[i];
		}
		
		if( this->fStore != fExtAlloc )
			delete [] this->fStore;
		
		this->fStore = newstore;
		this->fNElements = newsize;
		fNAlloc = realsize;
	}
}

template< class T, int NumExtAlloc >
int TPZManVector < T, NumExtAlloc >::ExpandSize (const int proposed) const
{
	return ( proposed > fNAlloc * 1.2 ? proposed : 
			static_cast< int >( fNAlloc * 1.2 ) );
}

#endif

//--| PZ |----------------------------------------------------------------------
