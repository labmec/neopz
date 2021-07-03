/**
 * @file
 * @brief Free store vector implementation.
 */
#ifndef PZMANVECTOR_H
#define PZMANVECTOR_H

#include "pzvec.h"
#include "pzerror.h"

/** @brief To allocate vector by default */
/** @ingroup util */
const int DEFAULTVEC_ALLOC = 10;

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
class TPZManVector : public TPZVec< T > {
public:
    /**
     * @brief Creates a vector of a given size.
     * @param size Size of the new vector.
     */
    /**
     * It will call the empty constructor on all objects of type T.
     */
    TPZManVector(const int64_t size = 0);

    /**
     * @brief Creates a vector of a given size, filling it.
     * @param size Size of the new vector.
     * @param copy Model object to initialize the other objects.
     */
    /**
     * It will call the empty constructor on all objects of type T
     * created. Copies the object copy to all elements.
     */
    TPZManVector(const int64_t size, const T& copy);


    TPZManVector(const TPZVec<T> & copy);

	/**
	 * @brief Creates a vector from a initializer list
	 * @param list: the initializer list, usually enclosed in curly brackets
	 */
	TPZManVector(const std::initializer_list<T>& list);

    /**
     * @brief Copy constructor.
     */
    TPZManVector(const TPZManVector< T, NumExtAlloc >& copy);
    /**
     * @brief Move constructor.
     */
    TPZManVector(TPZManVector< T, NumExtAlloc >&& copy);
    /** @brief Destructor. */
    /** Deletes the storage allocated. */
    virtual ~TPZManVector();
    
    /**
     * @brief Copy assignment operator.
     * @param copy Vector which will be copied.
     * @return Reference to the current object.
     */
    /**
     * It first deletes the allocated storage before allocating
     * storage for the copy only if necessary (when there is no
     * preallocated storage or when the current storage cannot hold
     * the copied vector
     */
    TPZManVector< T, NumExtAlloc >& operator=(const TPZManVector< T, NumExtAlloc >& copy);
    /**
     * @brief Move assignment operator.
     */
    TPZManVector< T, NumExtAlloc >& operator=(TPZManVector< T, NumExtAlloc >&& copy);
	/**
	 * @brief initializer list assignment operator
	 * @param list: list which will be assigned, usually wrapped in curly brackets
	 * @return Reference to the current object
	 */
	TPZManVector< T, NumExtAlloc >& operator=(const std::initializer_list<T>& list);

    /**
     * @brief Expands the allocated storage to fit the newsize parameter.
     * @param newsize Storage size which is requested.
     */
    /** Does nothing if the externally provided storage is larger than newsize. Does not change the size of the vector. */
    void Expand(const int64_t newsize);

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
    virtual void Resize(const int64_t newsize, const T& object) override;

    /**
     * @brief Resizes the vector object reallocating the storage if necessary.
     * @see TPZVec<T>::Resize(const int);
     * @param newsize Size of the vector.
     */
    /**
     * It copies the existing objects to the new storage. The new
     * members are not initialized.
     */
    virtual void Resize(const int64_t newsize) override;

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
    int64_t ExpandSize(const int64_t proposed) const;

    /** @brief Pointer to the externally allocated space. */
    T fExtAlloc[ NumExtAlloc ];
};

//--| IMPLEMENTATION |----------------------------------------------------------

template< class T, int NumExtAlloc >
TPZManVector< T, NumExtAlloc >::TPZManVector(const int64_t size) :
TPZVec<T>(0) // There is always some static allocation.
{
    /* If the size requested fits inside the size already provided
     * statically, provides access to that space, by setting some
     * TPZVec data.
     */
    if (size <= NumExtAlloc) {
        // Needed to make TPZVec::operator[] work properly.
        this->fStore = fExtAlloc;
        this->fNElements = size;
        // No memory was allocated by the constructor.
        this->fNAlloc = 0;
    } else // The size requested is bigger than the size already provided.
    {
        // Executes the allocation that would be done by TPZVec<T>(size).
        this->fStore = new T[ size ];
        this->fNElements = size;
        this->fNAlloc = size;
    }
}

template< class T, int NumExtAlloc >
TPZManVector< T, NumExtAlloc >::TPZManVector(const int64_t size, const T& copy) :
TPZVec<T>(0)
//fNExtAlloc( NumExtAlloc ) // There is always some static allocation.
{
    /* If the size requested fits inside the size already provided
     * statically, provides access to that space, by setting some
     * TPZVec data.
     */
    if (size <= NumExtAlloc) {
        // Needed to make TPZVec::operator[] work properly.
        this->fStore = fExtAlloc;
        this->fNElements = size;
        // No memory was allocated by the constructor.
        this->fNAlloc = 0;
    } else // The size requested is bigger than the size already provided.
    {
        // Executes the allocation that would be done by TPZVec<T>(size).
        this->fStore = new T[ size ];
        this->fNElements = size;
        this->fNAlloc = size;
    }

    for (int64_t i = 0; i < size; i++) {
        this->fStore[i] = copy;
    }
}


template< class T, int NumExtAlloc>
inline TPZManVector< T, NumExtAlloc >::TPZManVector(
        const TPZVec<T> & copy) {
    const int64_t size = copy.NElements();

    /* If the size requested fits inside the size already provided
     * statically, provides access to that space, by setting some
     * TPZVec data.
     */
    if (size <= (int64_t) (sizeof (fExtAlloc) / sizeof (T))) {
        // Needed to make TPZVec::operator[] work properly.
        this->fStore = fExtAlloc;
        this->fNElements = size;
        // No memory was allocated by the constructor.
        this->fNAlloc = 0;
    } else // The size requested is bigger than the size already provided.
    {
        // Executes the allocation that would be done by TPZVec<T>(size).
        this->fStore = new T[ size ];
        this->fNElements = size;
        this->fNAlloc = size;
    }

    for (int64_t i = 0; i < size; i++) {
        this->fStore[i] = copy[i];
    }
}

template< class T, int NumExtAlloc>
inline TPZManVector< T, NumExtAlloc >::TPZManVector(const std::initializer_list<T>& list)
{
	int size = list.size();
//	std::cout << "init pzmanvec" << std::endl;
	if (size <= (int64_t)(sizeof(fExtAlloc) / sizeof(T))) {
		// Needed to make TPZVec::operator[] work properly.
		this->fStore = fExtAlloc;
		this->fNElements = size;
		// No memory was allocated by the constructor.
		this->fNAlloc = 0;
	}
	else // The size requested is bigger than the size already provided.
	{
		// Executes the allocation that would be done by TPZVec<T>(size).
		this->fStore = new T[size];
		this->fNElements = size;
		this->fNAlloc = size;
	}

	auto it_end = list.end();
	T* aux = this->fStore;
	for (auto it = list.begin(); it != it_end; it++, aux++)
		*aux = *it;
}

template< class T, int NumExtAlloc>
inline TPZManVector< T, NumExtAlloc >::TPZManVector(const TPZManVector< T, NumExtAlloc >& copy) {
    const int64_t size = copy.NElements();

    /* If the size requested fits inside the size already provided
     * statically, provides access to that space, by setting some
     * TPZVec data.
     */
    if (size <= NumExtAlloc) {
        // Needed to make TPZVec::operator[] work properly.
        this->fStore = fExtAlloc;
        this->fNElements = size;
        // No memory was allocated by the constructor.
        this->fNAlloc = 0;
    } else // The size requested is bigger than the size already provided.
    {
        // Executes the allocation that would be done by TPZVec<T>(size).
        this->fStore = new T[ size ];
        this->fNElements = size;
        this->fNAlloc = size;
    }

    for (int64_t i = 0; i < size; i++) {
        this->fStore[i] = copy.fStore[i];
    }
}

template< class T, int NumExtAlloc>
inline TPZManVector< T, NumExtAlloc >::TPZManVector(TPZManVector< T, NumExtAlloc >&& rval) {

    /* If the size requested fits inside the size already provided
     * statically, provides access to that space, by setting some
     * TPZVec data.
     */
    if (auto size = rval.NElements();size <= NumExtAlloc){
        //we need to copy, unfortunately
        int i = 0;
        for(; i < size; i++) {fExtAlloc[i] = rval.fExtAlloc[i];}
        for(; i < NumExtAlloc; i++){fExtAlloc[i] = T();}
        this->fStore = fExtAlloc;
        this->fNElements = size;
        // No memory was allocated by the constructor.
        this->fNAlloc = 0;
    } else {// The size requested is bigger than the size already provided.
        this->fStore = rval.fStore;
        this->fNElements = size;
        this->fNAlloc = size;
    }
    rval.fStore = nullptr;
    rval.fNAlloc = 0;
    rval.fNElements = 0;
}

template< class T, int NumExtAlloc >
TPZManVector< T, NumExtAlloc >::~TPZManVector() {
    if (this->fStore == fExtAlloc) {
        this->fStore = nullptr;
    }
    this->fNAlloc = 0;
}

template< class T, int NumExtAlloc >
TPZManVector< T, NumExtAlloc >& TPZManVector< T, NumExtAlloc >::operator=(const TPZManVector< T, NumExtAlloc >& copy) {
    // Checking auto assignment.
    if (this != &copy) {
      const int64_t nel = copy.NElements();

      if (nel > this->fNAlloc && this->fStore && this->fStore != fExtAlloc) {
        delete[] this->fStore;
        this->fStore = 0;
        this->fNAlloc = 0;
      }

      if (nel <= NumExtAlloc) {
        if (this->fStore != fExtAlloc) {
          delete[] this->fStore;
        }
        this->fNAlloc = 0;
        this->fStore = fExtAlloc;
        this->fNElements = nel;
      } else if (this->fNAlloc >= nel) {
        this->fNElements = nel;
      } else {
        this->fStore = new T[nel];
        this->fNAlloc = nel;
        this->fNElements = nel;
      }

      for (int64_t i = 0; i < nel; i++) {
        this->fStore[i] = copy.fStore[i];
      }
    }
    return *this;
}

template< class T, int NumExtAlloc >
TPZManVector< T, NumExtAlloc >& TPZManVector< T, NumExtAlloc >::operator=(TPZManVector< T, NumExtAlloc >&& rval) {
    // Checking auto assignment.
    if (this != &rval) {
      const int64_t nel = rval.NElements();
      //let us dispose of previously allocated memory
      if (this->fStore && this->fStore != fExtAlloc) {
        delete[] this->fStore;
        this->fStore = nullptr;
        this->fNAlloc = 0;
      }
      if (nel <= NumExtAlloc) {
        //we need to copy, unfortunately
        int i = 0;
        for(; i < nel; i++) {fExtAlloc[i] = rval.fExtAlloc[i];}
        for(; i < NumExtAlloc; i++){fExtAlloc[i] = T();}
        this->fStore = fExtAlloc;
        this->fNAlloc = 0;
        this->fNElements = nel;
      } else {
        this->fStore = rval.fStore;
        this->fNAlloc = nel;
        this->fNElements = nel;
      }
    }
    rval.fStore = nullptr;
    rval.fNAlloc = 0;
    rval.fNElements = 0;
    return *this;
}

template< class T, int NumExtAlloc >
TPZManVector< T, NumExtAlloc >& TPZManVector< T, NumExtAlloc >::operator=(const std::initializer_list<T>& list) {
	size_t size = list.size();

	if (size > this->fNAlloc && this->fStore && this->fStore != this->fExtAlloc) {
		delete[] this->fStore;
		this->fStore = 0;
		this->fNAlloc = 0;
	}

	if (size <= NumExtAlloc) {
		if (this->fStore != fExtAlloc) {
			delete[]this->fStore;
		}
		this->fNAlloc = 0;
		this->fStore = fExtAlloc;
		this->fNElements = size;
	}
	else if (this->fNAlloc >= size) {
		this->fNElements = size;
	}
	else {
		this->fStore = new T[size];
		this->fNAlloc = size;
		this->fNElements = size;
	}

	auto it_end = list.end();
	T* aux = this->fStore;
	for (auto it = list.begin(); it != it_end; it++, aux++)
		*aux = *it;

	return *this;
}

template< class T, int NumExtAlloc >
void TPZManVector< T, NumExtAlloc >::Expand(const int64_t newsize) {
    // If newsize is negative then return.
    if (newsize <= this->fNAlloc || newsize <= NumExtAlloc) {
        return;
    } else {// the size is larger than the allocated memory
        T* newstore = new T[ newsize ];

        for (int64_t i = 0L; i < this->fNElements; i++) {
            newstore[i] = this->fStore[i];
        }

        if (this->fStore != fExtAlloc) {
            delete [] this->fStore;
        }

        this->fStore = newstore;
        this->fNAlloc = newsize;
    }
}

template< class T, int NumExtAlloc >
void TPZManVector< T, NumExtAlloc >::Shrink() {
    // Philippe : Porque NumExtAlloc <= this->fNAlloc????
    //    if(this->fNElements <= NumExtAlloc && NumExtAlloc <= this->fNAlloc) {
    if (this->fNElements <= NumExtAlloc) {
        if (this->fStore != fExtAlloc) {
            for (int64_t i = 0; i < this->fNElements; i++)
                fExtAlloc[i] = this->fStore[i];

            if (this->fStore)
                delete [] this->fStore;

            this->fStore = fExtAlloc;
            this->fNAlloc = NumExtAlloc;
        }
    } else if (this->fNAlloc != this->fNElements) { // then  fExtAlloc != this->fStore  because  NumExtAlloc != this->fNAlloc
        // Philippe : Memoria alocada externamente nao pode ser deletada
        //          if(fExtAlloc) delete[] fExtAlloc;
        T *newstore = 0;

        if (this->fNElements)
            newstore = new T[this->fNElements];

        for (int64_t i = 0; i < this->fNElements; i++)
            newstore[i] = this->fStore[i];

        if (this->fStore)
            delete[]this->fStore;

        this->fStore = newstore;
        this->fNAlloc = this->fNElements;

        // Philippe Isto eh um absurdo
        //          fExtAlloc = this->fStore;
        //          NumExtAlloc = this->fNAlloc;
    }
}

template< class T, int NumExtAlloc >
void TPZManVector< T, NumExtAlloc >::Resize(const int64_t newsize, const T& object) {
#ifndef PZNODEBUG
    if (newsize < 0) {
        PZError << "TManVec::Resize. Bad parameter newsize." << std::endl;
        PZError.flush();
        return;
    }
    if (newsize == this->fNElements)
        return;
#endif

    if (newsize <= this->fNAlloc) {
        for (int64_t i = this->fNElements; i < newsize; i++)
            this->fStore[i] = object;

        this->fNElements = newsize;
    } else if (newsize <= NumExtAlloc) { // that is, fExtAlloc != this->fStore. Moreover : this->fNElements <= this->fNAlloc
        // <= NumExtAlloc

        int64_t i;

        for (i = 0L; i < this->fNElements; i++) {
            fExtAlloc[i] = this->fStore[i];
        }

        for (; i < newsize; i++)
            fExtAlloc[i] = object;

        if (this->fStore != fExtAlloc) delete [] this->fStore;

        this->fStore = fExtAlloc;
        this->fNElements = newsize;
        this->fNAlloc = NumExtAlloc;
    } else { // the size is larger than the allocated memory, then this->fNElements
        // is always lower than newsize because fNElemets <=this->fNAllocs
        int64_t i, realsize = ExpandSize(newsize);

        T* newstore = new T[realsize];

        for (i = 0L; i < this->fNElements; i++) {
            newstore[i] = this->fStore[i];
        }

        for (; i < newsize; i++)
            newstore[i] = object;

        if (this->fStore != fExtAlloc)
            delete[]this->fStore;

        this->fStore = newstore;
        this->fNElements = newsize;
        this->fNAlloc = realsize;
    }
}

template< class T, int NumExtAlloc >
void TPZManVector< T, NumExtAlloc >::Resize(const int64_t newsize) {
#ifndef PZNODEBUG
    if (newsize < 0) {
        PZError << "TManVec::Resize. Bad parameter newsize." << std::endl;
        PZError.flush();
        return;
    }
#endif

    if (newsize == this->fNElements)
        return;

    if (newsize <= this->fNAlloc) {
        this->fNElements = newsize;
    } else if (newsize <= NumExtAlloc) { // that is, fExtAlloc != this->fStore
        if (this->fStore != fExtAlloc) {
            for (int64_t i = 0L; i < this->fNElements; i++) {
                fExtAlloc[i] = this->fStore[i];
            }
            delete [] this->fStore;
            this->fStore = fExtAlloc;
        }

        this->fNElements = newsize;
        this->fNAlloc = 0;
    } else { // the size is larger than the allocated memory, then this->fNElements
        // is always lower than newsize because fNElemets <=this->fNAllocs

        int64_t realsize = ExpandSize(newsize);
        T *newstore = new T[realsize];

        for (int64_t i = 0L; i < this->fNElements; i++) {
            newstore[i] = this->fStore[i];
        }

        if (this->fStore != fExtAlloc)
            delete [] this->fStore;

        this->fStore = newstore;
        this->fNElements = newsize;
        this->fNAlloc = realsize;
    }
}

template< class T, int NumExtAlloc >
int64_t TPZManVector < T, NumExtAlloc >::ExpandSize(const int64_t proposed) const {
    return ( proposed > this->fNAlloc * 1.2 ? proposed :
            static_cast<int64_t> (this->fNAlloc * 1.2));
}


//this three are often used
extern template class TPZManVector< float,3 >;
extern template class TPZManVector< double,3 >;
extern template class TPZManVector< long double,3 >;

extern template class TPZManVector< int >;
extern template class TPZManVector< int64_t >;
extern template class TPZManVector< int *>;
extern template class TPZManVector< char *>;
extern template class TPZManVector< float >;
extern template class TPZManVector< float * >;
extern template class TPZManVector< double >;
extern template class TPZManVector< double * >;
extern template class TPZManVector< long double >;
extern template class TPZManVector< long double * >;
#endif
