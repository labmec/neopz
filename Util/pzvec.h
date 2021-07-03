/**
 * @file pzvec.h
 * @brief Templated vector implementation.
 */

#ifndef TPZVEC_H
#define TPZVEC_H

#include "pzreal.h"
#include "pzerror.h"
#include <ostream>
#include <cstring>//memcpy
#include <iomanip>//std::setw

/// Overloading the operator <<
inline std::ostream &operator<<(std::ostream &out, const std::pair<int,int> &element)
{
	out << element.first << "|" << element.second;
	return out;
}

/**
 * @ingroup util
 * @brief This class implements a simple vector storage scheme for a templated class T. \ref util "Utility"
 */
/**
 * The copy constructor and operator= requires the
 * operator= to be implemented on the class T.
 */
template< class T >
class TPZVec {
public:
	/** @brief Creates a vector with size 0. */
	TPZVec();
	
	/**
	 * @brief Creates a vector of a given size.
	 * @param size Size of the new vector.
	 */
	/** It will call the empty constructor on all objects of type T created. */
	TPZVec(const int64_t size);
	
	/**
	 * @brief Creates a vector of a given size.
	 * @param size Size of the new vector.
	 * @param copy Model object to initialize the other objects.
	 */
	/**
	 * It will call the empty constructor on all objects of type T
	 * created copies the object copy to all elements.
	 */
	TPZVec(const int64_t size, const T& copy);   

	/**
	 * @brief Creates a vector from a initializer list
	 * @param list : initializer list, usually a list of elements in curly brackets
	 */
	TPZVec(const std::initializer_list<T>& list);
    
	/**
	 * @brief Creates a vector with copy constructor. \n
	 * will call the empty constructor on all objects of type T created
	 * @param copy : original vector
	 */
	TPZVec(const TPZVec<T> &copy);

    /**
     * @brief Creates a vector using the move constructor.
     * No memory is allocated.
     */
    TPZVec(TPZVec<T> &&rval);
    
	/** @brief destructor, will delete the storage allocated */
	virtual ~TPZVec();
	/**
	 * @brief will copy the vector into the current vector.
	 * @param copy vector which will be copied
	 * @return reference to the current object
	 */
	/** Will first delete the allocated storage before allocating storage for the copy */
	TPZVec<T> &operator=(const TPZVec<T> &copy);

    /**
     * @brief will move the vector into the current vector
     */
    TPZVec<T> &operator=(TPZVec<T> &&rval);
  
	void ReallocForNuma() {
        if (fNElements == 0 || fStore == NULL)
            return;
        T* newStore = new T[fNElements];
        memcpy((void*)newStore,(void*)fStore,fNElements*sizeof(T));
        delete [] fStore;
        fStore = newStore;
	}

	/**
	 * @brief create a new vector and stores in the current vector
	 * @param list: initializer list, usually a list in curly brackets
	 * @return reference to the current object
	 */
	/** Will first delete the allocated storage before allocating storage for the copy */
	TPZVec<T> &operator=(const std::initializer_list<T> &list);

	/**
	 * @brief Operator attribution. Fills the vector with a value of type T.
	 * @param a Element to fill the vector with.
	 * @return Reference to the current object.
	 */
	TPZVec<T>& operator=(const T& a);
	
	/**
	 * @brief Access operator, will perform bounds checking unless the variable PZNODEBUG is defined
	 * @param index element in the vector which will be acessed
	 * @return a reference to the element specified by index\
	 */
#if defined(WIN32) && !defined(_WIN64)
	T& operator[]( const int index ) const
	{
#ifdef PZDEBUG
		if( index < 0 || index >= fNElements )
		{
            PZError << __PRETTY_FUNCTION__ << " acessing element out of range.";
            PZError << "|" << std::endl;
            PZError << "+-> NElements = " << NElements() << std::endl;
            PZError << "|" << std::endl;
            PZError << "+-> Index = " << index << std::endl;
            DebugStop();
            exit( -1 );
		}
#endif
		return fStore[ index ];
	}
	T& operator[]( const unsigned int index ) const
	{
#ifdef PZDEBUG
		if( index >= (unsigned int)fNElements )
		{
            PZError << __PRETTY_FUNCTION__ << " acessing element out of range.";
            PZError << "|" << std::endl;
            PZError << "+-> NElements = " << NElements() << std::endl;
            PZError << "|" << std::endl;
            PZError << "+-> Index = " << index << std::endl;
            DebugStop();
            exit( -1 );
		}
#endif
		return fStore[ index ];
	}
#endif //WIN32
	T& operator[]( const int64_t index ) const
	{
#ifdef PZDEBUG
		if( index < 0 || index >= fNElements )
		{
            PZError << __PRETTY_FUNCTION__ << " acessing element out of range.";
            PZError << "|" << std::endl;
            PZError << "+-> NElements = " << NElements() << std::endl;
            PZError << "|" << std::endl;
            PZError << "+-> Index = " << index << std::endl;
            DebugStop();
            exit( -1 );
		}
#endif
		return fStore[ index ];
	}
	
	/** @brief Extraction operator. */
	template <class T2>
	friend std::ostream& operator<<( std::ostream& Out, const TPZVec< T2 >& v );
	
	/** @brief Casting operator. Returns The fStore pointer. */
	//operator T*() const { return fStore; }
  /** @{ */
  /** @brief Returns a pointer to the first element */
  T *begin() const;
  T *&begin();
  /** @{ */
  /** @brief Returns a pointer to the last+1 element */
  T *end() const;
	/**
	 * @brief Will fill the elements of the vector with a copy object.
	 * @param copy object which will be copied
	 * @param from first index which will be overwritten
	 * @param numelem number of elements which will be overwritten
	 */
	void Fill(const T& copy, const int64_t from=0, const int64_t numelem=-1);
	
	/**
	 * @brief Returns the number of elements of the vector
	 * @return number of elements used by the vector
	 */
	inline int64_t NElements() const { return fNElements; }

    /** @brief Returns number of elements dynamically for this object. */
    inline int64_t NAlloc() const { return fNAlloc; }
	/**
	 * @brief Returns the number of elements of the vector
	 * @return number of elements used by the vector
	 */
	inline int64_t size() const { return fNElements; }
	
	/**
	 * @brief Resizes the vector object reallocating the necessary storage,
	 * copying the existing objects to the new storage.
	 * @param newsize size of the vector
	 * @param object object used to initialize the new members
	 */
	virtual void Resize(const int64_t newsize,const T& object);
	
	/**
	 * @brief Resizes the vector object reallocating the necessary storage,
	 * copying the existing objects to the new storage. \n The new
	 * members are not initialized.
	 * @param newsize size of the vector
	 */
	virtual void Resize(const int64_t newsize);
	virtual void resize(const int64_t newsize)
	{
		Resize(newsize);
	}
	
	/**
	 * @brief Prints the structural information of the vector object to the
	 * output stream. \n This method will not print the objects
	 * themselves!
	 */
	void Print(std::ostream &out = std::cout);
    
    /** @brief Empty the vector, make its size zero */
    virtual void clear();
	
protected:
	/** @brief Allocated storage for the vector object */
	T* fStore;
	
	/** @brief Number of elements of the vector object */
	int64_t fNElements;
    
    /** @brief Number of elements allocated for this object. */
    int64_t fNAlloc;
    /*@orlandini note: for TPZVec<T> fNAlloc is always equal to
     * fNElements. This wont hold true at all times for derived types.
     */
};

template< class T >
inline TPZVec<T>::TPZVec() : fStore( nullptr ), fNElements( 0 ), fNAlloc(0)
{
	// NOTHING TO DO HERE!
}

template< class T >
TPZVec<T>::TPZVec( const int64_t size ) : fStore( nullptr )
{
#ifndef PZNODEBUG
	if( size < 0 )
	{
		PZError << __PRETTY_FUNCTION__;
        PZError << "Bad parameter size, then size = 0.\n";
		PZError.flush();
		fNElements = 0;
        fNAlloc = 0;
		return;
	}
#endif
	
	// If a positive value was requested, allocate it.
	if( size > 0 )
	{
		fStore = new T[ size ];
	}
	
	// Note that even 0 sized vectors are allowed.
	fNElements = size;
    fNAlloc = size;
    
}

template< class T >
TPZVec<T>::TPZVec( const int64_t size, const T& copy ) : fStore( nullptr )
{
#ifndef PZNODEBUG
	if( size < 0 )
	{
		PZError << __PRETTY_FUNCTION__;
        PZError << "Bad parameter size, then size = 0.\n";
		PZError.flush();
		fNElements = 0;
        fNAlloc = 0;
		return;
	}
#endif
	
	if( size )
	{
		fStore = new T[size];
	}
	
	fNElements = size;
    fNAlloc = size;
	
	for( int64_t i = 0; i < size; i++ )
	{
		fStore[i] = copy;
	}
}


template< class T >
TPZVec<T>::TPZVec(const std::initializer_list<T>& list) : fStore(nullptr) {
    const auto size = list.size();
	if (size > 0)
		fStore = new T[size];

	auto it = list.begin();
	auto it_end = list.end();
	T* aux = fStore;
	for (; it != it_end; it++, aux++)
		*aux = *it;

	fNElements = size;
    fNAlloc = size;
}

template< class T >
TPZVec<T>::TPZVec(const TPZVec<T> &copy){
	fStore = nullptr;
	
	if( copy.fNElements > 0 )
		fStore = new T[copy.fNElements];
	else
		fStore = nullptr;
	
	for(int64_t i=0; i<copy.fNElements; i++)
		fStore[i]=copy.fStore[i];
	
	fNElements = copy.fNElements;
    fNAlloc = copy.fNAlloc;
}

template< class T >
TPZVec<T>::TPZVec(TPZVec<T> &&rval) : fNElements(rval.fNElements){
    if(rval.fNAlloc){
        fStore = rval.fStore;
        fNAlloc = rval.fNAlloc;
        rval.fStore = nullptr;
        rval.fNAlloc = 0;
    }else{//rval did not allocate memory. let us just copy
        fStore = new T[fNElements];
        for(int64_t i=0; i<fNElements; i++)
            fStore[i]=rval.fStore[i];
    }
    
}

template<class T>
inline TPZVec<T>::~TPZVec() {
	if( fStore )
	{
		delete [] fStore;
	}
}

template< class T >
TPZVec<T> &TPZVec<T>::operator=(const TPZVec<T> &copy){
    if (this != &copy) {

      Resize(copy.NElements());

      for (int64_t i = 0; i < copy.fNElements; i++)
          fStore[i] = copy.fStore[i];
    }
    return *this;
}

template< class T >
TPZVec<T> &TPZVec<T>::operator=(TPZVec<T> &&rval){
	if(this != &rval){
        fNElements = rval.fNElements;
        if(fStore) delete [] fStore;
        if(rval.fNAlloc){
            fStore = rval.fStore;
        }else{
            fStore = new T[fNElements];
            for (int64_t i = 0; i < fNElements; i++)
                fStore[i] = rval.fStore[i];
        }
        fNAlloc = fNElements;//perhaps rval.fNalloc was 0
        rval.fStore = nullptr;
        rval.fNElements = 0;
        rval.fNAlloc = 0;
	}
	return *this;
}

template< class T >
TPZVec<T> &TPZVec<T>::operator=(const std::initializer_list<T> &list) {
	Resize(list.size());

	auto it = list.begin();
	auto it_end = list.end();
	auto *aux = fStore;

	for (; it != it_end; it++, aux++)
		*aux = *it;
	return *this;
}

// OPENED QUESTION: what to do with 0 size vectors??? Cantao (2002.01.09)
template< class T >
TPZVec<T>& TPZVec<T>::operator=( const T& a )
{
	T* end = fStore + fNElements;
	
	for( T* walk = fStore; walk < end; *walk++ = a );
	
	return *this;
}


template< class T >
void TPZVec<T>::Resize(const int64_t newsize,const T& object) {
#ifndef PZNODEBUG
	if(newsize<0) {
		PZError << "TPZVec::Resize. Bad parameter newsize." << std::endl;
		PZError.flush();
		return;
	}
#endif
	if(newsize == fNElements) return;
	T* newstore;
	if(newsize) newstore = new T[newsize];
	else newstore = 0;
	int64_t large = (fNElements < newsize) ? fNElements : newsize;
	int64_t i;
	for(i=0L; i<large; i++) {
		newstore[i] = fStore[i];
	}
	for(;i<newsize;i++) {   // then only to case : large=fNElement < newsize
		newstore[i] = object;
	}
    if(fStore) {
        delete[] fStore;
        fStore = nullptr;
    }
	fStore = newstore;
	fNElements = newsize;//cedric 20/11/99 e 29/04/00
    fNAlloc = newsize;
}

#include <limits>

template< class T >
void TPZVec<T>::Resize(const int64_t newsize) {
#ifndef PZNODEBUG
    //	int64_t nlongsize = 2147483647;
	if(newsize<0) {
		PZError << "TPZVec::Resize. Bad parameter newsize: " << newsize <<  std::endl;
		PZError.flush();
	}
#ifdef WIN32
	// Parece que o limite no windows é
	int sz = sizeof(T);
	int64_t nlongsize = 1704792168;
	if((newsize+1) > (1./sz)*nlongsize) {
		PZError << "TPZVec::Resize. Bad parameter newsize: " << newsize <<  std::endl;
		std::cout << "TPZVec::Resize. Bad parameter newsize: " << newsize <<  std::endl;
		PZError.flush();
        //		DebugStop();
        //		return;
	}
#endif
#endif
	
	if(newsize == fNElements) return;
	if (newsize == 0) {
		fNElements = 0;
        fNAlloc = 0;
		delete[] fStore;
		fStore = nullptr;
		return;
	}
	T *newstore = new T[newsize];
	int64_t large = (fNElements < newsize) ? fNElements : newsize;
	int64_t i;
	for(i=0L; i<large; i++) {
		newstore[i] = fStore[i];
	}
	if(fStore) delete[] fStore;
	fStore = newstore;
	fNElements = newsize;
    fNAlloc = newsize;
}

template<class T>
void TPZVec<T>::clear()
{
    this->Resize(0);
}

template<class T>
T *TPZVec<T>::begin() const {
    return fStore;
}

template<class T>
T *&TPZVec<T>::begin() {
    return fStore;
}
template<class T>
T *TPZVec<T>::end() const {
    return fStore+fNElements;
}

template< class T >
void TPZVec<T>::Fill(const T& copy, const int64_t from, const int64_t numelem){
#ifndef PZNODEBUG
	if(numelem<0 && numelem != -1) {
		PZError << "TPZVec::Fill" << std::endl
		<< "It's negative parameter numelem, then numelem = "
		<< fNElements << std::endl;
		
		PZError.flush();
	}
#endif
	
	int64_t first = (from < 0) ? 0 : from;
	int64_t nel = numelem;
	first = (first > fNElements) ? fNElements : first;
	if (nel < 0) nel = fNElements;
	int64_t last = (from+nel > fNElements) ? fNElements : from+nel;
	
	for(int64_t i=first; i<last; i++)
		fStore[i] = copy;
}

template< class T >
inline void TPZVec<T>::Print(std::ostream &out)
{
	out << std::endl << "Number of elements = " << fNElements;
}

template <class T>
std::ostream& operator<<( std::ostream& Out, const TPZVec< T >& v )
{
	std::streamsize width = Out.width();
	
	const char* sep = ( width == 0 ? ", " : "" );
	
	int64_t size = v.NElements();
	
	if(size) Out << std::setw(width) << v.fStore[0];
	
	for( int64_t ii = 1; ii < size; ii++ )
	{
	    Out << std::setw( width ) << sep << v.fStore[ ii ];
	}
	
	return Out;
}

inline std::ostream& operator<<( std::ostream& Out, const TPZVec< std::pair<double,double> >& v )
{
    
	Out << "{";
    
	int64_t size = v.NElements();
	
    if(size > 0)
    {
        for( int64_t ii = 0; ii < size; ii++ )
        {
            Out << "{" << v[ii].first << ',' << v[ii].second << "}";
            if( ii < size-1) Out << ",";
        }
    }
	
    Out << "}";
	return Out;
}

extern template class TPZVec<float>;
extern template class TPZVec<float * >;
extern template class TPZVec<double>;
extern template class TPZVec<double * >;
extern template class TPZVec<long double>;
extern template class TPZVec<long double * >;
extern template class TPZVec<int>;
extern template class TPZVec<int64_t>;
extern template class TPZVec<int64_t *>;
extern template class TPZVec<int *>;
extern template class TPZVec<char *>;
extern template class TPZVec<void *>;
extern template class TPZVec<char>;

class TPZGeoEl;
class TPZGeoNode;
struct TPZGeoNodeBC;
struct TPZGeoElBC;
extern template class TPZVec<TPZGeoEl *>;
extern template class TPZVec<TPZGeoNode *>;
extern template class TPZVec<TPZGeoNodeBC *>;
extern template class TPZVec<TPZGeoElBC *>;


class TPZCompEl;
class TPZConnect;
struct TPZConnectBC;
extern template class TPZVec<TPZCompEl *>;
extern template class TPZVec<TPZConnect *>;
extern template class TPZVec<TPZConnectBC *>;

class TPZMaterial;
class TPZBndCond;
extern template class TPZVec<TPZMaterial *>;
extern template class TPZVec<TPZBndCond *>;

struct TPZElementMatrix;
extern template class TPZVec<TPZElementMatrix *>;

class TPZGraphEl;
class TPZGraphNode;
class TPZCompMesh;
extern template class TPZVec<TPZGraphEl *>;
extern template class TPZVec<TPZGraphNode *>;
extern template class TPZVec<TPZCompMesh *>;
#endif
