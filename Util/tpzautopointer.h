/**
 * @file
 * @brief Contains declaration of the TPZAutoPointer class.
 */

#ifndef TPZAUTOPOINTER_H
#define TPZAUTOPOINTER_H

#include <stdlib.h>
#include <atomic>
#include <mutex>
#include <stdexcept>
// #include <iostream>//useful for debugging

/**
 * \addtogroup util
 * @{
 */

// /**
//  * @brief Increment and Decrement actions are mutexed by this mutex
//  */

namespace pzinternal{
    extern std::recursive_mutex g_ap_mut;
}

/**
 * @brief This class implements a reference counter mechanism to administer a dynamically allocated object. \ref util "Utility"
 * @author Philippe R. B. Devloo
 */
template<class T>
class TPZAutoPointer {
    
	/** @brief Counter struct */
	template<class T2>
    struct TPZReference
    {
        /** @brief Pointer to T2 object */
        T2 *fPointer;
        /** @brief Reference counter*/
        std::atomic_int *fCounter;
        
        TPZReference()
        {
            fPointer = nullptr;
            fCounter = new std::atomic_int{};
            fCounter->store(1);
        }
        
        TPZReference(T2 *pointer)
        {
            fPointer = pointer;
            fCounter = new std::atomic_int{};
            fCounter->store(1);
        }

        TPZReference(const TPZReference &) = delete;        
        TPZReference(TPZReference&&rval) = delete;
        TPZReference& operator=(const TPZReference &copy) = delete;
        TPZReference& operator=(TPZReference &&rval) = delete;
        
        ~TPZReference()
        {
            if(fPointer) {
                delete fPointer;
            }
            fPointer = nullptr;
            
            if(fCounter) {
                delete fCounter;
            }
            fCounter = nullptr;
        }
        

        /** @brief Increment the counter */
        bool Increment()
        {
//            std::mutex *mut = get_ap_mutex((void*) this);
//            std::unique_lock<std::mutex> lck(*mut);//already locks
            fCounter->fetch_add(1);
//            lck.unlock();
            return true;
        }
        /** @brief Decrease the counter. The return value will specify if
            it is needed to delete the reference.*/
        bool Decrease()
        {
            if(!fCounter){
                throw std::logic_error("TPZAutoPointer::Decrease() called without fCounter");
            }
            bool should_delete = false;
            fCounter->fetch_sub(1);
            
            if((*fCounter) == 0) should_delete = true;
            else if((*fCounter) < 0){
                throw std::logic_error(
                    "Invalid value for ref counter of TPZAutoPointer");
            }
            return should_delete;
        }
        
    };
    
	/** @brief The object which contains the pointer and the reference count */
	TPZReference<T> *fRef;
    
public:
	/** @brief Creates an reference counted null pointer */
	TPZAutoPointer() : fRef(new TPZReference<T>())
	{
		
	}
    
	/** @brief The destructor will delete the administered pointer if its reference count is zero */
	~TPZAutoPointer()
	{        
        if (fRef && fRef->Decrease()){
            Release();
        }
	}
    
	/** @brief This method will create an object which will administer the area pointed to by obj */
	TPZAutoPointer(T *obj) : fRef(new TPZReference<T>(obj))
	{

	}
    
	/** @brief Share the pointer of the copy */
	TPZAutoPointer(const TPZAutoPointer<T> &copy) :
        fRef(copy.fRef)
	{
		fRef->Increment();
	}
    /** @brief Acquire the pointer of the rval */
	TPZAutoPointer(TPZAutoPointer<T> &&rval) :
        fRef(rval.fRef)
	{
		rval.fRef = nullptr;
	}
        
	/** @brief Move assignment operator */

	TPZAutoPointer &operator=(TPZAutoPointer<T> &&rval)
    {
        if (fRef && fRef->Decrease()){
            Release();
        }
        fRef = rval.fRef;
        rval.fRef = nullptr;
        return *this;
    }
    
	/** @brief Assignment operator */
	TPZAutoPointer &operator=(const TPZAutoPointer<T> &copy)
	{
		if(copy.fRef == fRef) return *this;		
        else if (fRef && fRef->Decrease()){
            Release();
        }
		fRef = copy.fRef;
        copy.fRef->Increment();
        
		return *this;
	}

    /** @brief Raw ptr assignment operator */

	TPZAutoPointer &operator=(T *rval)
    {
        if (fRef && fRef->Decrease()){
            Release();
        }
        fRef = new TPZReference<T>(rval);
        return *this;
    }
    
	/** @brief User-defined conversion from TPZAutoPointer<T> to &T */
	operator T&()
	{
		return *(fRef->fPointer);
	}
        
	/** @brief Returns the referenced object */
	T& operator *() const
	{
		return *(fRef->fPointer);
	}
        
	// /** @brief Returns the referenced object */
	// T& operator *()
	// {
	// 	return *(fRef->fPointer);
	// }
    
	/** @brief Returns the pointer for referenced object */
	T *operator->() const
	{
		return fRef->fPointer;
	}
	// T *operator->()
	// {
		// return fRef->fPointer;
    // }
    
	inline void ReallocForNuma([[maybe_unused]] int node)
	{
		return;
	}    
    
	/** @brief Returns if pointer was attributed */
	operator bool() const{
		return (fRef && fRef->fPointer != nullptr);
	}
	operator bool() {
		return (fRef && fRef->fPointer != nullptr);
	}
    
	/** @brief Returns the counter value */
	int Count() const
	{
        if(!fRef) return 0;
		return *(fRef->fCounter);
	}
    template<typename R, typename T2>
    friend TPZAutoPointer<R> TPZAutoPointerDynamicCast(TPZAutoPointer<T2> in);

private:
    /** @brief Method for deleting the reference*/
    inline void Release(){
        std::scoped_lock<std::recursive_mutex> lck(
                    pzinternal::g_ap_mut);
        //just checking if nobody else deleted fRef
        //or another fRef pointing to the same address
        if(fRef && fRef->fCounter) {
            // std::cout<<__PRETTY_FUNCTION__<<'\n';
            // std::cout<<"fref address: "<<fRef<<std::endl;
            // std::cout<<"fcount address: "<<fRef->fCounter<<std::endl;
            // std::cout<<"fcount value: "<<*(fRef->fCounter)<<std::endl;
            // std::cout<<"fpointer address: "<<fRef->fPointer<<std::endl;
            delete fRef;}
        fRef = nullptr;
    }
};

// this is the reason why the counter is dynamically allocated
template<typename R, typename T>
TPZAutoPointer<R> TPZAutoPointerDynamicCast(TPZAutoPointer<T> in) {
    TPZAutoPointer<R> rv;
        rv.fRef->fPointer = p;
		delete rv.fRef->fCounter;
		rv.fRef->fCounter = in.fRef->fCounter;
    if (R* p; (p = dynamic_cast<R*> (in.fRef->fPointer)) ) {
        rv.fRef->Increment();
    }
    return rv;
}

/** @} */
/** The following templates will be properly instantiated in
the file tpzautopointer.cpp.

These extern template instantiations are just to avoid each translation
units creating their own symbols (and afterwards throwing them away)*/
template <class T>
class TPZMatrix;


extern template class TPZAutoPointer<TPZMatrix<float>>;
extern template class TPZAutoPointer<TPZMatrix<double>>;
extern template class TPZAutoPointer<TPZMatrix<long double>>;
#include <complex>
extern template class TPZAutoPointer<TPZMatrix<std::complex<float> >>;
extern template class TPZAutoPointer<TPZMatrix<std::complex<double> >>;
extern template class TPZAutoPointer<TPZMatrix<std::complex<long double> >>;

template<class T>
class TPZFunction;
extern template class TPZAutoPointer<TPZFunction<float>>;
extern template class TPZAutoPointer<TPZFunction<double>>;
extern template class TPZAutoPointer<TPZFunction<long double>>;
extern template class TPZAutoPointer<TPZFunction<std::complex<float> >>;
extern template class TPZAutoPointer<TPZFunction<std::complex<double> >>;
extern template class TPZAutoPointer<TPZFunction<std::complex<long double> >>;

class TPZGeoMesh;
extern template class TPZAutoPointer<TPZGeoMesh>;
class TPZCompMesh;
extern template class TPZAutoPointer<TPZCompMesh>;
class TPZRefPattern;
extern template class TPZAutoPointer<TPZRefPattern>;

#endif
