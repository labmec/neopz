/**
 * @file
 * @brief Contains declaration of the TPZAutoPointer class.
 */

#ifndef TPZAUTOPOINTER_H
#define TPZAUTOPOINTER_H

#include <stdlib.h>
#include <atomic>
#include <thread>
#include <chrono>
#include <iostream>
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

//namespace pzinternal{
//    extern std::recursive_mutex g_ap_mut;
//extern std::mutex g_diag_mut;
//}

/**
 * @brief This class implements a reference counter mechanism to administer a dynamically allocated object. \ref util "Utility"
 * @author Philippe R. B. Devloo
 */
template<class T>
class TPZAutoPointer {
    
	/** @brief Counter struct */
    struct TPZReference
    {
        /** @brief Pointer to T object */
        T *fPointer{nullptr};
        /** @brief Reference counter*/
        std::atomic_int fCounter{0};
        
        TPZReference()
        {
            fPointer = nullptr;
            fCounter.store(1);
        }
        
        TPZReference(T *pointer)
        {
            fPointer = pointer;
            fCounter.store(1);
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
        }
        

        /** @brief Increment the counter */
        bool Increment()
        {
//            std::mutex *mut = get_ap_mutex((void*) this);
//            std::unique_lock<std::mutex> lck(*mut);//already locks
            fCounter.fetch_add(1);
//            lck.unlock();
            return true;
        }
        /** @brief Decrease the counter. The return value will specify if
            it is needed to delete the reference.*/
        bool Decrease()
        {
            bool should_delete = false;
            int result = fCounter.fetch_sub(1);
            result--;
//            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            // at this point the object may already have been deleted by another thread
//            {
//                const std::lock_guard<std::mutex> lock(pzinternal::g_diag_mut);
//                std::cout <<  "thread id " <<  std::this_thread::get_id() <<
//                " result " << result << std::endl;
//            }
            if((result) == 0) should_delete = true;
            else if((result) < 0){
                throw std::logic_error(
                    "Invalid value for ref counter of TPZAutoPointer");
            }
            return should_delete;
        }
        
    };
    
	/** @brief The object which contains the pointer and the reference count */
	TPZReference *fRef{nullptr};
    
public:
	/** @brief Creates an reference counted null pointer */
	TPZAutoPointer() : fRef(new TPZReference())
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
	TPZAutoPointer(T *obj) : fRef(new TPZReference(obj))
	{

	}

    /** @brief This method will create an object which will administer the area pointed to by obj */
    template<class T2>
	TPZAutoPointer(T2 *obj)
	{
        if constexpr (!std::is_base_of_v<T, T2>){
            static_assert(!sizeof(T2*),"Incompatible types when creating TPZAutoPointer");
        }
        fRef = new TPZReference((T*)obj);
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
        fRef = new TPZReference(rval);
        return *this;
    }

    /** @brief Share the pointer of the rval of derived type T2 */
    template<class T2>
	TPZAutoPointer(const TPZAutoPointer<T2> &other)
	{
        if constexpr (!std::is_base_of_v<T, T2>){
            throw std::logic_error(
                    "Incompatible types when creating TPZAutoPointer");
        }
        fRef = (TPZReference*) other.fRef;
        fRef->Increment();
	}
    //! Assignment operator from TPZAutoPointer of derived type T2
    template<class T2>
    TPZAutoPointer &operator=(const TPZAutoPointer<T2> &other)
    {
        if constexpr (!std::is_base_of_v<T, T2>){
            static_assert(!sizeof(T2*),"Incompatible types when creating TPZAutoPointer");
        }
        if (fRef && fRef->Decrease()){
            Release();
        }
        fRef = (TPZReference *) other.fRef;
        fRef->Increment();
        return *this;
    }
    //! Assignment operator from  derived type T2
    template<class T2>
    TPZAutoPointer &operator=(T2 *other)
    {
        if constexpr (!std::is_base_of_v<T, T2>){
            static_assert(!sizeof(T2*),"Incompatible types when creating TPZAutoPointer");
        }
        if (fRef && fRef->Decrease()){
            Release();
        }
        fRef = new TPZReference(other);
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
		return fRef->fCounter;
	}
    template<typename R, typename T2>
    friend TPZAutoPointer<R> TPZAutoPointerDynamicCast(TPZAutoPointer<T2> in);
    template<class T2>
    friend class TPZAutoPointer;
private:
    /** @brief Method for deleting the reference*/
    inline void Release(){
        if(fRef && ! fRef->fCounter) {
            delete fRef;
        }
        fRef = nullptr;
    }
};

// this is the reason why the counter is dynamically allocated
template<typename R, typename T>
TPZAutoPointer<R> TPZAutoPointerDynamicCast(TPZAutoPointer<T> in) {
    TPZAutoPointer<R> rv;
    if (R* p; (p = dynamic_cast<R*> (in.fRef->fPointer)) ) {
        delete rv.fRef;
        rv.fRef = (typename TPZAutoPointer<R>::TPZReference *) in.fRef;
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
