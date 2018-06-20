/**
 * @file
 * @brief Contains declaration of the TPZAutoPointer class which has Increment and Decrement actions are mutexed by this mutex.
 */

#ifndef TPZAUTOPOINTER_H
#define TPZAUTOPOINTER_H

#include "pz_pthread.h"

/**
 * \addtogroup util
 * @{
 */

/**
 * @brief Increment and Decrement actions are mutexed by this mutex
 */
#define AP_MUTEX_ARRAY_SZ 512

/* Define PROFILE_AP_MUTEXES if you want to profile the autopointer mutexes use. */
//#define PROFILE_AP_MUTEXES

#ifdef PROFILE_AP_MUTEXES
extern uint64_t ap_mutex_accesses[];
#endif

#define AP_MUTEX_HASH_1         \
addr = (addr >> 32) ^ addr;   \
addr = (addr >> 16) ^ addr;   \
addr = (addr >> 8)  ^ addr;   \
addr = (addr >> 4)  ^ addr;   \
i = (unsigned) (addr % AP_MUTEX_ARRAY_SZ)

#define AP_MUTEX_HASH_2         \
addr = (addr >> 8)  ^ addr;   \
addr = (addr >> 4)  ^ addr;   \
i = (unsigned) (addr % AP_MUTEX_ARRAY_SZ)

extern pthread_mutex_t gAutoPointerMutexArray[];
inline pthread_mutex_t* get_ap_mutex(void* obj)
{
	unsigned i;
	uint64_t addr = (uint64_t) obj;
	//  AP_MUTEX_HASH_1;
	AP_MUTEX_HASH_2;
#ifdef PROFILE_AP_MUTEXES
	ap_mutex_accesses[i]++;
#endif
	return &(gAutoPointerMutexArray[i]);
}

#include <stdlib.h>



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
        int *fCounter;
        
        TPZReference()
        {
            fPointer = 0;
            fCounter = new int;
            (*fCounter) = 1;
        }
        
        TPZReference(T2 *pointer)
        {
            fPointer = pointer;
            fCounter = new int;
            (*fCounter) = 1;
        }
        
        ~TPZReference()
        {
            if(fPointer) {
                delete fPointer;
            }
            fPointer = 0;
            if(fCounter) {
                delete fCounter;
            }
            fCounter = 0;
        }
        
        void ReallocForNuma(int node_id)
        {
            // -2: Do not realloc
            // -1: Realloc to the node that is currently running this thread.
            // node_id>=0 : Realloc to node_id
            
            if (node_id == -2) return;
            
            //TODO: Currently relying on first-touch policy to implement case -1
            
            T2* old = fPointer;
            fPointer = (T2*) fPointer->Clone();
            delete old;
        }
        

        /** @brief Increment the counter */
        bool Increment()
        {
            if(PZ_PTHREAD_MUTEX_LOCK(get_ap_mutex((void*) this), __PRETTY_FUNCTION__))
                return false;
            (*fCounter)++;
            PZ_PTHREAD_MUTEX_UNLOCK(get_ap_mutex((void*) this), __PRETTY_FUNCTION__);
            return true;
        }
        /** @brief Decrease the counter. If the counter is zero, delete myself */
        bool Decrease()
        {
            bool should_delete = false;
            if(PZ_PTHREAD_MUTEX_LOCK(get_ap_mutex((void*) this), __PRETTY_FUNCTION__))
                return false;
            (*fCounter)--;
            
            if((*fCounter) <= 0) should_delete = true;
            
            PZ_PTHREAD_MUTEX_UNLOCK(get_ap_mutex((void*) this), __PRETTY_FUNCTION__);
            if(should_delete)
            {
                delete this;
            }
            return true;
        }
        
    };
    
	/** @brief The object which contains the pointer and the reference count */
	TPZReference<T> *fRef;
    
public:
	/** @brief Creates an reference counted null pointer */
	TPZAutoPointer()
	{
		fRef = new TPZReference<T>();
	}
    
	/** @brief The destructor will delete the administered pointer if its reference count is zero */
	~TPZAutoPointer()
	{
            if (fRef){
		fRef->Decrease();
            }
	}
    
	/** @brief This method will create an object which will administer the area pointed to by obj */
	TPZAutoPointer(T *obj)
	{
		fRef = new TPZReference<T>(obj);
	}
    
	/** @brief Share the pointer of the copy */
	TPZAutoPointer(const TPZAutoPointer<T> &copy)
	{
		fRef = copy.fRef;
		fRef->Increment();
	}
        
	/** @brief Move assignment operator */
	TPZAutoPointer &operator=(TPZAutoPointer<T> &&copy){
            if (fRef) {
                fRef->Decrease();
            }
            fRef = copy.fRef;
            copy.fRef = NULL;
            return *this;
        }
        
	/** @brief Assignment operator */
	TPZAutoPointer &operator=(const TPZAutoPointer<T> &copy)
	{
		if(copy.fRef == fRef) return *this;
		copy.fRef->Increment();
		fRef->Decrease();
		fRef = copy.fRef;
		return *this;
	}
    
	/** @brief Returns the referenced object */
	operator T&()
	{
		return *(fRef->fPointer);
	}
        
	/** @brief Returns the referenced object */
	T& operator *() const
	{
		return *(fRef->fPointer);
	}
        
	/** @brief Returns the referenced object */
	T& operator *()
	{
		return *(fRef->fPointer);
	}
    
	/** @brief Returns the pointer for referenced object */
	T *operator->() const
	{
		return fRef->fPointer;
	}
	T *operator->()
	{
		return fRef->fPointer;
	}
    
	void ReallocForNuma(int node)
	{
		fRef->ReallocForNuma(node);
	}
    
	/** @brief Returns if pointer was attributed */
	operator bool() const{
		return (this->fRef->fPointer != 0);
	}
	operator bool() {
		return fRef->fPointer != 0;
	}
    
	/** @brief Returns the counter value */
	int Count() const
	{
		return *(fRef->fCounter);
	}
	int Count()
	{
		return *(fRef->fCounter);
	}
    template<typename R, typename T2>
    friend TPZAutoPointer<R> TPZAutoPointerDynamicCast(TPZAutoPointer<T2> in);
};
        
template<typename R, typename T>
TPZAutoPointer<R> TPZAutoPointerDynamicCast(TPZAutoPointer<T> in) {
    TPZAutoPointer<R> rv;
    R* p;
    if ( (p = dynamic_cast<R*> (in.operator->())) ) {
        rv.fRef->fPointer = dynamic_cast<R*> (in.fRef->fPointer);
		delete rv.fRef->fCounter;
		rv.fRef->fCounter = in.fRef->fCounter;
        rv.fRef->Increment();
    }
    return rv;
}

/** @} */

#endif
