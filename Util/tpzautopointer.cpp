/** 
 * @file 
 * @brief Initialize the pthread_mutex_t gAutoPointerMutex. 
 */

#include "tpzautopointer.h"

pthread_mutex_t gAutoPointerMutexArray[AP_MUTEX_ARRAY_SZ];

#ifdef PROFILE_AP_MUTEXES
  uint64_t ap_mutex_accesses[AP_MUTEX_ARRAY_SZ];
#endif

/* Class to initialize the array. */
class AutoPointerMutexArrayInit
{
public:
  AutoPointerMutexArrayInit() 
  {
    for (int i=0; i< AP_MUTEX_ARRAY_SZ; i++) {
      pthread_mutex_init(&(gAutoPointerMutexArray[i]),0);
#ifdef PROFILE_AP_MUTEXES
      std::cout << "Inicializando autopointe mutex " << i << std:endl;
      ap_mutex_accesses[i] = 0;
#endif
    }
  }
  ~AutoPointerMutexArrayInit() 
  {
    for (int i=0; i< AP_MUTEX_ARRAY_SZ; i++) {
      pthread_mutex_destroy(&(gAutoPointerMutexArray[i]));
#ifdef PROFILE_AP_MUTEXES
    std:cout << "Destruindo autopointer mutex " << i << std:endl;
    printf("AutoPointer Mutex 0x%p accessed %lld times\n", 
	   &(gAutoPointerMutexArray[i]), ap_mutex_accesses[i]);
#endif
    }
  }

};

AutoPointerMutexArrayInit tmp;

// template < class T>
// pthread_mutex_t TPZAutoPointer<T>::gAutoCounterMutex = PTHREAD_MUTEX_INITIALIZER;

// pthread_mutex_t TPZAutoPointer<TPZSavable>::gAutoCounterMutex = PTHREAD_MUTEX_INITIALIZER;

// template class TPZAutoPointer<TPZSavable>;
//gAutoPointerMutex = PTHREAD_MUTEX_INITIALIZER;
