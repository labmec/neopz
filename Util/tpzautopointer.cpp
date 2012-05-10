/** 
 * @file 
 * @brief Initialize the pthread_mutex_t gAutoPointerMutex. 
 */
//
// C++ Implementation: tpzautopointer
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpzautopointer.h"

pthread_mutex_t gAutoPointerMutex = PTHREAD_MUTEX_INITIALIZER;

pthread_mutex_t gAutoPointerMutexArray[AP_MUTEX_ARRAY_SZ];

#ifdef PROFILE_AP_MUTEXES
  unsigned long long ap_mutex_accesses[AP_MUTEX_ARRAY_SZ];
#endif

/* Class to initialize the array. */
class AutoPointerMutexArrayInit
{
public:
  AutoPointerMutexArrayInit() 
  {
    for (int i=0; i< AP_MUTEX_ARRAY_SZ; i++) {
      pthread_mutex_init(&(gAutoPointerMutexArray[i]),0);
      printf("Inicializando autopointe mutex %d\n", i);
#ifdef PROFILE_AP_MUTEXES
      ap_mutex_accesses[i] = 0;
#endif
    }
  }
  ~AutoPointerMutexArrayInit() 
  {
    for (int i=0; i< AP_MUTEX_ARRAY_SZ; i++) {
      pthread_mutex_destroy(&(gAutoPointerMutexArray[i]));
      printf("Destruindo autopointer mutex %d\n", i);
#ifdef PROFILE_AP_MUTEXES
      printf("AutoPointer Mutex 0x%p accessed %lld times\n", 
	     &(gAutoPointerMutexArray[i]), ap_mutex_accesses[i]);
#endif
    }
  }

};

AutoPointerMutexArrayInit tmp;

// template < class T>
// pthread_mutex_t TPZAutoPointer<T>::gAutoCounterMutex = PTHREAD_MUTEX_INITIALIZER;

// pthread_mutex_t TPZAutoPointer<TPZSaveable>::gAutoCounterMutex = PTHREAD_MUTEX_INITIALIZER;

// template class TPZAutoPointer<TPZSaveable>;
//gAutoPointerMutex = PTHREAD_MUTEX_INITIALIZER;
