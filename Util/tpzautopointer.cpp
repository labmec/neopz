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
// template < class T>
// pthread_mutex_t TPZAutoPointer<T>::gAutoCounterMutex = PTHREAD_MUTEX_INITIALIZER;

// pthread_mutex_t TPZAutoPointer<TPZSaveable>::gAutoCounterMutex = PTHREAD_MUTEX_INITIALIZER;

// template class TPZAutoPointer<TPZSaveable>;
//gAutoPointerMutex = PTHREAD_MUTEX_INITIALIZER;
