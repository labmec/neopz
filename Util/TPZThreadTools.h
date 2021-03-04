//---------------------------------------------------------------------------
#ifndef TPZThreadToolsH
#define TPZThreadToolsH

#define using_pthread 
#ifdef using_pthread
#include "pthread.h"
#include "TPZSemaphore.h"
#endif
#ifdef WINDOWS_THREADS
#ifndef NOMINMAX
#define NOMINMAX // Preventing the redefinition of min and max as macros
#endif
	#include <windows.h>
#endif
#ifdef EMBARCADERO_THREADS
	#include "System.hpp"
	#include "SyncObjs.hpp"
	#include "tpzautopointer.h"
	#include "TSWXEmbarcaderoThread.h"
#endif


#ifdef using_pthread
typedef  pthread_mutex_t pz_critical_section_t;
typedef  pthread_mutex_t pz_mutex_t;
typedef TPZSemaphore pz_semaphore_t;
typedef pthread_t pz_thread_t;
#endif
#ifdef WINDOWS_THREADS
typedef CRITICAL_SECTION pz_critical_section_t;
typedef HANDLE pz_mutex_t;
typedef HANDLE pz_semaphore_t;
typedef HANDLE pz_thread_t;
#endif
#ifdef EMBARCADERO_THREADS
typedef TPZAutoPointer<TCriticalSection> pz_critical_section_t;
typedef TPZAutoPointer< TMutex > pz_mutex_t;
typedef TPZAutoPointer<TSemaphore> pz_semaphore_t;
typedef TPZAutoPointer<TThread> pz_thread_t;
#endif

// namespace tht means "ThreadTools"
namespace tht{

  ///Critical sections
  void InitializeCriticalSection(pz_critical_section_t & cs);

  void DeleteCriticalSection(pz_critical_section_t & cs);

  void EnterCriticalSection(pz_critical_section_t & cs);

  void LeaveCriticalSection(pz_critical_section_t & cs);

  ///Mutexes
  void InitializeMutex(pz_mutex_t & m);

  void DeleteMutex(pz_mutex_t & m);

  void MutexLock(pz_mutex_t & m);

  void MutexUnlock(pz_mutex_t & m);

  ///Semaforos
  int InitializeSemaphore(pz_semaphore_t & s  );

  void DeleteSemaphore(pz_semaphore_t & s  );

	void SemaphorePost( pz_semaphore_t & s  );

  void SemaphoreWait( pz_semaphore_t & s );

  ///Threads
	void CreateThread( pz_thread_t & thread, void* (*function)(void* param), void* data);

  void ThreadWaitFor(pz_thread_t & thread);

}///namespace
#endif
