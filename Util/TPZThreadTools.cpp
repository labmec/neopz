#include "TPZThreadTools.h"
#include "pzerror.h"
#include "pz_pthread.h"
/// namespace tht means "ThreadTools"
namespace tht
{

	//// critical section ///

	void InitializeCriticalSection(pz_critical_section_t & cs)
	{
#ifdef using_pthread
		pthread_mutexattr_t * nullPtr = NULL;
		pthread_mutex_init(&cs, nullPtr);
#endif
#ifdef WINDOWS_THREADS
		::InitializeCriticalSection(&cs);
#endif
#ifdef EMBARCADERO_THREADS
		cs = new TCriticalSection();
#endif
	}

	void DeleteCriticalSection(pz_critical_section_t & cs)
	{
#ifdef using_pthread
		pthread_mutex_destroy(&cs);
#endif
#ifdef WINDOWS_THREADS
		::DeleteCriticalSection(&cs);
#endif
#ifdef EMBARCADERO_THREADS
		// ::DeleteCriticalSection( &cs );
		///nothing here: autopointer will destroy it
#endif
	}

	void EnterCriticalSection(pz_critical_section_t & cs)
	{
#ifdef using_pthread
		pthread_mutex_lock(&cs);
#endif
#ifdef WINDOWS_THREADS
		::EnterCriticalSection(&cs);
#endif
#ifdef EMBARCADERO_THREADS
		// ::EnterCriticalSection(&cs);
		cs->Enter();
#endif
	}

	void LeaveCriticalSection(pz_critical_section_t & cs)
	{
#ifdef using_pthread
		pthread_mutex_unlock(&cs);
#endif
#ifdef WINDOWS_THREADS
		::LeaveCriticalSection(&cs);
#endif
#ifdef EMBARCADERO_THREADS
		// ::LeaveCriticalSection(&cs);
		cs->Leave();
#endif
	}

	//// mutex ///

	void InitializeMutex(pz_mutex_t & m)
	{
#ifdef using_pthread
		pthread_mutexattr_t * nullPtr = NULL;
		pthread_mutex_init(&m, nullPtr);
#endif
#ifdef WINDOWS_THREADS
		m = CreateMutex(NULL, false, NULL);
		if(!m)
		{
			MessageBox(NULL, "Error", "Fail to create Windows Mutex",
				MB_OK | MB_ICONERROR);
			DebugStop();
		}
#endif
#ifdef EMBARCADERO_THREADS
		m = new TMutex(false);
#endif
	}

	void DeleteMutex(pz_mutex_t & m)
	{
#ifdef using_pthread
		pthread_mutex_destroy(&m);
#endif
#ifdef WINDOWS_THREADS
		CloseHandle(m);
#endif
#ifdef EMBARCADERO_THREADS
		///nothing here: autopointer will destroy it
#endif
	}

	void MutexLock(pz_mutex_t & m)
	{
#ifdef using_pthread
		pthread_mutex_lock(&m);
#endif
#ifdef WINDOWS_THREADS
		WaitForSingleObject(m, INFINITE);
#endif
#ifdef EMBARCADERO_THREADS
		m->Acquire();
#endif
	}

	void MutexUnlock(pz_mutex_t & m)
	{
#ifdef using_pthread
		pthread_mutex_unlock(&m);
#endif
#ifdef WINDOWS_THREADS
		::ReleaseMutex(m);
#endif
#ifdef EMBARCADERO_THREADS
		m->Release();
#endif
	}

	//// semaphore ///

	int InitializeSemaphore(pz_semaphore_t & s)
	{
#ifdef using_pthread
//		int sem_result = sem_init(&s, 0, 0);
//		return sem_result;
		// In TPZSemaphore nothing needs to be done (remember to create an object TPZSemaphore and not a pointer)
#endif
#ifdef WINDOWS_THREADS
		s = ::CreateSemaphore(NULL, 0, 1, NULL);
		if(!s)
		{
			MessageBox(NULL, "Error", "Fail to create Windows Mutex",
				MB_OK | MB_ICONERROR);
			DebugStop();
			return 1;
		}
		return 0;
#endif
#ifdef EMBARCADERO_THREADS
		int max = INT_MAX;
		System::UnicodeString name = "";
		s = new TSemaphore(NULL, 0, max, name, false);
		return 0;
#endif
		return 0;
	}

	void DeleteSemaphore(pz_semaphore_t & s)
	{
#ifdef using_pthread
//		sem_destroy(&s);
		// as it is and object (and not a pointer) it will be destroied at the end of the scope
		
#endif
#ifdef WINDOWS_THREADS
		CloseHandle(s);
#endif
#ifdef EMBARCADERO_THREADS
		///nothing to do here. Autopointer will handle this.
#endif
	}

	void SemaphorePost(pz_semaphore_t & s)
	{
#ifdef using_pthread
//		sem_post(&s);
		s.Post();
		//(remember to create an object TPZSemaphore and not a pointer)
#endif
#ifdef WINDOWS_THREADS
		int64_t pCount = 0;
		::ReleaseSemaphore(s, 1, &pCount);
#endif
#ifdef EMBARCADERO_THREADS
		s->Release();
#endif
	}

	void SemaphoreWait(pz_semaphore_t & s)
	{
#ifdef using_pthread
		//sem_wait(&s);
		s.Wait();
		//(remember to create an object TPZSemaphore and not a pointer)
#endif
#ifdef WINDOWS_THREADS
		WaitForSingleObject(s, INFINITE);
#endif
#ifdef EMBARCADERO_THREADS
		s->Acquire();
#endif
	}

	//// thread ///

	void CreateThread(pz_thread_t & thread, void* (*function)(void* param),
		void* data)
	{
#ifdef using_pthread
		pthread_attr_t * attr = NULL;
		int res = pthread_create(&thread, attr, function, data);
		if (res)
			DebugStop();
#endif
#ifdef WINDOWS_THREADS
		thread = ::CreateThread(NULL, 0, (uint64_t (__stdcall* )(void *))function, data, 0, NULL);
#endif
#ifdef EMBARCADERO_THREADS
		thread = new TSWXEmbarcaderoThread(function, data);
#endif
	}

	void ThreadWaitFor(pz_thread_t & thread)
	{
#ifdef using_pthread
		pthread_join(thread, NULL);
#endif
#ifdef WINDOWS_THREADS
		WaitForSingleObject(thread, INFINITE);
#endif
#ifdef EMBARCADERO_THREADS
		thread->WaitFor();
#endif
	}

}
///namespace
