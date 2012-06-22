/** 
 * @file
 * @brief Contains the implementation of the methods to TPZSemaphore class.
 */

#include "TPZSemaphore.h"

#include <fstream>

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.util.semaphore"));
#endif

TPZSemaphore::TPZSemaphore(void)
{
	fCounter = 0;
        PZ_PTHREAD_MUTEX_INIT(&fMutex, NULL, "TPZSemaphore::TPZSemaphore()");
        PZ_PTHREAD_COND_INIT(&fCond, NULL, "TPZSemaphore::TPZSemaphore()");
}

TPZSemaphore::~TPZSemaphore(void)
{
        PZ_PTHREAD_MUTEX_DESTROY(&fMutex, "TPZSemaphore::~TPZSemaphore()");
        PZ_PTHREAD_COND_DESTROY(&fCond, "TPZSemaphore::~TPZSemaphore()");
}

TPZSemaphore::TPZSemaphore(int initCount)
{
	fCounter = initCount;
        PZ_PTHREAD_MUTEX_INIT(&fMutex, NULL, "TPZSemaphore::TPZSemaphore(int)");
        PZ_PTHREAD_COND_INIT(&fCond, NULL, "TPZSemaphore::TPZSemaphore(int)");
}

void TPZSemaphore::Wait()
{
        PZ_PTHREAD_MUTEX_LOCK(&fMutex, "TPZSemaphore::Wait()");
	if (fCounter > 0)
	{
		fCounter--;
#ifdef LOG4CXX
		if(logger->isDebugEnabled())
		{
			std::stringstream sout;
#ifdef VC
			sout << "THREAD IN SEMAPHORE WAIT: " << pthread_self().x <<" " <<__LINE__ << std::endl;
#else
			sout << "THREAD IN SEMAPHORE WAIT: " << pthread_self() <<" " <<__LINE__ << std::endl;
#endif
			sout << "FCOUNTER VALUE : " << fCounter << std::endl;
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		PZ_PTHREAD_MUTEX_UNLOCK(&fMutex, "TPZSemaphore::Wait()");
		return;
	}
	while (fCounter == 0) {
                PZ_PTHREAD_COND_WAIT(&fCond, &fMutex, "TPZSemaphore::Wait()");
		if (fCounter > 0) {
			fCounter--;
#ifdef LOG4CXX
			if(logger->isDebugEnabled())
			{
				std::stringstream sout;
#ifdef VC
				sout << "THREAD IN SEMAPHORE AFTER WAIT: " << pthread_self().x <<" " << __LINE__ << std::endl;
#else
				sout << "THREAD IN SEMAPHORE AFTER WAIT: " << pthread_self() <<" " << __LINE__ << std::endl;
#endif
				sout << "FCOUNTER VALUE : " << fCounter << std::endl;
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			PZ_PTHREAD_MUTEX_UNLOCK(&fMutex, "TPZSemaphore::Wait()");
			return;
		}
	}
	
}

void TPZSemaphore::Post()
{
	PZ_PTHREAD_MUTEX_LOCK(&fMutex, "TPZSemaphore::Post()");
	fCounter++;
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
#ifdef VC
		sout << "THREAD IN SEMAPHORE POST: " << pthread_self().x <<" " << __LINE__ << std::endl;
#else
		sout << "THREAD IN SEMAPHORE POST: " << pthread_self() <<" " << __LINE__ << std::endl;
#endif
		sout << "FCOUNTER VALUE : " << fCounter << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	PZ_PTHREAD_COND_SIGNAL(&fCond, "TPZSemaphore::Post()");
	PZ_PTHREAD_MUTEX_UNLOCK(&fMutex, "TPZSemaphore::Post()");
	return;
}
