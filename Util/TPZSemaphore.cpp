/** 
 * @file
 * @brief Contains the implementation of the methods to TPZSemaphore class.
 */
//#include "StdAfx.h"
#include "TPZSemaphore.h"
#include <fstream>

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.util.semaphore"));
#endif

TPZSemaphore::TPZSemaphore(void)
{
	fCounter = 0;
	pthread_mutex_init(&fMutex, NULL);
	pthread_cond_init(&fCond, NULL);
}

TPZSemaphore::~TPZSemaphore(void)
{
	pthread_mutex_destroy(&fMutex);
	pthread_cond_destroy(&fCond);
}

TPZSemaphore::TPZSemaphore(int initCount)
{
	fCounter = initCount;
	pthread_mutex_init(&fMutex, NULL);
	pthread_cond_init(&fCond, NULL);
}

void TPZSemaphore::Wait()
{
	pthread_mutex_lock(&fMutex);
	if (fCounter > 0)
	{
		fCounter--;
#ifdef LOG4CXX
		if(logger->isDebugEnabled())
		{
			std::stringstream sout;
			sout << "THREAD IN SEMAPHORE WAIT: " << pthread_self() <<" " <<__LINE__ << std::endl;
			sout << "FCOUNTER VALUE : " << fCounter << std::endl;
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		pthread_mutex_unlock(&fMutex);
		return;
	}
	while (fCounter == 0) {
		pthread_cond_wait(&fCond, &fMutex);
		if (fCounter > 0) {
			fCounter--;
#ifdef LOG4CXX
			if(logger->isDebugEnabled())
			{
				std::stringstream sout;
				sout << "THREAD IN SEMAPHORE AFTER WAIT: " << pthread_self() <<" " << __LINE__ << std::endl;
				sout << "FCOUNTER VALUE : " << fCounter << std::endl;
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			pthread_mutex_unlock(&fMutex);
			return;
		}
	}
	
}

void TPZSemaphore::Post()
{
	pthread_mutex_lock(&fMutex);
	fCounter++;
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "THREAD IN SEMAPHORE POST: " << pthread_self() <<" " << __LINE__ << std::endl;
		sout << "FCOUNTER VALUE : " << fCounter << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	pthread_cond_signal(&fCond);
	pthread_mutex_unlock(&fMutex);
	return;
}
