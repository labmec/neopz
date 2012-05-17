/**
 * @file
 * @brief Contains declaration of the TPZSemaphore class which implements semaphore to threads.
 */

#ifndef TPZSEMAPHOREHPP
#define TPZSEMAPHOREHPP

#include "pz_pthread.h"

/**
 * @ingroup util
 * @brief Implements semaphore to threads. \ref util "Utility"
 */
class TPZSemaphore
{
private:
	/** @brief Counter of the times the semaphore is locked */
	int fCounter;
	/** @brief Mutex for the thread */
	pthread_mutex_t fMutex;
	/** @brief Condition for the thread must to be waiting */
	pthread_cond_t fCond;
	
public:
	/** @brief Default constructor */
	TPZSemaphore(void);
	/** @brief Default destructor */
	~TPZSemaphore(void);
	/** @brief Constructor with initial number for the counter */
	TPZSemaphore(int initCount);
	
	void Wait();
	void Post();
};

#endif
