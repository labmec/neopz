/**
 * @file
 * @brief Contains declaration of the TPZSemaphore class which implements semaphore to threads.
 */
#ifndef TPZSEMAPHOREHPP
#define TPZSEMAPHOREHPP

#include <pthread.h>

/**
 * @ingroup util
 * @brief Implements semaphore to threads. \ref util "Utility"
 */
class TPZSemaphore
{
private:
	int fCounter;
	pthread_mutex_t fMutex;
	pthread_cond_t fCond;
	
public:
	TPZSemaphore(void);
	~TPZSemaphore(void);
	
	TPZSemaphore(int initCount);
	
	void Wait();
	void Post();
};

#endif