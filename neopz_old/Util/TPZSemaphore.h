#ifndef TPZSEMAPHOREHPP
#define TPZSEMAPHOREHPP

#include <pthread.h>

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