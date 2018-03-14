#ifndef PZTHREADPOOL_H
#define PZTHREADPOOL_H

#include <pzvec.h>
#include <pz_pthread.h>

class PZThreadPool
{
public:
    static PZThreadPool &globalInstance() const;
    void setMaxNumberOfThreads(const int numberOfThreads);

private:
    PZThreadPool();
    TPZVec<pthread_t> threads;
};

#endif // PZTHREADPOOL_H
