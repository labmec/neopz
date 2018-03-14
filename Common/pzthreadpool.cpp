#include "pzthreadpool.h"

PZThreadPool::PZThreadPool()
{

}

PZThreadPool &PZThreadPool::globalInstance() {
    static PZThreadPool globalIntstance;
    return globalIntstance;
}

void PZThreadPool::setMaxNumberOfThreads(const int numberOfThreads) {
    // TODO add more thread to the threads vec if numberOfThreads greater than the actual vec.size()
    // TODO remove threads from the vec if numberOfThreads is less than the actual value
}
