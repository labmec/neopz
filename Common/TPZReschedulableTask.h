#ifndef TPZRESCHEDULABLETASK_H
#define TPZRESCHEDULABLETASK_H

#include <future>
#include <mutex>
#include "TPZTask.h"

class TPZReschedulableTask : public TPZTask {
public:
    friend class TPZThreadPool;

    std::shared_future<void> GetFuture();
    

protected:
    TPZReschedulableTask(const int priority, TPZAutoPointer<std::packaged_task<void(void)>> task);
    
    // Should be called in a thread-safe context
    virtual void startInternal();
    
    virtual void start();

    enum EProcessingState {
        CREATED,
        SCHEDULED,
        STARTED,
        FINISHED
    };
    
    EProcessingState mState;
    std::mutex mStateMutex;
private:


    std::shared_future<void> mFuture;
    std::condition_variable mCondition;
};

#endif // TPZRESCHEDULABLETASK_H
