#include "pzerror.h"
#include "TPZReschedulableTask.h"
#include <utility>
#include "pzthreadpool.h"
#include "pzlog.h"
#include <mutex>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("TPZThreadPool"));
#endif

TPZReschedulableTask::TPZReschedulableTask(const int priority, TPZAutoPointer<std::packaged_task<void(void) >> task) :
TPZTask(priority, task),
mState(EProcessingState::CREATED) {

}

std::shared_future<void> TPZReschedulableTask::GetFuture() {
    return mFuture;
}

void TPZReschedulableTask::startInternal() {
    mState = TPZReschedulableTask::STARTED;
    TPZTask::start();
    mState = TPZReschedulableTask::FINISHED;
}

void TPZReschedulableTask::start() {
    bool executed = false;
    {
        std::unique_lock<std::mutex> lock(mStateMutex);
        if (mState <= TPZReschedulableTask::SCHEDULED) {
            this->startInternal();
            executed = true;
        }
    }
    if (executed) {
        mCondition.notify_all();
    }
}
