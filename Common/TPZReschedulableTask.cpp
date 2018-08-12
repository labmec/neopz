#include "pzerror.h"
#include "TPZReschedulableTask.h"
#include <utility>
#include "TPZThreadPool.h"
#include "pzlog.h"
#include <mutex>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("TPZThreadPool"));
#endif

TPZReschedulableTask::TPZReschedulableTask(const int priority, TPZAutoPointer<std::packaged_task<void(void) >> task, TPZTaskGroup *taskGroup) :
TPZTask(priority, task, taskGroup) {

}

std::shared_future<void> TPZReschedulableTask::GetFuture() {
    return mFuture;
}

void TPZReschedulableTask::startInternal() {
    TPZTask::start();
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
