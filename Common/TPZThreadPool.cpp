#include "pzerror.h"

#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("TPZThreadPool");
#endif

#include "TPZThreadPool.h"

#include <exception>
#include <iostream>
#include "TPZTask.h"
#include <mutex>


#include "TPZReschedulableTask.h"

void TPZThreadPool::updatePriorities() {
    if (mTasksQueue.size() != 0) {
        mMaxPriority = mTasksQueue.top()->priority();
    } else {
        mMaxPriority = std::numeric_limits<int>::min();
        mMinPriority = std::numeric_limits<int>::max();
    }
}

void TPZThreadPool::threadsLoop() {
    while (true) {
        TPZAutoPointer<TPZTask> task;
        std::function<void(void) > thread_join_task;
        {
            std::unique_lock<std::mutex> lock(mTasksQueue.mMutex);
            mTaskAvailableCond.wait(lock, [this] {
                return mStop || mThreadsToDelete != 0 || mTasksQueue.size() != 0;
            });
            if (mStop) {
                while (mTasksQueue.size() != 0) {
                    task = mTasksQueue.popTop();
                    task->Cancel();
                }
                return;
            }
            if (mTasksQueue.size() == 0 || !mTasksQueue.top()->mSystemTask) {
                {
                    if ((mThreadsToDelete != 0)) { // It may seem odd to check this so soon, but mind the cost of locking a mutex
                        std::unique_lock<std::mutex> lock(mThreadsMutex);
                        if ((mThreadsToDelete != 0)) {
                            typename std::thread::id thread_id = std::this_thread::get_id();
                            thread_join_task = [this, thread_id] {
                                std::unique_lock<std::mutex> lock(mThreadsMutex);
                                unsigned int num_threads = ActualThreadCount();
                                for (unsigned int i = 0; i < num_threads; ++i) {
                                    if (mThreads[i].get_id() == thread_id) {
                                        mThreads[i].join();
                                        if (i != num_threads - 1) {
                                            mThreads[i].swap(mThreads[num_threads - 1]);
                                        }
                                        mThreads.pop_back();
                                        --mZombieThreads;
                                        //std::cout << "Deleted thread " << thread_id << " totalizing " << mThreads.size()-mZombieThreads << " active and " << mThreads.size() << " total." << std::endl;
                                        return;
                                    }
                                }
                            };
                            --mThreadsToDelete;
                            ++mZombieThreads;
                            //std::cout << "Scheduling deletion of thread " << thread_id << " totalizing " << mThreads.size()-mZombieThreads << " active and " << mThreads.size() << " total." << std::endl;
                        }
                    }
                }
            }
            if (!thread_join_task) {
                task = mTasksQueue.popTop();
                updatePriorities();
            }
        }
        if (thread_join_task) {
            runSystemTask(std::numeric_limits<int>::max(), thread_join_task);
            return ;
        }
        if (task) {
            task->start();
        }
    }
}

void TPZThreadPool::SetNumThreads(const unsigned numThreads) {
    std::unique_lock<std::mutex> lock(mThreadsMutex);
    unsigned int num_threads_before = threadCount();
    int threads_to_create = numThreads - num_threads_before;
    if (threads_to_create < 0) {
        mThreadsToDelete -= threads_to_create;
        for (unsigned int i = 0; i < -threads_to_create; ++i) {
            mTaskAvailableCond.notify_one();
        }
    }
    for (int i = 0; i < threads_to_create; ++i) {
        //std::cout << "Creating thread " << i + 1 << " of " << threads_to_create << " totalizing " << threadCount() + 1 << " active and " << mThreads.size() +1 << " total." << std::endl;
        mThreads.emplace_back([this] {
            threadsLoop();
        });
    }
}

TPZThreadPool::TPZThreadPool() : mThreadsToDelete(0), mZombieThreads(0), mStop(false), mTasksQueue(), mTaskAvailableCond(),
mMinPriority(std::numeric_limits<int>::max()),
mMaxPriority(std::numeric_limits<int>::min()){
    //SetNumThreads(std::thread::hardware_concurrency());
}

std::shared_future<void> TPZThreadPool::run(const int priority, TPZAutoPointer<std::packaged_task<void(void)> >& task, TPZTaskGroup *taskGroup) {
    std::shared_future<void> fut = task->get_future().share();
    if (threadCount() != 0) {
        appendTaskToQueue(priority, task, false, taskGroup);
    } else {
        (*task)();
    }
    return fut;
}

void TPZThreadPool::reschedule(const int priority, TPZAutoPointer<TPZReschedulableTask> &task) {
    TPZAutoPointer<TPZTask> autoPointerTask(TPZAutoPointerDynamicCast<TPZTask>(task));
    {
        std::unique_lock<std::mutex> lock(mTasksQueue.mMutex);
        mTasksQueue.remove(autoPointerTask);
    }
    task->mPriority = priority;
    appendTaskToQueue(autoPointerTask);
}

TPZThreadPool::~TPZThreadPool() {
    {
        std::unique_lock<std::mutex> lock(mTasksQueue.mMutex);
        mStop = true;
    }
    mTaskAvailableCond.notify_all();
    {
        std::unique_lock<std::mutex> lock(mThreadsMutex);
        for (auto &thread : mThreads) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    }
}

int TPZThreadPool::maxPriority() const {
    return mMaxPriority;
}

int TPZThreadPool::minPriority() const {
    return mMinPriority;
}

int TPZThreadPool::threadCount() const {
    return ActualThreadCount() - mThreadsToDelete - mZombieThreads;
}

int TPZThreadPool::ActualThreadCount() const {
    return mThreads.size();
}

void TPZThreadPool::appendTaskToQueue(TPZAutoPointer<TPZTask> &task) {
    std::unique_lock<std::mutex> lock(mTasksQueue.mMutex);
    mTasksQueue.addItem(task);
    mTaskAvailableCond.notify_one();
}

TPZAutoPointer<TPZTask> TPZThreadPool::appendTaskToQueue(const int priority, TPZAutoPointer<std::packaged_task<void (void)>> &task, bool system_task, TPZTaskGroup *taskGroup) {
    TPZAutoPointer<TPZTask> newTask(new TPZTask(priority, task, taskGroup));
    newTask->mSystemTask = system_task;
    appendTaskToQueue(newTask);
    return newTask;
}

void TPZThreadPool::checkForMaxAndMinPriority(const int priority) {
    if (priority > mMaxPriority) {
        mMaxPriority = priority;
    }
    if (priority < mMinPriority) {
        mMinPriority = priority;
    }
}

std::shared_future<void> TPZThreadPool::runNow(TPZAutoPointer<TPZReschedulableTask> &task) {
    std::unique_lock<std::mutex> lock(task->mStateMutex);
    switch (task->mState) {
        case TPZReschedulableTask::EProcessingState::CREATED:
            task->mFuture = task->mTask->get_future().share();
            if (threadCount() != 0) {
                TPZAutoPointer<TPZTask> tpztask = TPZAutoPointerDynamicCast<TPZTask>(task);
                tpztask->mPriority = std::numeric_limits<int>::max();
                appendTaskToQueue(tpztask);
                task->mState = TPZReschedulableTask::EProcessingState::SCHEDULED;
                task->mCondition.wait(lock);
            } else {
                task->startInternal();
                task->mCondition.notify_all();
            }
            break;
        case TPZReschedulableTask::SCHEDULED:
            TPZThreadPool::globalInstance().reschedule(std::numeric_limits<int>::max(), task);
            task->mCondition.wait(lock);
            break;
        case TPZReschedulableTask::STARTED:
            task->mCondition.wait(lock);
            break;
        case TPZReschedulableTask::FINISHED:
            break;
    }
    return task->mFuture;
}

void TPZThreadPool::run(TPZAutoPointer<TPZReschedulableTask> &task) {
    std::unique_lock<std::mutex> lock(task->mStateMutex);
    switch (task->mState) {
        case TPZReschedulableTask::EProcessingState::CREATED:
        {
            task->mFuture = task->mTask->get_future().share();
            if (threadCount() != 0) {
                TPZAutoPointer<TPZTask> tpztask = TPZAutoPointerDynamicCast<TPZTask>(task);
                appendTaskToQueue(tpztask);
                task->mState = TPZReschedulableTask::EProcessingState::SCHEDULED;
            } else {
                task->startInternal();
                task->mCondition.notify_all();
            }
            break;
        }
        case TPZReschedulableTask::EProcessingState::SCHEDULED:
        case TPZReschedulableTask::EProcessingState::STARTED:
        case TPZReschedulableTask::EProcessingState::FINISHED:
            break;
    }
}


TPZThreadPool &TPZThreadPool::globalInstance() {
    static TPZThreadPool globalIntstance;
    return globalIntstance;
}
