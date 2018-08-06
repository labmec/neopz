#include "pzthreadpool.h"

#include <iostream>

/// if there can only be one instance of threadpool, these data objects should be members of the threadpool class
// in order to make them local scope (phil) tagged the static keyword to it

static TPZPriorityQueue<TPZTask*, std::vector<TPZTask*>, TPZTaskOrdering> globalTasksQueue;
static std::condition_variable globalTaskAvailableCond;
static int globalMinPriority = std::numeric_limits<int>::max();
static int globalMaxPriority = std::numeric_limits<int>::min();

void TPZThreadPool::updatePriorities() {
    if (globalTasksQueue.size() != 0) {
        globalMaxPriority = globalTasksQueue.top()->priority();
    } else {
        globalMaxPriority = std::numeric_limits<int>::min();
        globalMinPriority = std::numeric_limits<int>::max();
    }
}

void TPZThreadPool::threadsLoop() {
    while (true) {
        TPZAutoPointer<TPZTask> task;
        {
            std::unique_lock<std::mutex> lock(globalTasksQueue.mMutex);
            globalTaskAvailableCond.wait(lock, [this] {
                return mStop || mThreadsToDelete != 0 || globalTasksQueue.size() != 0;
            });
            if (mStop) {
                while (globalTasksQueue.size() != 0) {
                    task = globalTasksQueue.popTop(); // combined with TPZAutoPointer, deletes the tasks.
                }
                return;
            }
            if (globalTasksQueue.size() == 0 || !globalTasksQueue.top()->mSystemTask) {
                std::function<void(void) > thread_join_task;
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
                if (thread_join_task) {
                    runSystemTask(std::numeric_limits<int>::max(), thread_join_task);
                    return;
                }
            }
            task = globalTasksQueue.popTop();
            updatePriorities();
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
            globalTaskAvailableCond.notify_one();
        }
    }
    for (int i = 0; i < threads_to_create; ++i) {
        //std::cout << "Creating thread " << i + 1 << " of " << threads_to_create << " totalizing " << threadCount() + 1 << " active and " << mThreads.size() +1 << " total." << std::endl;
        mThreads.emplace_back([this] {
            threadsLoop();
        });
    }
}

TPZThreadPool::TPZThreadPool() : mThreadsToDelete(0), mZombieThreads(0), mStop(false) {
    SetNumThreads(std::thread::hardware_concurrency());
}

TPZThreadPool::~TPZThreadPool() {
    {
        std::unique_lock<std::mutex> lock(globalTasksQueue.mMutex);
        mStop = true;
    }
    globalTaskAvailableCond.notify_all();
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
    return globalMaxPriority;
}

int TPZThreadPool::minPriority() const {
    return globalMinPriority;
}

int TPZThreadPool::threadCount() const {
    return ActualThreadCount() - mThreadsToDelete - mZombieThreads;
}

int TPZThreadPool::ActualThreadCount() const {
    return mThreads.size();
}

void TPZThreadPool::appendTaskToQueue(const int priority, const TPZAutoPointer<std::packaged_task<void ()> > task, bool system_task) {
    TPZTask *newTask = new TPZTask(priority, task);
    newTask->mSystemTask = system_task;
    globalTasksQueue.addItem(newTask);
    globalTaskAvailableCond.notify_one();
}

void TPZThreadPool::checkForMaxAndMinPriority(const int priority) {
    if (priority > globalMaxPriority) {
        globalMaxPriority = priority;
    }
    if (priority < globalMinPriority) {
        globalMinPriority = priority;
    }
}

TPZThreadPool &TPZThreadPool::globalInstance() {
    static TPZThreadPool globalIntstance;
    return globalIntstance;
}
