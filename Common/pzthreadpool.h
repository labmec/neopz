#ifndef PZTHREADPOOL_H
#define PZTHREADPOOL_H

#include <vector>
#include "pzpriorityqueue.h"
#include <limits>

#include <cstdarg>
#include <functional>

#include "tpzautopointer.h"

#include <thread>
#include <future>
#include <condition_variable>
#include <mutex>

// Helper class for ordering the tasks that the user have requested

class TPZTask {
public:

    TPZTask(const int priority, TPZAutoPointer<std::packaged_task<void(void) >> task) :
    mSystemTask(false),
    mPriority(priority),
    mTask(task) {

    }

    int priority() const {
        return mPriority;
    }

    void start() const {
        (*mTask)();
    }

    friend class TPZTaskOrdering;
    friend class TPZThreadPool;
    
private:
    bool mSystemTask;
    const int mPriority;
    TPZAutoPointer<std::packaged_task<void(void) >> mTask;
};

// Simple struct needed by std::priority_queue for ordering the items

struct TPZTaskOrdering {

    bool operator()(const TPZTask *lhs, const TPZTask *rhs) {
        if (lhs->mSystemTask) {
            if (rhs->mSystemTask) {
                return lhs->priority() < rhs->priority();
            } else {
                return false;
            }
        } else if (rhs->mSystemTask){
                return true;
        } else {
                return lhs->priority() < rhs->priority();
        }
    }
};

class TPZThreadPool {
public:
    static TPZThreadPool &globalInstance();

    template<typename... Args>
    std::future<void> run(const int priority, std::function<void(Args...) > func, Args... args) {
        checkForMaxAndMinPriority(priority);
        TPZAutoPointer < std::packaged_task<void(void) >> task = new std::packaged_task<void(void) >(std::bind(func, args...));
        if (threadCount() != 0){
            appendTaskToQueue(priority, task, false);
        } else {
            (*task)();
        }
        return task->get_future();
    }
    
    void SetNumThreads(const unsigned numThreads);
    int maxPriority() const;
    int minPriority() const;
    int threadCount() const;

private:
    TPZThreadPool();
    int ActualThreadCount() const;
    void threadsLoop();
    void updatePriorities();
    void appendTaskToQueue(const int priority, const TPZAutoPointer<std::packaged_task<void(void) >> task, const bool system_task);
    void checkForMaxAndMinPriority(const int priority);
    ~TPZThreadPool();
    
    template<typename... Args>
    std::future<void> runSystemTask(const int priority, std::function<void(Args...) > func, Args... args) {
        checkForMaxAndMinPriority(priority);
        TPZAutoPointer < std::packaged_task<void(void) >> task = new std::packaged_task<void(void) >(std::bind(func, args...));
        appendTaskToQueue(priority, task, true);
        return task->get_future();
    }
    
    
    std::vector<std::thread> mThreads;
    std::mutex mThreadsMutex;
    unsigned int mThreadsToDelete;
    unsigned int mZombieThreads;
    bool mStop;
    

};

#endif // PZTHREADPOOL_H
