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

/// Helper class for ordering the tasks that the user have requested
// Object of this class will be created by the TPZThreadPool class "automatically"
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

    friend struct TPZTaskOrdering;
    friend class TPZThreadPool;
    
private:
    bool mSystemTask;
    const int mPriority;
    TPZAutoPointer<std::packaged_task<void(void) >> mTask;
};

/// Simple struct needed by std::priority_queue for ordering the items
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

/// class created to administer tasks the will be executed asynchronously
class TPZThreadPool {
public:
    /// return a reference to the one and only instance of TPZThreadPool
    static TPZThreadPool &globalInstance();

    /// submit a task to be executed by threadpool
    // the return value of the std::future<void> object allows to block the calling thread if necessary
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
    
    /// sets the number of threads to be executed simultaneously
    void SetNumThreads(const unsigned numThreads);
    int maxPriority() const;
    int minPriority() const;
    /// return the number of threads currently executing
    int threadCount() const;

private:
    TPZThreadPool();
    int ActualThreadCount() const;
    void threadsLoop();
    void updatePriorities();
    void appendTaskToQueue(const int priority, const TPZAutoPointer<std::packaged_task<void(void) >> task, const bool system_task);
    void checkForMaxAndMinPriority(const int priority);
    ~TPZThreadPool();
    
    /// submit and process a "maximum priority" system task
    // the only use of this method is to run a thread_join task when the number of threads was decreased
    template<typename... Args>
    std::future<void> runSystemTask(const int priority, std::function<void(Args...) > func, Args... args) {
        checkForMaxAndMinPriority(priority);
        TPZAutoPointer < std::packaged_task<void(void) >> task = new std::packaged_task<void(void) >(std::bind(func, args...));
        appendTaskToQueue(priority, task, true);
        return task->get_future();
    }
    
    /// vector of thread objects
    std::vector<std::thread> mThreads;
    /// one mutex to sincronize access to the data structures
    std::mutex mThreadsMutex;
    unsigned int mThreadsToDelete;
    unsigned int mZombieThreads;
    bool mStop;
    

};

#endif // PZTHREADPOOL_H
