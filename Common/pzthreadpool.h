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
    mPriority(priority),
    mTask(task) {

    }

    int priority() const {
        return mPriority;
    }

    void start() const {
        (*mTask)();
    }

private:
    const int mPriority;
    TPZAutoPointer<std::packaged_task<void(void) >> mTask;
};

// Simple struct needed by std::priority_queue for ordering the items

struct TPZOrderGreaterToMin {

    bool operator()(const TPZTask *lhs, const TPZTask *rhs) {
        return lhs->priority() < rhs->priority();
    }
};

class TPZThreadPool {
public:
    static TPZThreadPool &globalInstance();

    template<typename... Args>
    std::future<void> run(const int priority, std::function<void(Args...) > func, Args... args) {
        checkForMaxAndMinPriority(priority);
        TPZAutoPointer < std::packaged_task<void(void) >> task = new std::packaged_task<void(void) >(std::bind(func, args...));
        appendTaskToQueue(priority, task);
        return task->get_future();
    }

    int maxPriority() const;
    int minPriority() const;
    int threadCount() const;

private:
    TPZThreadPool();
    void threadsLoop();
    void updatePriorities();
    void appendTaskToQueue(const int priority, const TPZAutoPointer<std::packaged_task<void(void) >> task);
    void checkForMaxAndMinPriority(const int priority);
    ~TPZThreadPool();
    
    std::vector<std::thread> mThreads;
    bool stop;

};

#endif // PZTHREADPOOL_H
