#ifndef PZTHREADPOOL_H
#define PZTHREADPOOL_H

#include <vector>
#include <pthread.h>
#include "pzpriorityqueue.h"
#include <cstdarg>
#include <functional>
#include <thread>
#include <chrono>
#include <random>

// Helper class for ordering the tasks that the user have requested
class TPZTask {
public:
    TPZTask(const int priority, std::function<void(void)> func) :
    mPriority(priority),
    mFunc(func)
    {

    }

    int priority() const { return mPriority; }
    void start() const { mFunc(); }

private:
    int mPriority;
    std::function<void(void)> mFunc;
};

// Simple struct needed by std::priority_queue for ordering the items
struct TPZOrderGreaterToMin {
    bool operator()(const TPZTask *lhs, const TPZTask *rhs) {
        return lhs->priority() < rhs->priority();
    }
};

TPZPriorityQueue<TPZTask*, std::vector<TPZTask*>, TPZOrderGreaterToMin> globalTasksQueue;
pthread_cond_t globalTaskAvailableCond;
pthread_mutex_t globalMutex;
pthread_mutex_t globalQueueSizeMutex;



class TPZThreadPool
{
public:
    static TPZThreadPool &globalInstance() {
        static TPZThreadPool globalIntstance;
        return globalIntstance;
    }

    void addTask(const int priority, const std::function<void(void)> &func) {
        TPZTask *newTask = new TPZTask(priority, func);
        globalTasksQueue.addItem(newTask);
        pthread_cond_signal(&globalTaskAvailableCond);
    }

    template<typename... Args>
    void run(const int priority, std::function<void(Args...)> func, Args... args) {
        std::function<void(void)> taskFunc = std::bind(func, args...);
        addTask(priority, taskFunc);
    }

private:
    std::vector<pthread_t*> mThreads;

    TPZThreadPool() {
        // setup the thread control structures
        pthread_cond_init(&globalTaskAvailableCond, NULL);
        pthread_mutex_init(&globalMutex, NULL);

        // Finds the optimal number of threads
        unsigned numThreads = std::thread::hardware_concurrency();
        // Creates the threads
        mThreads.resize(numThreads);
        for(int i = 0; i < numThreads; i++) {
            mThreads[i] = new pthread_t();
            pthread_create(mThreads[i], NULL, threadsLoop, NULL);
        }
    }

    static void *threadsLoop(void*) {
        while(true) {
            pthread_mutex_lock(&globalTasksQueue.mMutex);
            if (globalTasksQueue.size() == 0) {
                pthread_cond_wait(&globalTaskAvailableCond, &globalTasksQueue.mMutex);
            }
            TPZTask *task = NULL;
            if (globalTasksQueue.size() != 0) {
                task = globalTasksQueue.popTop();
            }
            pthread_mutex_unlock(&globalTasksQueue.mMutex);
            if (task) {
                task->start();
                delete task;
            }
        }
    }

    ~TPZThreadPool() {
        for(int i = 0;i < mThreads.size(); i++) {
            pthread_cancel(*mThreads[i]);
            delete mThreads[i];
        }
    }
};

#endif // PZTHREADPOOL_H
