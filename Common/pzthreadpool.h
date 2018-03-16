#ifndef PZTHREADPOOL_H
#define PZTHREADPOOL_H

#include <vector>
#include "pz_pthread.h"

#include <queue>
#include <cstdarg>
#include <functional>
#include <thread>

std::priority_queue<PZTask, std::vector<PZTask>, CompareLessByPriority> globalTasksQueue;
pthread_cond_t globalTaskAvailableCond;
pthread_mutex_t globalMutex;

static void threadsLoop() {
    while(true) {
        if(globalTasksQueue.size() == 0)
            pthread_cond_wait(&globalTaskAvailableCond, &globalMutex);
        PZTask task = globalTasksQueue.top();
        task.start();
    }
}

class PZThreadPool
{
public:
    static PZThreadPool &globalInstance() {
        static PZThreadPool globalIntstance;
        return globalIntstance;
    }

    template<typename... Args>
    void run(const int priority, std::function func, Args... args) {
        if(sizeof...(args) > 0) {
            func = std::bind(func, args);
        }
        PZTask<Args...> newTask(priority, func);
        globalTasksQueue.push(newTask);
        pthread_cond_signal(&globalTaskAvailableCond);
    }

    // Helper class for ordering the tasks that the user have requested
    template<typename... Args>
    class PZTask {
    public:
        PZTask(std::function<void(Args...)> func, const int priority) {
            mPriority = priority;
            mFunc = func;
        }

        int priority() const { return mPriority; }
        void start() const { mFunc(); }

    private:
        const int mPriority;
        std::function<void(Args...)> mFunc;
    };

    // Simple struct needed by std::priority_queue for ordering the PZTasks
    struct CompareLessByPriority {
        bool operator()(const PZTask  &lhs, const PZTask  &rhs) {
            return lhs.priority() < rhs.priority();
        }
    };

private:
    std::vector<pthread_t*> mThreads;

    PZThreadPool() {
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

    ~PZThreadPool() {
        for(int i = 0;i < mThreads.size(); i++) {
            pthread_cancel(mThreads[i]);
            delete mThreads[i];
        }
    }
};

#endif // PZTHREADPOOL_H
