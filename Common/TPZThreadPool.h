#ifndef TPZTHREADPOOL_H
#define TPZTHREADPOOL_H

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
#include <memory>
#include "TPZTask.h"
#include "TPZReschedulableTask.h"
#include "TPZTaskGroup.h"

class TPZThreadPool {
public:
    static TPZThreadPool &globalInstance();

    std::shared_future<void> run(const int priority, TPZAutoPointer<std::packaged_task<void(void) > > &task, TPZTaskGroup *taskGroup = NULL);

    void run(TPZAutoPointer<TPZReschedulableTask> &task);

    void reschedule(const int priority, TPZAutoPointer<TPZReschedulableTask> &task);
    std::shared_future<void> runNow(TPZAutoPointer<TPZReschedulableTask> &task);

    template<typename... Args>
    std::shared_future<void> run(const int priority, TPZTaskGroup *taskGroup, std::function<void(Args...) > func, Args... args) {
        checkForMaxAndMinPriority(priority);
        TPZAutoPointer < std::packaged_task<void(void) >> task(new std::packaged_task<void(void) >(std::bind(func, args...)));
        return run(priority, task, taskGroup);
    }

    template<typename... Args>
    TPZAutoPointer<TPZReschedulableTask> runReschedulable(const int priority, TPZTaskGroup *taskGroup, std::function<void(Args...) > func, Args... args) {
        checkForMaxAndMinPriority(priority);
        TPZAutoPointer < std::packaged_task<void(void) >> task(new std::packaged_task<void(void) >(std::bind(func, args...)));
        TPZAutoPointer<TPZReschedulableTask> rescTask(new TPZReschedulableTask(priority, task, taskGroup));
        run(rescTask);
        return rescTask;
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
    void appendTaskToQueue(TPZAutoPointer<TPZTask> &task);
    TPZAutoPointer<TPZTask> appendTaskToQueue(const int priority, TPZAutoPointer<std::packaged_task<void(void) >> &task, const bool system_task, TPZTaskGroup *taskGroup = NULL);
    void checkForMaxAndMinPriority(const int priority);
    ~TPZThreadPool();

    template<typename... Args>
    std::shared_future<void> runSystemTask(const int priority, std::function<void(Args...) > func, Args... args) {
        checkForMaxAndMinPriority(priority);
        TPZAutoPointer < std::packaged_task<void(void) >> task(new std::packaged_task<void(void) >(std::bind(func, args...)));
        std::shared_future<void> fut = task->get_future().share();
        appendTaskToQueue(priority, task, true, NULL);
        return fut;
    }


    std::vector<std::thread> mThreads;
    std::mutex mThreadsMutex;
    unsigned int mThreadsToDelete;
    unsigned int mZombieThreads;
    bool mStop;
    TPZPriorityQueue<TPZAutoPointer<TPZTask>, std::vector<TPZAutoPointer<TPZTask>>, TPZTaskOrdering> mTasksQueue;
    std::condition_variable mTaskAvailableCond;
    int mMinPriority;
    int mMaxPriority;

};

#endif // TPZTHREADPOOL_H
