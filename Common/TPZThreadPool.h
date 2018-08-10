#ifndef TPZTHREADPOOL_H
#define TPZTHREADPOOL_H

#include <vector>
#include "TPZPriorityQueue.h"
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

/// Administers tasks that will be executed asynchronously
class TPZThreadPool {
public:
    /// @return a reference to the one and only instance of TPZThreadPool
    static TPZThreadPool &globalInstance();
    
    /// submits a task to be executed by TPZThreadPool
    // @return a std::future<void> object that allows to block the calling thread if necessary (by calling wait())
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

    /// sets the number of threads to be executed simultaneously
    void SetNumThreads(const unsigned numThreads);
    int maxPriority() const;
    int minPriority() const;
    
    /// @return the number of threads currently available
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

    /// Submits and processes a "maximum priority" system task
    // The only use of this method is to run a thread_join task when the number of threads was decreased
    template<typename... Args>
    std::shared_future<void> runSystemTask(const int priority, std::function<void(Args...) > func, Args... args) {
        checkForMaxAndMinPriority(priority);
        TPZAutoPointer < std::packaged_task<void(void) >> task(new std::packaged_task<void(void) >(std::bind(func, args...)));
        std::shared_future<void> fut = task->get_future().share();
        appendTaskToQueue(priority, task, true, NULL);
        return fut;
    }

    /// vector of thread objects
    std::vector<std::thread> mThreads;
    /// one mutex to synchronize access to the data structures
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
