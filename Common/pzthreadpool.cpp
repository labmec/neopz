#include "pzthreadpool.h"

TPZPriorityQueue<TPZTask*, std::vector<TPZTask*>, TPZOrderGreaterToMin> globalTasksQueue;
std::condition_variable globalTaskAvailableCond;
int globalMinPriority = std::numeric_limits<int>::max();
int globalMaxPriority = std::numeric_limits<int>::min();

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
                return stop || globalTasksQueue.size() != 0;
            });
            if (stop) {
                return;
            }
            task = globalTasksQueue.popTop();
            updatePriorities();
        }
        if (task) {
            task->start();
        }
    }
}

TPZThreadPool::TPZThreadPool() : stop(false) {
    const unsigned numThreads = std::thread::hardware_concurrency();
    for (unsigned int i = 0; i < numThreads; ++i) {
        mThreads.emplace_back([this]{
            threadsLoop();
        });
    }
}

TPZThreadPool::~TPZThreadPool() {
    {
        std::unique_lock<std::mutex> lock(globalTasksQueue.mMutex);
        stop = true;
    }
    globalTaskAvailableCond.notify_all();
    for (auto &thread: mThreads) {
        thread.join();
    }
    while (globalTasksQueue.size() != 0){
        TPZTask *task = globalTasksQueue.popTop();
        delete task;
    }
}

int TPZThreadPool::maxPriority() const {
    return globalMaxPriority;
}

int TPZThreadPool::minPriority() const {
    return globalMinPriority;
}

int TPZThreadPool::threadCount() const {
    return mThreads.size();
}

void TPZThreadPool::appendTaskToQueue(const int priority, const TPZAutoPointer<std::packaged_task<void ()> > task) {
    TPZTask *newTask = new TPZTask(priority, task);
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
