#include "pzthreadpool.h"

TPZPriorityQueue<TPZTask*, std::vector<TPZTask*>, TPZOrderGreaterToMin> globalTasksQueue;
std::condition_variable globalTaskAvailableCond;
int globalMinPriority = std::numeric_limits<int>::max();
int globalMaxPriority = std::numeric_limits<int>::min();

void updatePriorities()
{
    if(globalTasksQueue.size() != 0) {
        globalMaxPriority = globalTasksQueue.top()->priority();
    } else {
        globalMaxPriority = std::numeric_limits<int>::min();
        globalMinPriority = std::numeric_limits<int>::max();
    }
}

void threadsLoop() {
    while(true) {
		TPZAutoPointer<TPZTask> task;
		{
			std::unique_lock<std::mutex> lock(globalTasksQueue.mMutex);
			if (globalTasksQueue.size() == 0) {
				globalTaskAvailableCond.wait(lock);
			}
			if (globalTasksQueue.size() != 0) {
				task = globalTasksQueue.popTop();
				updatePriorities();
			}
		}
		if (task) {
			task->start();
		}
    }
}

TPZThreadPool::TPZThreadPool() {
    const unsigned numThreads = std::thread::hardware_concurrency();
    mThreads.resize(numThreads);
	for (unsigned int i = 0; i < numThreads; ++i) {
		mThreads[i] = std::thread(threadsLoop);
	}
}

TPZThreadPool::~TPZThreadPool() {
	for (unsigned int i = 0; i < threadCount(); ++i) {
		mThreads[i].join();
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

void TPZThreadPool::appendTaskToQueue(const int priority, const TPZAutoPointer<std::packaged_task<void ()> > task)  {
    TPZTask *newTask = new TPZTask(priority, task);
    globalTasksQueue.addItem(newTask);
    globalTaskAvailableCond.notify_one();
}

void TPZThreadPool::checkForMaxAndMinPriority(const int priority) {
    if(priority > globalMaxPriority) {
        globalMaxPriority = priority;
    }
    if(priority < globalMinPriority){
        globalMinPriority = priority;
    }
}

TPZThreadPool &TPZThreadPool::globalInstance() {
    static TPZThreadPool globalIntstance;
    return globalIntstance;
}
