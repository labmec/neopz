#ifndef PZPRIORITYQUEUE_H
#define PZPRIORITYQUEUE_H

#include <queue>
#include <pthread.h>
#include <iostream>
template <class T,  class Container = std::vector<T>, class Compare = std::less<typename Container::value_type>>
class TPZPriorityQueue
{
public:
    TPZPriorityQueue() {
        pthread_mutex_init(&mMutex, NULL);
    }

    void addItem(const T &item) {
        pthread_mutex_lock(&mMutex);
        mQueue.push(item);
        pthread_mutex_unlock(&mMutex);
    }

    int size() const { return mQueue.size(); }

    friend class TPZThreadPool;

private:
    pthread_mutex_t mMutex;
    std::priority_queue<T, Container, Compare> mQueue;

    T popTop() {
        T highestItem = mQueue.top();
        mQueue.pop();
        return highestItem;
    }

};

#endif // PZPRIORITYQUEUE_H
