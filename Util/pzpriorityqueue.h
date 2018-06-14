#ifndef PZPRIORITYQUEUE_H
#define PZPRIORITYQUEUE_H

#include <queue>
#include <mutex>
#include <iostream>
template <class T,  class Container = std::vector<T>, class Compare = std::less<typename Container::value_type>>
class TPZPriorityQueue
{
public:
    TPZPriorityQueue() :
    mMutex(),
    mQueue() {

    }

    void addItem(const T &item) {
        mQueue.push(item);
    }

    T popTop() {
        T highestItem = mQueue.top();
        mQueue.pop();
        return highestItem;
    }

    T top() const {
        T top = mQueue.top();
        return top;
    }

    int size() const { return mQueue.size(); }

    std::mutex mMutex;

private:

    std::priority_queue<T, Container, Compare> mQueue;
};

#endif // PZPRIORITYQUEUE_H
