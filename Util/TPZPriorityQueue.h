#ifndef TPZPRIORITYQUEUE_H
#define TPZPRIORITYQUEUE_H

#include <queue>
#include <mutex>
#include <iostream>
#include <algorithm>

template <class T, class Container = std::vector<T>, class Compare = std::less<typename Container::value_type>>
class TPZPriorityQueue : public std::priority_queue<T, Container, Compare> {
public:

    TPZPriorityQueue() : std::priority_queue<T, Container, Compare>(),
    mMutex() {

    }

    void addItem(const T &item) {
        this->push(item);
    }
    
    T popTop() {
        T highestItem = this->top();
        this->pop();
        return highestItem;
    }

    bool remove(T& value) {
        auto it = std::find(this->c.begin(), this->c.end(), value);
        if (it != this->c.end()) {
            this->c.erase(it);
            std::make_heap(this->c.begin(), this->c.end(), this->comp);
            return true;
        } else {
            return false;
        }
    }
    
    std::mutex mMutex;

private:

};

#endif // TPZPRIORITYQUEUE_H
