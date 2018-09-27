#ifndef TPZPRIORITYQUEUE_H
#define TPZPRIORITYQUEUE_H

#include <queue>
#include <mutex>
#include <iostream>
#include <algorithm>

template <class T, class Container = std::vector<T>, class Compare = std::less<typename Container::value_type>>
class TPZPriorityQueue {
public:

    TPZPriorityQueue() {

    }
    
    TPZPriorityQueue(const TPZPriorityQueue &other) {
        this->c = other.c;
        this->comp = other.comp;
    };
    
    TPZPriorityQueue& operator = (const TPZPriorityQueue& other){
        if (this != &other){
            this->c = other.c;
            this->comp = other.comp;
        }
        return *this;
    }

    void addItem(const T &item) {
        this->c.push_back(item);
        std::sort(this->c.begin(), this->c.end(), this->comp);
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
            std::sort(this->c.begin(), this->c.end(), this->comp);
            return true;
        } else {
            return false;
        }
    }

    void remove(const typename Container::size_type begin, const typename Container::size_type end) {
        this->c.erase(this->c.begin() + begin, this->c.begin() + end);
    }
    
    T top() {
        return c.operator[](0);
    }
    
    typename Container::size_type size() const {
        return c.size();
    }
    
    const T &getItem(const typename Container::size_type index) const {
        return this->c.operator[](index);
    }
    
    void pop() {
        c.erase(c.begin());
    }
    
    void pop_back(const typename Container::size_type count) {
        c.erase(c.end()-count, c.end());
    }
    
    void push(T &item) {
        c.push_back(item);
    }
    
    void push(const T &item) {
        c.push_back(item);
    }
    
    mutable std::mutex mMutex;

protected :

    Container c;
    Compare comp;
    
};

#endif // TPZPRIORITYQUEUE_H
