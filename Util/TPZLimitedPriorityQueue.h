/* 
 * File:   TPZLimitedPriorityQueue.h
 * Author: thiago
 *
 * Created on 10 de Agosto de 2018, 18:12
 */

#ifndef TPZLIMITEDPRIORITYQUEUE_H
#define TPZLIMITEDPRIORITYQUEUE_H

#include "TPZPriorityQueue.h"
#include <vector>

template <class T, class Compare = std::less<typename std::vector<T>::value_type>>
class TPZLimitedPriorityQueue : public TPZPriorityQueue<T, std::vector<T>, Compare> {
public:
    TPZLimitedPriorityQueue(const typename std::vector<T>::size_type limit) : TPZPriorityQueue<T, std::vector<T>, Compare>(), limit(limit) {
    }
    TPZLimitedPriorityQueue(const TPZLimitedPriorityQueue& orig) = default;
        
    TPZLimitedPriorityQueue& operator = (const TPZLimitedPriorityQueue& other) {
        TPZPriorityQueue<T, std::vector<T>, Compare>::operator =(other);
        limit = other.limit;
        return *this;
    }
    
    void addItem(const T &item) {
        this->push(item);
        std::sort(this->c.begin(), this->c.end(), this->comp);
        if (this->c.size() > this->limit){
            this->c.pop_back();
        }
    }
    
    virtual ~TPZLimitedPriorityQueue() = default;
private:
    typename std::vector<T>::size_type limit;
};

#endif /* TPZLIMITEDPRIORITYQUEUE_H */

