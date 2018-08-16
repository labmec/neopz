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
    
    void addItem(const T &item) {
        this->push(item);
        if (this->c.size() > this->limit){
            this->c.pop_back();
            std::make_heap(this->c.begin(), this->c.end(), this->comp);
        }
    }
    
    virtual ~TPZLimitedPriorityQueue() = default;
private:
    const typename std::vector<T>::size_type limit;
};

#endif /* TPZLIMITEDPRIORITYQUEUE_H */

