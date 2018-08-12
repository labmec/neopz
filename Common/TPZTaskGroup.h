/* 
 * File:   TPZTaskGroup.h
 * Author: thiago
 *
 * Created on 1 de Agosto de 2018, 17:28
 */

#ifndef TPZTASKGROUP_H
#define TPZTASKGROUP_H

#include "tpzautopointer.h"
#include "TPZTask.h"
#include <mutex>
#include "pzstack.h"
#include <set>

class TPZTaskGroup {
public:
    TPZTaskGroup();
    TPZTaskGroup(const TPZTaskGroup& orig);
    void Wait();
    long unsigned int Active();
    virtual ~TPZTaskGroup();
    
    friend class TPZTask;
private:
    void Notify(TPZTask *task);
    void RegisterTask(TPZTask *task);
    std::mutex fMutex;
    std::condition_variable fObservers;
    std::set<TPZTask*> fPendingTasks;
};

#endif /* TPZTASKGROUP_H */

