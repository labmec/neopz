/* 
 * File:   TPZTaskGroup.cpp
 * Author: thiago
 * 
 * Created on 1 de Agosto de 2018, 17:28
 */

#include "TPZTaskGroup.h"
#include <mutex>
#include <condition_variable>
		
TPZTaskGroup::TPZTaskGroup() {
}

TPZTaskGroup::TPZTaskGroup(const TPZTaskGroup& orig) {
}

TPZTaskGroup::~TPZTaskGroup() {
}

void TPZTaskGroup::Notify(TPZTask *task) {
    std::unique_lock<std::mutex> lock(fMutex);
    fPendingTasks.erase(task);
    if (fPendingTasks.empty()){
        fObservers.notify_all();
    }
}

void TPZTaskGroup::RegisterTask(TPZTask *task) {
    std::unique_lock<std::mutex> lock(fMutex);
    fPendingTasks.insert(task);
}

void TPZTaskGroup::Wait() {
    std::unique_lock<std::mutex> lock(fMutex);
    fObservers.wait(lock, [this] {
        return fPendingTasks.empty();
    });
}

long unsigned int TPZTaskGroup::Active() {
   return fPendingTasks.size();
}
