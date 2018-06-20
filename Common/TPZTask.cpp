/* 
 * File:   TPZTask.cpp
 * Author: Thiago
 * 
 * Created on 19 de Junho de 2018, 14:08
 */

#include "TPZTask.h"

TPZTask::TPZTask(const int priority, TPZAutoPointer<std::packaged_task<void(void)>> &task) :
mSystemTask(false),
mPriority(priority),
mTask(task) {

}

int TPZTask::priority() const {
    return mPriority;
}

void TPZTask::start() {
    (*mTask)();
}

TPZTask::~TPZTask() {
}