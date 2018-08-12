/* 
 * File:   TPZTask.cpp
 * Author: Thiago
 * 
 * Created on 19 de Junho de 2018, 14:08
 */

#include "TPZTask.h"

#include "TPZTaskGroup.h"

TPZTask::TPZTask(const int priority, TPZAutoPointer<std::packaged_task<void(void)>> &task, TPZTaskGroup *taskGroup) :
mSystemTask(false),
mPriority(priority),
mTask(task),
mState(EProcessingState::CREATED),
mTaskGroup(taskGroup) {
    if (taskGroup){
        taskGroup->RegisterTask(this);
    }
}

int TPZTask::priority() const {
    return mPriority;
}

void TPZTask::start() {
    mState = EProcessingState::STARTED;
    (*mTask)();
    mState = EProcessingState::FINISHED;
    if (mTaskGroup) {
        mTaskGroup->Notify(this);
    }
}

void TPZTask::Cancel() {
    if (mTaskGroup){
        mTaskGroup->Notify(this);
    }
}


TPZTask::~TPZTask() {
}