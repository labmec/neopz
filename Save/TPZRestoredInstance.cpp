#include "TPZRestoredInstance.h"
#include <stddef.h>  // for NULL
#include "pzerror.h" // for DebugStop
#include "TPZSavable.h"

template <class T> class TPZVec;

TPZRestoredInstance::TPZRestoredInstance() : mpInstance(NULL), is_already_read(false) {
}

TPZRestoredInstance::TPZRestoredInstance(TPZSavable *instance) : mpInstance(instance), is_already_read(false) {
}

TPZRestoredInstance::~TPZRestoredInstance() {

}

void TPZRestoredInstance::SetInstance(TPZSavable *obj) {
    mpInstance = obj;
}

TPZSavable *TPZRestoredInstance::GetPointerToMyObj() const {
    return mpInstance;
}

TPZAutoPointer<TPZSavable> TPZRestoredInstance::GetAutoPointerToMyObj(){
#ifdef PZDEBUG
        if (mSharedPtrToInstance) {
            DebugStop();
        }
#endif
    if (!mAutoPointerToInstance){
        mAutoPointerToInstance = TPZAutoPointer<TPZSavable>(mpInstance);
    }
    return mAutoPointerToInstance;
}

std::shared_ptr<TPZSavable> TPZRestoredInstance::GetSharedPtrToMyObj(){
#ifdef PZDEBUG
        if (mAutoPointerToInstance) {
            DebugStop();
        }
#endif
    if (!mSharedPtrToInstance){
        mSharedPtrToInstance = std::shared_ptr<TPZSavable>(mpInstance);
    }
    return mSharedPtrToInstance;
}

TPZVec<int> &TPZRestoredInstance::MyPointersVec() {
    return mPointersVec;
}

void TPZRestoredInstance::ResetReadStatus() {
    is_already_read = false;
}

bool TPZRestoredInstance::IsAlreadyRead() const {
    return is_already_read;
}

void TPZRestoredInstance::SetRead() {
    is_already_read = true;
}
