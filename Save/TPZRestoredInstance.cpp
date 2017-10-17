#include "TPZRestoredInstance.h"
#include <stddef.h>  // for NULL
#include "pzerror.h" // for DebugStop

class TPZSavable;
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
    if (!mAutoPointerToInstance){
        mAutoPointerToInstance = TPZAutoPointer<TPZSavable>(mpInstance);
    }
    return mAutoPointerToInstance;
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