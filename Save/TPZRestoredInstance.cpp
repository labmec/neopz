#include "TPZRestoredInstance.h"
#include <stddef.h>  // for NULL
#include "pzerror.h" // for DebugStop
class TPZSaveable;
template <class T> class TPZVec;

TPZRestoredInstance::TPZRestoredInstance() {
    mpInstance = NULL;
}
TPZRestoredInstance::TPZRestoredInstance(TPZSaveable *instance) {
    mpInstance = instance;
}

void TPZRestoredInstance::SetInstance(TPZSaveable *obj) {
    mpInstance = obj;
}

TPZSaveable *TPZRestoredInstance::GetPointerToMyObj() const {
    return mpInstance;
}

TPZAutoPointer<TPZSaveable> TPZRestoredInstance::GetAutoPointerToMyObj(){
    if (!mAutoPointerToInstance){
        mAutoPointerToInstance = TPZAutoPointer<TPZSaveable>(mpInstance);
    }
    return mAutoPointerToInstance;
}

TPZVec<int> &TPZRestoredInstance::MyPointersVec() {
    return mPointersVec;
}