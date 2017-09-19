#include "TPZRestoredInstance.h"
#include <stddef.h>  // for NULL
#include "pzerror.h"
#include "TPZContBufferedStream.h" // for DebugStop
class TPZSaveable;
template <class T> class TPZVec;

TPZRestoredInstance::TPZRestoredInstance() {
    mpInstance = NULL;
    mpStream = new TPZContBufferedStream();
}

TPZRestoredInstance::TPZRestoredInstance(TPZSaveable *instance) {
    mpInstance = instance;
    mpStream = new TPZContBufferedStream();
}

TPZRestoredInstance::~TPZRestoredInstance() {
    delete mpStream;
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