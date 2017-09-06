#include "TPZRestoredInstance.h"
#include <stddef.h>  // for NULL
#include "pzerror.h" // for DebugStop
class TPZSaveable;
template <class T> class TPZVec;

TPZRestoredInstance::TPZRestoredInstance() {
    mpInstance = NULL;
    mAssignedPointers = false;
}
TPZRestoredInstance::TPZRestoredInstance(TPZSaveable *instance) {
    mpInstance = instance;
    mAssignedPointers = false;
}
void TPZRestoredInstance::SetAsRestored() {
    mAssignedPointers = true;
}
bool TPZRestoredInstance::HaveIBeenRestored() const {
    return mAssignedPointers;
}

void TPZRestoredInstance::SetPointerToMyObj(TPZSaveable *obj) {
    mpInstance = obj;
}

TPZSaveable *TPZRestoredInstance::GetPointerToMyObj() const {
    return mpInstance;
}

TPZVec<int> &TPZRestoredInstance::MyPointersVec() {
#ifdef PZDEBUG
    if (mpInstance == NULL)
        DebugStop(); // vector should be empty if thats the case
#endif
    return mPointersVec;
}
