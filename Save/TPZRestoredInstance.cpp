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
void TPZRestoredInstance::SetPtrsAsRestored() {
    mAssignedPointers = true;
}
bool TPZRestoredInstance::HaveMyPtrsBeenAssigned() const {
    return mAssignedPointers;
}

void TPZRestoredInstance::SetInstance(TPZSaveable *obj) {
    mpInstance = obj;
}

TPZSaveable *TPZRestoredInstance::GetPointerToMyObj() const {
    return mpInstance;
}

TPZVec<int> &TPZRestoredInstance::MyPointersVec() {
    return mPointersVec;
}
