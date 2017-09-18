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

void TPZRestoredInstance::ReadFromStream(TPZStream &stream, size_t nBytes){
    char temp[nBytes];
    stream.Read(temp, nBytes);
    mStream.Write(temp, nBytes);
}

void TPZRestoredInstance::SetObjId(const long unsigned int &objId) {
    this->mObjId = objId;
}

long unsigned int TPZRestoredInstance::GetObjId() const {
    return mObjId;
}

void TPZRestoredInstance::SetClassId(const int &classId) {
    this->mClassId = classId;
}

int TPZRestoredInstance::GetClassId() const {
    return mClassId;
}
