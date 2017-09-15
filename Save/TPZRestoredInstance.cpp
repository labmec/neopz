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

void TPZRestoredInstance::SetFileVersionInfo(std::map<std::string, long unsigned int> &fileVersionInfo){
    mFileVersionInfo = fileVersionInfo;
}

void TPZRestoredInstance::ReadFromStream(TPZStream &stream, size_t nBytes){
    char temp[nBytes];
    stream.Read(temp, nBytes);
    mStream.Write(temp, nBytes);
}

void TPZRestoredInstance::SetObjId(long unsigned int ObjId) {
    this->mObjId = ObjId;
}

long unsigned int TPZRestoredInstance::GetObjId() const {
    return mObjId;
}

void TPZRestoredInstance::SetClassId(int mClassId) {
    this->mClassId = mClassId;
}

int TPZRestoredInstance::GetClassId() const {
    return mClassId;
}