#ifndef TPZRESTOREOBJ_H
#define TPZRESTOREOBJ_H

#include <ostream>       // for operator<<
#include "pzmanvector.h"
#include "TPZContBufferedStream.h" // for TPZManVector
class TPZSaveable;
template <class T> class TPZVec;

class TPZRestoredInstance {
  public:
    TPZRestoredInstance();
    TPZRestoredInstance(TPZSaveable *);
    void SetPtrsAsRestored();
    bool HaveMyPtrsBeenAssigned() const;
    void SetInstance(TPZSaveable *);
    TPZSaveable *GetPointerToMyObj() const;
    TPZVec<int> &MyPointersVec();
    void SetFileVersionInfo(std::map<std::string, long unsigned int> &fileVersionInfo);
    void ReadFromStream(TPZStream &stream, size_t nBytes);
    void SetObjId(long unsigned int objId);
    long unsigned int GetObjId() const;
    void SetClassId(int mClassId);
    int GetClassId() const;
  protected:
    TPZSaveable *mpInstance;
    bool mAssignedPointers;
    TPZManVector<int, 1> mPointersVec;
    TPZContBufferedStream mStream;
    std::map<std::string, long unsigned int> mFileVersionInfo;
    
    long unsigned int mObjId;
    int mClassId;
};

#endif // TPZRESTOREOBJ_H
