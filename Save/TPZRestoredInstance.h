#ifndef TPZRESTOREOBJ_H
#define TPZRESTOREOBJ_H

#include <ostream>       // for operator<<
#include "pzmanvector.h" // for TPZManVector
#include "pzvec.h" // for TPZVec
#include "tpzautopointer.h"
#include <memory>

class TPZSavable;
class TPZContBufferedStream;

class TPZRestoredInstance {
  public:
    TPZRestoredInstance();
    TPZRestoredInstance(TPZSavable *);
    virtual ~TPZRestoredInstance();
    void SetInstance(TPZSavable *);
    TPZSavable *GetPointerToMyObj() const;
    TPZAutoPointer<TPZSavable> GetAutoPointerToMyObj();
    std::shared_ptr<TPZSavable> GetSharedPtrToMyObj();
    TPZVec<int> &MyPointersVec();
    void SetObjId(const uint64_t &objId);
    uint64_t GetObjId() const;
    void SetClassId(const int &classId);
    int GetClassId() const;
    void ResetReadStatus();
    bool IsAlreadyRead() const;
    void SetRead();
  protected:
    TPZSavable *mpInstance;
    TPZManVector<int, 1> mPointersVec;
    TPZAutoPointer<TPZSavable> mAutoPointerToInstance;
    std::shared_ptr<TPZSavable> mSharedPtrToInstance;
    bool is_already_read;
};

#endif // TPZRESTOREOBJ_H
