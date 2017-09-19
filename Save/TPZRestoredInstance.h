#ifndef TPZRESTOREOBJ_H
#define TPZRESTOREOBJ_H

#include <ostream>       // for operator<<
#include "pzmanvector.h" // for TPZManVector
#include "pzvec.h" // for TPZVec
#include "tpzautopointer.h"
#include "TPZSaveable.h"

class TPZSaveable;
class TPZContBufferedStream;

class TPZRestoredInstance {
  public:
    TPZRestoredInstance();
    TPZRestoredInstance(TPZSaveable *);
    void SetInstance(TPZSaveable *);
    TPZSaveable *GetPointerToMyObj() const;
    TPZAutoPointer<TPZSaveable> GetAutoPointerToMyObj();
    TPZVec<int> &MyPointersVec();
    void SetObjId(const long unsigned int &objId);
    long unsigned int GetObjId() const;
    void SetClassId(const int &classId);
    int GetClassId() const;
    virtual ~TPZRestoredInstance();
  protected:
    TPZSaveable *mpInstance;
    TPZManVector<int, 1> mPointersVec;
    TPZContBufferedStream *mpStream;
    TPZAutoPointer<TPZSaveable> mAutoPointerToInstance;
};

#endif // TPZRESTOREOBJ_H
