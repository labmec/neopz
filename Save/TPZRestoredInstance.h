#ifndef TPZRESTOREOBJ_H
#define TPZRESTOREOBJ_H

#include <ostream>       // for operator<<
#include "pzmanvector.h" // for TPZManVector
#include "TPZContBufferedStream.h" 
class TPZSaveable;
template <class T> class TPZVec;

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
  protected:
    TPZSaveable *mpInstance;
    TPZManVector<int, 1> mPointersVec;
    TPZContBufferedStream mStream;
    TPZAutoPointer<TPZSaveable> mAutoPointerToInstance;
};

#endif // TPZRESTOREOBJ_H
