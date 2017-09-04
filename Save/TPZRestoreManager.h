#include "TPZSaveable.h"
#include "pzmanvector.h"
#include "TPZRestoreObj.h"

class TPZRestoreManager {
  public:
    TPZSaveable *AssignPointers(const int &id);
    void AddObject(TPZSaveable *, const int &id);

  protected:
    TPZRestoreManager();
    TPZManVector<TPZRestoreObj, 10> objVec;
};
