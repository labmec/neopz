#include "TPZSaveable.h"
#include "pzmanvector.h"

class TPZRestoreObj {
  public:
    TPZRestoreObj();
    TPZRestoreObj(TPZSaveable *);
    void SetAsRestored();
    bool HaveIBeenRestored() const;
    void SetPointerToMyObj(TPZSaveable *);
    TPZSaveable *GetPointerToMyObj() const;
    TPZVec<int> &MyPointersVec();

  protected:
    TPZSaveable *me_irl;
    bool have_i_had_my_pointers_properly_assigned;
    TPZManVector<int, 1> my_pointers_ids;
    
};
