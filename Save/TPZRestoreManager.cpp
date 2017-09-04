#include "TPZRestoreManager.h"

TPZSaveable *
TPZRestoreManager::AssignPointers(const int &id) {
    TPZRestoreObj &obj = objVec[id];
    if (!obj.HaveIBeenRestored()) {
        obj.SetAsRestored();
        obj.GetPointerToMyObj()->AssignPointers(obj.MyPointersVec());
    }
    
    return obj.GetPointerToMyObj();
}

void TPZRestoreManager::AddObject(TPZSaveable *obj, const int &id){
    const int oldSize = objVec.size();
#ifdef PZDEBUG
    if (oldSize!= id) {
        DebugStop();//objects should be added at the same order as they were saved in the first place.
    }
#endif
    objVec.Resize(oldSize + 1);
    objVec[oldSize].SetPointerToMyObj(obj);    
}
