#include "TPZRestoredInstance.h"
#include <stddef.h>   // for NULL
#include "pzerror.h"  // for DebugStop
class TPZSaveable;
template <class T> class TPZVec;

TPZRestoredInstance::TPZRestoredInstance(){
    me_irl = NULL;
    have_i_had_my_pointers_properly_assigned = false;
}
TPZRestoredInstance::TPZRestoredInstance(TPZSaveable *instance){
    me_irl = instance;
    have_i_had_my_pointers_properly_assigned = false;
}
void TPZRestoredInstance::SetAsRestored() {
    have_i_had_my_pointers_properly_assigned = true;
}
bool TPZRestoredInstance::HaveIBeenRestored() const {
    return have_i_had_my_pointers_properly_assigned;
}

void TPZRestoredInstance::SetPointerToMyObj(TPZSaveable *obj){
    me_irl = obj;
}

TPZSaveable *TPZRestoredInstance::GetPointerToMyObj() const{
    if(me_irl == NULL) DebugStop(); //array entry yet to be initialized.
    return me_irl;
}

TPZVec<int> &TPZRestoredInstance::MyPointersVec() {
    return my_pointers_ids;
}
