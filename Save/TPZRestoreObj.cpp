#include "TPZRestoreObj.h"

TPZRestoreObj::TPZRestoreObj(){
    me_irl = NULL;
    have_i_had_my_pointers_properly_assigned = false;
}
TPZRestoreObj::TPZRestoreObj(TPZSaveable *instance){
    me_irl = instance;
    have_i_had_my_pointers_properly_assigned = false;
}
void TPZRestoreObj::SetAsRestored() {
    have_i_had_my_pointers_properly_assigned = true;
}
bool TPZRestoreObj::HaveIBeenRestored() const {
    return have_i_had_my_pointers_properly_assigned;
}

void TPZRestoreObj::SetPointerToMyObj(TPZSaveable *obj){
    me_irl = obj;
}

TPZSaveable *TPZRestoreObj::GetPointerToMyObj() const{
    if(me_irl == NULL) DebugStop(); //array entry yet to be initialized.
    return me_irl;
}

TPZVec<int> &TPZRestoreObj::MyPointersVec() {
    return my_pointers_ids;
}
