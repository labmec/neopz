//
//  TPZBndCondWithMem_impl.h
//  pz
//
//  Created by Omar Durán on 6/12/19.
//

#include "TPZBndCondWithMem.h"
template <class TMEM>
void TPZBndCondWithMem<TMEM>::Print(std::ostream &out)
{
    out << this->Name();
    TPZBndCond::Print(out);
}

template <class TMEM>
TMEM & TPZBndCondWithMem<TMEM>::MemItem(const int i) const{
    return fMemory.get()->operator [](i);
}

template <class TMEM>
int TPZBndCondWithMem<TMEM>::ClassId() const{
    return Hash("TPZBndCondWithMem") ^ TPZBndCond::ClassId() << 1 ^ TMEM().ClassId() << 2;
}

template <class TMEM>
void TPZBndCondWithMem<TMEM>::Write(TPZStream &buf, int withclassid) const{
    TPZBndCond::Write(buf, withclassid);
    int updatemem = fUpdateMem;
    buf.Write(&updatemem);
    fDefaultMem.Write(buf, 0);
    TPZPersistenceManager::WritePointer(fMemory.get(), &buf);
}

template <class TMEM>
void TPZBndCondWithMem<TMEM>::Read(TPZStream &buf, void *context){
    TPZBndCond::Read(buf, context);
    int updatemem;
    buf.Read(&updatemem);
    if (updatemem) {
        fUpdateMem = true;
    } else {
        fUpdateMem = false;
    }
    fDefaultMem.Read(buf, 0);
    fMemory = std::dynamic_pointer_cast<TPZAdmChunkVector<TMEM> >(TPZPersistenceManager::GetSharedPointer(&buf));
}

template <class TMEM>
std::shared_ptr<TPZAdmChunkVector<TMEM>> & TPZBndCondWithMem<TMEM>::GetMemory(){
    return fMemory;
}

template <class TMEM>
void TPZBndCondWithMem<TMEM>::SetMemory(std::shared_ptr<TPZAdmChunkVector<TMEM>> & memory){
    fMemory = memory;
}

template <class TMEM>
int TPZBndCondWithMem<TMEM>::PushMemItem(int sourceIndex){
    int index = fMemory->AllocateNewElement();
    if (sourceIndex < 0) {
        this->ResetMemItem(index);
    } else {
        this->MemItem(index) = this->MemItem(sourceIndex);
    }
    return index;
}

template <class TMEM>
void TPZBndCondWithMem<TMEM>::FreeMemItem(int index){
    fMemory->SetFree(index);
}

template <class TMEM>
void TPZBndCondWithMem<TMEM>::SetDefaultMem(TMEM & defaultMem){
    fDefaultMem = defaultMem;
}

template <class TMEM>
void TPZBndCondWithMem<TMEM>::SetUpdateMem(bool update){
    fUpdateMem = update;
}

template <class TMEM>
void TPZBndCondWithMem<TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    TPZBndCond::Contribute(data, weight, ek, ef);
}

template <class TMEM>
void TPZBndCondWithMem<TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> ek_fake(ef.Rows(),ef.Rows());
    this->Contribute(data, weight, ek_fake, ef);
}

template <class TMEM>
void TPZBndCondWithMem<TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    TPZBndCond::Contribute(datavec, weight, ek, ef);
}

template <class TMEM>
void TPZBndCondWithMem<TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> ek_fake(ef.Rows(),ef.Rows());
    this->Contribute(datavec, weight, ek_fake, ef);
}

#define IMPLEMENTTPZBNDCONDWITHMEM(TMEM) \
\
template class \
TPZRestoreClass< TPZBndCondWithMem<TMEM> >;
