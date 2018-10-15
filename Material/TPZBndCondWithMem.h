//
//  TPZBndCondWithMem.h
//  pz
//
//  Created by Omar Durán on 10/13/18.
//

#ifndef TPZBndCondWithMem_h
#define TPZBndCondWithMem_h

#include <stdio.h>
#include <memory>
#include "TPZMaterial.h"
#include "pzbndcond.h"

template <class TMEM>
class TPZBndCondWithMem : public TPZBndCond {
    
public:
    
    TPZBndCondWithMem() : TPZBndCond(){
        DebugStop();
    }
    

    TPZBndCondWithMem(int matid) : TPZBndCond(matid){
        DebugStop();
    }
    
    /** @brief Default destructor */
    ~TPZBndCondWithMem(){
        DebugStop();
    }
    
    TPZBndCondWithMem(TPZMaterial * material,int matid,int type,TPZFMatrix<STATE> &val1,TPZFMatrix<STATE> &val2) : TPZBndCond(material,matid,type,val1,val2),fMemory(new TPZAdmChunkVector<TMEM>()), fDefaultMem(), fUpdateMem(0) {
    }
    
    TPZBndCondWithMem(TPZBndCondWithMem<TMEM> & bc) : TPZBndCond(bc), fMemory(new TPZAdmChunkVector<TMEM>()), fDefaultMem(), fUpdateMem(0) {

    }
    
    TPZBndCondWithMem(TPZBndCondWithMem<TMEM> &bc, TPZMaterial * ref) : TPZBndCond(bc,ref), fMemory(new TPZAdmChunkVector<TMEM>()), fDefaultMem(bc.fDefaultMem), fUpdateMem(bc.fUpdateMem) {
         *fMemory = *bc.fMemory;
    }
    
    virtual std::string Name() {
        return "TPZBndCondWithMem<TMEM>";
    }
    
    virtual void Print(std::ostream &out);
    
    virtual TMEM &MemItem(const int i) const;
    
    /** @brief Unique identifier for serialization purposes */
    virtual int ClassId() const;
    
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    virtual void Read(TPZStream &buf, void *context);
    
    std::shared_ptr<TPZAdmChunkVector<TMEM>> & GetMemory();
    
    void SetMemory(std::shared_ptr<TPZAdmChunkVector<TMEM>> & memory);
    
    virtual int PushMemItem(int sourceIndex = -1);

    virtual void FreeMemItem(int index);
    
    void ResetMemItem(int index) {
        this->MemItem(index) = fDefaultMem;
    }
    
    void ResetMemory() {
        int nmem = fMemory->NElements();
        for (unsigned int im = 0; im < nmem; im++) {
            ResetMemItem(im);
        }
    }
    
    virtual void SetDefaultMem(TMEM & defaultMem);
    
    virtual void SetUpdateMem(bool update = 1);
    
    
    void UpdateBCValues(TPZVec<TPZMaterialData> &datavec){
        DebugStop();
    }
    
    /// Contribute methods

    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
protected:
    
    /** @brief Shared pointer to material memory */
    std::shared_ptr<TPZAdmChunkVector<TMEM>> fMemory;
    
    /** @brief Default memory settings */
    TMEM fDefaultMem;
    
    /** @brief Flag to indicate whether the memory data are to be updated in an assemble loop */
    bool fUpdateMem;
    
};

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

#endif /* TPZBndCondWithMem_h */
