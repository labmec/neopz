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
#include "pzadmchunk.h"

template <class TMEM>
class TPZBndCondWithMem : public TPZBndCond {
    
public:
    
    TPZBndCondWithMem() : TPZBndCond(), fMemory(new TPZAdmChunkVector<TMEM>()), fDefaultMem(), fUpdateMem(false){

    }
    

    TPZBndCondWithMem(int matid) : TPZBndCond(matid), fMemory(new TPZAdmChunkVector<TMEM>()), fDefaultMem(), fUpdateMem(false){

    }
    
    /** @brief Default destructor */
    ~TPZBndCondWithMem(){

    }
    
    TPZBndCondWithMem(TPZMaterial * material,int matid,int type,TPZFMatrix<STATE> &val1,TPZFMatrix<STATE> &val2) : TPZBndCond(material,matid,type,val1,val2),fMemory(new TPZAdmChunkVector<TMEM>()), fDefaultMem(), fUpdateMem(false) {
    }
    
    TPZBndCondWithMem(TPZBndCondWithMem<TMEM> & bc) : TPZBndCond(bc), fMemory(new TPZAdmChunkVector<TMEM>()), fDefaultMem(), fUpdateMem(false) {

    }
    
    TPZBndCondWithMem(TPZBndCondWithMem<TMEM> &bc, TPZMaterial * ref) : TPZBndCond(bc,ref), fMemory(new TPZAdmChunkVector<TMEM>()), fDefaultMem(bc.fDefaultMem), fUpdateMem(bc.fUpdateMem) {
         *fMemory = *bc.fMemory;
    }
    
    virtual std::string Name()  override {
        return "TPZBndCondWithMem<TMEM>";
    }
    
    virtual void Print(std::ostream &out) override;
    
    virtual TMEM &MemItem(const int i) const;
    
    /** @brief Unique identifier for serialization purposes */
    int ClassId() const override;
    
    void Write(TPZStream &buf, int withclassid) const override;
    
    void Read(TPZStream &buf, void *context) override;
    
    std::shared_ptr<TPZAdmChunkVector<TMEM>> & GetMemory();
    
    void SetMemory(std::shared_ptr<TPZAdmChunkVector<TMEM>> & memory);
    
    virtual int PushMemItem(int sourceIndex = -1) override;

    virtual void FreeMemItem(int index) override;
    
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
    
    void Clone(std::map<int, TPZMaterial * > &matvec)  override {
        int matid = Id();
        
        TPZMaterial * refmaterial = Material();
        TPZMaterial * newrefmaterial = NULL;
        int refmatid = 0;
        if(refmaterial) {
            refmaterial->Clone(matvec);
            refmatid = refmaterial->Id();
            newrefmaterial = matvec[refmatid];
        }
        std::map<int, TPZMaterial * >::iterator matit;
        matit = matvec.find(matid);
        if(matit == matvec.end())
        {
            TPZMaterial * newmat = (new TPZBndCondWithMem<TMEM>(*this, newrefmaterial));
            matvec[matid] = newmat;
        }
    }
    
    /// Contribute methods

    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override {
        DebugStop();
    }
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override {
        DebugStop();
    }
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override {
        DebugStop();
    }
    
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef) override {
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
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

#endif /* TPZBndCondWithMem_h */
