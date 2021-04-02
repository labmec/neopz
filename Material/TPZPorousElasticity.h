
//
//  TPZPorousElasticity.h
//  pz
//
//  Created by Omar Durán on 4/13/18.
//

#ifndef TPZPorousElasticity_h
#define TPZPorousElasticity_h

#include <stdio.h>
#include "TPZMaterial.h"
#include <iostream>


class TPZPorousElasticity : public TPZMaterial {
    

public:
    
    virtual int ClassId() const override;
    
    TPZPorousElasticity();
    
    TPZPorousElasticity(int matid);
    
    TPZPorousElasticity &operator=(const TPZPorousElasticity &copy);
    
    virtual ~TPZPorousElasticity();
    
    /** @brief Copy constructor */
    TPZPorousElasticity(const TPZPorousElasticity &cp){
        DebugStop();
    }
    
    
    void Print(std::ostream & out) override {
        DebugStop();
    }
    
    std::string Name()  override { return "TPZPorousElasticity"; }
    
    int Dimension() const  override {return 2;}
    
    virtual int NStateVariables() const override{
        DebugStop();
        return 0;
    }
    
    void SetParameters()
    {
        DebugStop();
    }
    
    void FillDataRequirements(TPZMaterialData &data) override {
        DebugStop();
    }
    
    void FillBoundaryConditionDataRequirement (int type, TPZMaterialData &data) override {
        DebugStop();
    }
    

    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override {
        DebugStop();
    }
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override {
        DebugStop();
    }
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
        DebugStop();
    }
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
        DebugStop();
    }
    
    int VariableIndex(const std::string &name)  override {
        DebugStop();
        return 0;
    }
    
    int NSolutionVariables(int var) override {
        DebugStop();
        return 0;
    }
    
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override {
        DebugStop();
    }
    
    /**
     * Save the element data to a stream
     */
    void Write(TPZStream &buf, int withclassid) const override {
        DebugStop();
    }
    
    /**
     * Read the element data from a stream
     */
    void Read(TPZStream &buf, void *context) override {
        DebugStop();
    }
    
};


#endif /* TPZPorousElasticity_h */
