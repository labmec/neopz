
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
    
    virtual int ClassId() const;
    
    TPZPorousElasticity();
    
    TPZPorousElasticity(int matid);
    
    TPZPorousElasticity &operator=(const TPZPorousElasticity &copy);
    
    virtual ~TPZPorousElasticity();
    
    /** @brief Copy constructor */
    TPZPorousElasticity(const TPZPorousElasticity &cp){
        DebugStop();
    }
    
    
    void Print(std::ostream & out){
        DebugStop();
    }
    
    std::string Name() { return "TPZPorousElasticity"; }
    
    int Dimension() const {return 2;}
    
    int NStateVariables(){
        DebugStop();
        return 0;
    }
    
    void SetParameters()
    {
        DebugStop();
    }
    
    void FillDataRequirements(TPZMaterialData &data){
        DebugStop();
    }
    
    void FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data){
        DebugStop();
    }
    

    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    int VariableIndex(const std::string &name){
        DebugStop();
        return 0;
    }
    
    int NSolutionVariables(int var){
        DebugStop();
        return 0;
    }
    
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
        DebugStop();
    }
    
    void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right) {
        DebugStop();
    }
    
    /**
     * Save the element data to a stream
     */
    void Write(TPZStream &buf, int withclassid) const{
        DebugStop();
    }
    
    /**
     * Read the element data from a stream
     */
    void Read(TPZStream &buf, void *context){
        DebugStop();
    }
    
};


#endif /* TPZPorousElasticity_h */
