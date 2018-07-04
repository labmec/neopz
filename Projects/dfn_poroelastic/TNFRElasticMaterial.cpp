//
//  TNFRElasticMaterial.cpp
//  pz
//
//  Created by Omar Durán on 7/2/18.
//

#include "TNFRElasticMaterial.h"


/// Constructor
TNFRElasticMaterial::TNFRElasticMaterial(){
    
}

/// Constructor based on material identifier
TNFRElasticMaterial::TNFRElasticMaterial(int id) : TPZMatWithMem<TNRFElasticMemory,TPZMaterial>(id) {
    
}

/// Destructor
TNFRElasticMaterial::~TNFRElasticMaterial(){
    
}

/// material dimension
int TNFRElasticMaterial::Dimension() const {
    return 2;
}

int TNFRElasticMaterial::NStateVariables(){
    return Dimension();
}

/// setting data on domain
void TNFRElasticMaterial::FillDataRequirements(TPZMaterialData &data){
    DebugStop();
}

///  setting data on boundary
void TNFRElasticMaterial::FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data){
    DebugStop();
}

/// jacobian contribution
void TNFRElasticMaterial::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    DebugStop();
}

/// residual contribution
void TNFRElasticMaterial::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
    DebugStop();
}

void TNFRElasticMaterial::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    DebugStop();
}

/// post-processing
void TNFRElasticMaterial::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    DebugStop();
}

