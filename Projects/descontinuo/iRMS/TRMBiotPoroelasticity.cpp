//
//  TRMBiotPoroelasticity.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMBiotPoroelasticity.h"


/** @brief Default constructor */
TRMBiotPoroelasticity::TRMBiotPoroelasticity() : TPZMatWithMem<TRMMemory, TPZMaterial>()
{
    fdimension = 0;
}

/** @brief Constructor based on a material id */
TRMBiotPoroelasticity::TRMBiotPoroelasticity(int matid, int dimension) : TPZMatWithMem<TRMMemory, TPZMaterial>(matid)
{
    fdimension = dimension;
}

/** @brief Default desconstructor */
TRMBiotPoroelasticity::~TRMBiotPoroelasticity(){
    
}

/** @brief Copy constructor $ */
TRMBiotPoroelasticity::TRMBiotPoroelasticity(const TRMBiotPoroelasticity& other){
    this->fdimension    = other.fdimension;
    this->fSimulationData    = other.fSimulationData;
}

/** @brief Copy assignemnt operator $ */
TRMBiotPoroelasticity& TRMBiotPoroelasticity::operator = (const TRMBiotPoroelasticity& other){
    
    if (this != & other) // prevent self-assignment
    {
        this->fdimension    = other.fdimension;
        this->fSimulationData    = other.fSimulationData;
    }
    return *this;
}


void TRMBiotPoroelasticity::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
    
}

void TRMBiotPoroelasticity::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TRMBiotPoroelasticity::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

int TRMBiotPoroelasticity::VariableIndex(const std::string &name) {
    if (!strcmp("u", name.c_str())) return 0;
    if (!strcmp("div_u", name.c_str())) return 1;
    if (!strcmp("s", name.c_str())) return 2;
    return TPZMatWithMem::VariableIndex(name);
}

int TRMBiotPoroelasticity::NSolutionVariables(int var) {
    switch(var) {
        case 0:
            return fdimension; // vector
        case 1:
            return 1; // Scalar
        case 2:
            return fdimension*fdimension; // Tensor
    }
    return TPZMatWithMem::NSolutionVariables(var);
}


void TRMBiotPoroelasticity::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    DebugStop();
}

void TRMBiotPoroelasticity::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    DebugStop();
}

void TRMBiotPoroelasticity::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    DebugStop();
}


void TRMBiotPoroelasticity::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    DebugStop();
}

void TRMBiotPoroelasticity::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    DebugStop();
}

void TRMBiotPoroelasticity::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    DebugStop();
}

void TRMBiotPoroelasticity::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
    DebugStop();
}

void TRMBiotPoroelasticity::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    DebugStop();
}


int TRMBiotPoroelasticity::ClassId() const{
    return -63786378;
}

// -------------------------------------------------------------------------------------------

void TRMBiotPoroelasticity::Write(TPZStream &buf, int withclassid) const{
    
    TPZMaterial::Write(buf, withclassid);
    
}

// -------------------------------------------------------------------------------------------

void TRMBiotPoroelasticity::Read(TPZStream &buf, void *context) {
    TPZMaterial::Read(buf, context);
    
}

// Update element memory by copying the n+1 data to the n data
void TRMBiotPoroelasticity::UpdateMemory()
{
    DebugStop();
}
