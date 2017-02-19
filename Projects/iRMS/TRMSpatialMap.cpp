//
//  TRMSpatialMap.cpp
//  PZ
//
//  Created by Omar on 2/19/17.
//
//

#include "TRMSpatialMap.h"


TRMSpatialMap::TRMSpatialMap() : TPZDiscontinuousGalerkin ()
{
    
    /** @brief material dimension */
    fdimension = 0;
    
}

TRMSpatialMap::TRMSpatialMap(int matid, int dimension) : TPZDiscontinuousGalerkin (matid)
{
    
    /** @brief material dimension */
    fdimension = dimension;
}


TRMSpatialMap::TRMSpatialMap(const TRMSpatialMap &mat) : TPZDiscontinuousGalerkin (mat)
{
    this->fdimension        = mat.fdimension;
}

TRMSpatialMap::~TRMSpatialMap()
{
    
}

void TRMSpatialMap::FillDataRequirements(TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
}


void TRMSpatialMap::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

int TRMSpatialMap::VariableIndex(const std::string &name) {
    if (!strcmp("kappa", name.c_str())) return 0;
    if (!strcmp("phi", name.c_str())) return 1;
    return -1;
}

int TRMSpatialMap::NSolutionVariables(int var) {
    switch(var) {
        case 0:
            return 3; // Vector
        case 1:
            return 1; // Scalar
    }
    return -1;
}

void TRMSpatialMap::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout) {
    
    Solout.Resize(this->NSolutionVariables(var));
    switch(var) {
        case 0:
        {
            Solout[0] = data.sol[0][0];
            Solout[1] = data.sol[0][1];
            Solout[2] = data.sol[0][2];
        }
            break;
        case 1:
        {
            Solout[0] = data.sol[0][3];
        }
            break;
        default:
        {
            DebugStop();
        }
    }
    
}


int TRMSpatialMap::ClassId() const {
    return -637806378114538;
}

// -------------------------------------------------------------------------------------------

void TRMSpatialMap::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    
}

// -------------------------------------------------------------------------------------------

void TRMSpatialMap::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    
}

