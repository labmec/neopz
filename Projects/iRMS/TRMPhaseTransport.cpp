//
//  TRMPhaseTransport.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMPhaseTransport.h"

TRMPhaseTransport::TRMPhaseTransport() : TPZMatWithMem<TRMPhaseMemory, TPZDiscontinuousGalerkin>()
{
    
}

TRMPhaseTransport::TRMPhaseTransport(int matid) : TPZMatWithMem<TRMPhaseMemory, TPZDiscontinuousGalerkin>(matid)
{
    
}


TRMPhaseTransport::TRMPhaseTransport(const TRMPhaseTransport &mat) : TPZMatWithMem<TRMPhaseMemory, TPZDiscontinuousGalerkin>(mat)
{
    
}

TRMPhaseTransport::~TRMPhaseTransport()
{
    
}

void TRMPhaseTransport::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

void TRMPhaseTransport::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TRMPhaseTransport::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

int TRMPhaseTransport::VariableIndex(const std::string &name) {
    if (!strcmp("p", name.c_str())) return 0;
    if (!strcmp("u", name.c_str())) return 1;
    if (!strcmp("div_u", name.c_str())) return 2;
    if (!strcmp("s_a", name.c_str())) return 3;
    //    if (!strcmp("AWeightedPressure", name.c_str())) return 3;
    //    if (!strcmp("ABulkVelocity", name.c_str())) return 4;
    //    if (!strcmp("ADivOfBulkVeclocity", name.c_str())) return 5;
    return TPZMatWithMem::VariableIndex(name);
}

int TRMPhaseTransport::NSolutionVariables(int var) {
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return 3; // Vector
        case 2:
            return 1; // Scalar
        case 3:
            return 1; // Scalar
    }
    return TPZMatWithMem::NSolutionVariables(var);
}

void TRMPhaseTransport::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    DebugStop();
//    switch (fSimulationData->SystemType().size()) {
//        case 1:
//        {
//            Solution_a(datavec, var, Solout);
//        }
//            break;
//        case 2:
//        {
//            Solution_ab(datavec, var, Solout);
//        }
//            break;
//        case 3:
//        {
//            DebugStop();
//        }
//            break;
//        default:
//        {
//            DebugStop();
//        }
//            break;
//    }
    
}

// Jacobian contribution
void TRMPhaseTransport::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    DebugStop();
//    switch (fSimulationData->SystemType().size()) {
//        case 1:
//        {
//            Contribute_a(datavec, weight, ek, ef);
//        }
//            break;
//        case 2:
//        {
//            Contribute_ab(datavec, weight, ek, ef);
//        }
//            break;
//        case 3:
//        {
//            DebugStop();
//        }
//            break;
//        default:
//        {
//            DebugStop();
//        }
//            break;
//    }
    
}


void TRMPhaseTransport::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
//    switch (fSimulationData->SystemType().size()) {
//        case 1:
//        {
//            Contribute_a(datavec, weight, ef);
//        }
//            break;
//        case 2:
//        {
//            Contribute_ab(datavec, weight, ef);
//        }
//            break;
//        case 3:
//        {
//            DebugStop();
//        }
//            break;
//        default:
//        {
//            DebugStop();
//        }
//            break;
//    }
}

//void TRMPhaseTransport::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
//{
//    DebugStop();
//}
//
//void TRMPhaseTransport::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef)
//{
//    std::cout << " This method should be called only for Capillary pressure terms " << std::endl;
//    DebugStop();
//}
//
//
//void TRMPhaseTransport::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
//{
//    DebugStop();
//}
//
//void TRMPhaseTransport::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
//{
//    std::cout << " This method should be called only for Capillary pressure terms " << std::endl;
//    DebugStop();
//}
//
//
//void TRMPhaseTransport::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
//{
//    DebugStop();
//}


int TRMPhaseTransport::ClassId() const {
    return -637806378;
}

// -------------------------------------------------------------------------------------------

void TRMPhaseTransport::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    
}

// -------------------------------------------------------------------------------------------

void TRMPhaseTransport::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    
}

// Update element memory by copying the n+1 data to the n data
void TRMPhaseTransport::UpdateMemory()
{
    DebugStop();
//    long nel = fMemory.NElements();
//    for (long el=0; el<nel; el++) {
//        fMemory[el].UpdateSolutionMemory();
//    }
}

