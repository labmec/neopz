//
//  TRMSegregatedAnalysis.cpp
//  PZ
//
//  Created by Omar on 7/18/16.
//
//

#include "TRMSegregatedAnalysis.h"


TRMSegregatedAnalysis::TRMSegregatedAnalysis() : TPZAnalysis() {
    
    /** @brief define the simulation data */
    fSimulationData = NULL;
    
    /** @brief define the transfer matrices */
    fTransfer = NULL;
    
    /** @brief define the parabolic system */
    fParabolic = NULL;
    
    /** @brief define the hyperbolic system */
    fHyperbolic = NULL;
    
    /** @brief Residue error for flux - pressure */
    ferror_flux_pressure = 0.0;
    
    /** @brief Residue error for saturations */
    ferror_saturation = 0.0;
    
    /** @brief Correction variation for flux - pressure */
    fdx_norm_flux_pressure = 0.0;
    
    /** @brief Correction variation for saturations */
    fdx_norm_saturation = 0.0;
    
}

TRMSegregatedAnalysis::~TRMSegregatedAnalysis(){
    
}

/** @brief Copy constructor $ */
TRMSegregatedAnalysis::TRMSegregatedAnalysis(const TRMSegregatedAnalysis &copy)
{
    fSimulationData         = copy.fSimulationData;
    fTransfer               = copy.fTransfer;
    fParabolic              = copy.fParabolic;
    fHyperbolic             = copy.fHyperbolic;
    ferror_flux_pressure    = copy.ferror_flux_pressure;
    ferror_saturation       = copy.ferror_saturation;
    fdx_norm_flux_pressure  = copy.fdx_norm_flux_pressure;
    fdx_norm_saturation     = copy.fdx_norm_saturation;
    
}

/** @brief Copy assignemnt operator $ */
TRMSegregatedAnalysis & TRMSegregatedAnalysis::operator=(const TRMSegregatedAnalysis &other)
{
    if (this != & other) {  // prevent self-assignment
        
        fSimulationData         = other.fSimulationData;
        fTransfer               = other.fTransfer;
        fParabolic              = other.fParabolic;
        fHyperbolic             = other.fHyperbolic;
        ferror_flux_pressure    = other.ferror_flux_pressure;
        ferror_saturation       = other.ferror_saturation;
        fdx_norm_flux_pressure  = other.fdx_norm_flux_pressure;
        fdx_norm_saturation     = other.fdx_norm_saturation;
    }
    return *this;
}

/** @brief Resize and fill residue and solution vectors */
void TRMSegregatedAnalysis::AdjustVectors(){
    
    fParabolic->AdjustVectors();
    fHyperbolic->AdjustVectors();
}

void TRMSegregatedAnalysis::SegregatedIteration(){

    fParabolic->ExcecuteOneStep();
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    this->UpdateFluxes_at_n();
    
}

void TRMSegregatedAnalysis::ExcecuteOneStep(bool IsFrozenQ){

    if (IsFrozenQ) {
        this->SegregatedIteration();
    }
    
    this->UpdateMemory_at_n();
    fHyperbolic->ExcecuteOneStep();

}

/** @brief Update memory using the Transfer object at state n */
void TRMSegregatedAnalysis::UpdateMemory_at_n(){
    
    fTransfer->s_To_Transport_Memory(fHyperbolic->Meshvec()[0], fHyperbolic->Mesh(),0);// sa
//    fTransfer->s_To_Transport_Memory(fHyperbolic->Meshvec()[0], fSpaceGenerator->TransportMesh().operator->(),1);// sb

    fTransfer->Reciprocal_Memory_Transfer(fParabolic->Mesh(), fHyperbolic->Mesh());
    if (fSimulationData->IsThreePhaseQ()) {
        DebugStop(); 
    }
    
}

/** @brief Update memory using the Transfer object at state n */
void TRMSegregatedAnalysis::UpdateFluxes_at_n(){

    fParabolic->UpdateMemory_at_n();
    fTransfer->un_To_Transport_Mesh(fParabolic->Mesh(), fHyperbolic->Mesh(),true);
    fTransfer->un_To_Transport_Mesh(fParabolic->Mesh(), fHyperbolic->Mesh(),false);
    
    if (fSimulationData->IsThreePhaseQ()) {
        DebugStop();
    }
    
}

/** @brief Update memory using the Transfer object */
void TRMSegregatedAnalysis::UpdateMemory(){

    DebugStop();
//    Mesh()->LoadSolution(fX);
//    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, Mesh());
//    
//    // Volumetric update
//    fTransfer->u_To_Mixed_Memory(fmeshvec[0], Mesh());
//    fTransfer->p_To_Mixed_Memory(fmeshvec[1], Mesh());
    
}

void TRMSegregatedAnalysis::PostProcessStep(bool IsFrozenQ){
    
    if(IsFrozenQ){
        fParabolic->PostProcessStep();
    }
    
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    fHyperbolic->PostProcessStep();
    
}