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

void TRMSegregatedAnalysis::SegregatedIteration(bool IsActiveQ){
    

    this->UpdateMemory();
//    this->UpdateMemory_at_n();    
    if (IsActiveQ) {
        fParabolic->ExcecuteOneStep();
        if (fSimulationData->IsOnePhaseQ()) {
            return;
        }
        this->UpdateFluxes_at_n();
    }
    
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    this->UpdateMemory_at_n();
    
    fHyperbolic->ExcecuteOneStep();
    
    this->UpdateMemory_at_n();

    
}

void TRMSegregatedAnalysis::ExcecuteOneStep(bool IsActiveQ){

   this->SegregatedIteration(IsActiveQ);

}

/** @brief Update memory using the Transfer object at state n */
void TRMSegregatedAnalysis::UpdateMemory_at_n(){
    
    fSimulationData->SetCurrentStateQ(true);
    Parabolic()->UpdateMemory();
    
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    Hyperbolic()->UpdateMemory();
    fTransfer->Reciprocal_Memory_Transfer(fParabolic->Mesh(), fHyperbolic->Mesh());

}

/** @brief Update memory using the Transfer object */
void TRMSegregatedAnalysis::UpdateMemory(){
    
    fSimulationData->SetCurrentStateQ(false);
    Parabolic()->UpdateMemory();
    
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    Hyperbolic()->UpdateMemory();
    fTransfer->Reciprocal_Memory_Transfer(fParabolic->Mesh(), fHyperbolic->Mesh());
    
}

/** @brief Update memory using the Transfer object at state n */
void TRMSegregatedAnalysis::UpdateFluxes_at_n(){

    fParabolic->UpdateMemory_at_n();
    fTransfer->un_To_Transport_Mesh(fParabolic->Mesh(), fHyperbolic->Mesh(),true);
    fTransfer->un_To_Transport_Mesh(fParabolic->Mesh(), fHyperbolic->Mesh(),false);
    
    
}

void TRMSegregatedAnalysis::PostProcessStep(bool IsActiveQ){
    
    if(IsActiveQ){
        fParabolic->PostProcessStep();
    }
    
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    fHyperbolic->PostProcessStep();
    
}