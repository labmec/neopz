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
    
//    STATE epsilon_res = this->SimulationData()->epsilon_res();
//    STATE epsilon_cor = this->SimulationData()->epsilon_cor();
//    int n  =   this->SimulationData()->n_corrections();
    
//    for (int k = 1; k <= n; k++) {

//        this->SegregatedIteration(IsActiveQ);

//        if(ferror < epsilon_res || fdx_norm < epsilon_cor)
//        {
//            std::cout << "Hyperbolic:: Converged with iterations:  " << k << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
//            if (k == 1 && dt_max > dt && dt_up > 1.0) {
//                dt *= dt_up;
//                if(dt_max < dt ){
//                    fSimulationData->Setdt(dt_max);
//                }
//                else{
//                    fSimulationData->Setdt(dt);
//                }
//                std::cout << "Hyperbolic:: Increasing time step to " << fSimulationData->dt()/86400.0 << "; (day): " << std::endl;
//            }
//            
//            fX = fX_n;
//            return;
//        }
//        
//        if(k == n  && dt > dt_min && dt_down < 1.0){
//            dt *= dt_down;
//            if(dt_min > dt ){
//                fSimulationData->Setdt(dt_min);
//            }
//            else{
//                fSimulationData->Setdt(dt);
//            }
//            std::cout << "Hyperbolic:: Decreasing time step to " << fSimulationData->dt()/86400.0 << "; (day): " << std::endl;
//            std::cout << "Hyperbolic:: Restarting current time step correction " << std::endl;
//            
//            this->SimulationData()->SetCurrentStateQ(false);
//            this->LoadSolution(fX);
//            
//            this->UpdateMemory();
//            this->AssembleResidual();
//            fR = this->Rhs();
//            
//            this->SimulationData()->SetCurrentStateQ(true);
//            this->UpdateMemory_at_n();
//            fX_n = fX;
//            k = 1;
//        }
//        
//        
//    }
//    
//    std::cout << "Hyperbolic:: Exit max iterations with min dt:  " << fSimulationData->dt()/86400.0 << "; (day) " << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
    

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