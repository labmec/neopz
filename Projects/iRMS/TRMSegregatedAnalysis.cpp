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

    //    this->UpdateMemory();
    this->UpdateMemory_at_n(); // @omar:: It is time to verify
    fParabolic->ExcecuteOneStep();

    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    this->UpdateFluxes_at_n();
    
//    this->UpdateMemory_at_n();

    fHyperbolic->ExcecuteOneStep();
    
//    this->UpdateMemory_at_n();    

    
}

void TRMSegregatedAnalysis::ExcecuteOneStep(){

    
    REAL dt_min    = fSimulationData->dt_min();
    REAL dt_max    = fSimulationData->dt_max();
    REAL dt_up     = fSimulationData->dt_up();
    REAL dt_down   = fSimulationData->dt_down();
    REAL dt        = fSimulationData->dt();
    
    REAL epsilon_res = this->SimulationData()->epsilon_res();
    REAL epsilon_cor = this->SimulationData()->epsilon_cor();
    int n  =   this->SimulationData()->n_corrections();
    
    ferror_flux_pressure = 1.0;
    ferror_saturation = 1.0;
    fdx_norm_flux_pressure = 1.0;
    fdx_norm_saturation = 1.0;
    
    bool IsConverged_eQ = false;
    bool IsConverged_dQ = false;
    bool IsConverged_iQ = false;
    
    for (int k = 1; k <= n; k++) {

        this->SegregatedIteration();
        
        ferror_flux_pressure = fParabolic->error_norm();
        ferror_saturation = fHyperbolic->error_norm();
        
        fdx_norm_flux_pressure = fParabolic->dx_norm();
        fdx_norm_saturation = fHyperbolic->dx_norm();
        
        IsConverged_eQ = (ferror_flux_pressure < epsilon_res) &&  (ferror_saturation < epsilon_res);
        IsConverged_dQ = (fdx_norm_flux_pressure < epsilon_cor) &&  (fdx_norm_saturation < epsilon_cor);
        IsConverged_iQ = (fParabolic->k_ietrarions() <= 5) &&  (fHyperbolic->k_ietrarions() <= 10);
        
        if((IsConverged_eQ || IsConverged_dQ) &&  IsConverged_iQ)
        {
            std::cout << "Segregated:: Converged with iterations:  " << k << "; error: " << ferror_flux_pressure + ferror_saturation <<  "; dx: " << fdx_norm_flux_pressure + fdx_norm_saturation << std::endl;
            
            // update Time value
            REAL current_time = fSimulationData->t() + fSimulationData->dt();
            fSimulationData->SetTime(current_time);
            
            if (k <= 1 && dt_max > dt && dt_up > 1.0) {
                dt *= dt_up;
                if(dt_max < dt ){
                    fSimulationData->Setdt(dt_max);
                }
                else{
                    fSimulationData->Setdt(dt);
                }
                std::cout << "Segregated:: Increasing time step to " << fSimulationData->dt()/86400.0 << "; (day): " << std::endl;
            }
            
            this->UpdateGlobalSolution();
            return;
        }
        
        if(k == n  && dt > dt_min && dt_down < 1.0){
            dt *= dt_down;
            if(dt_min > dt ){
                fSimulationData->Setdt(dt_min);
            }
            else{
                fSimulationData->Setdt(dt);
            }
            std::cout << "Segregated:: Decreasing time step to " << fSimulationData->dt()/86400.0 << "; (day): " << std::endl;
            std::cout << "Segregated:: Restarting current time step correction " << std::endl;
            
            this->KeepGlobalSolution();
            k = 1;
        }
        
        
    }
    
    // update Time value
    REAL current_time = fSimulationData->t() + fSimulationData->dt();
    fSimulationData->SetTime(current_time);
    
    std::cout << "Segregated:: Exit max iterations with min dt:  " << fSimulationData->dt()/86400.0 << "; (day) " << "; error: " << ferror_flux_pressure + ferror_saturation <<  "; dx: " << fdx_norm_flux_pressure + fdx_norm_saturation << std::endl;

}

/** @brief Update memory using the Transfer object at REAL n */
void TRMSegregatedAnalysis::UpdateMemory_at_n(){
    
    fSimulationData->SetCurrentStateQ(true);
    fParabolic->UpdateMemory_at_n();
    
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    Hyperbolic()->UpdateMemory_at_n();
    fTransfer->Reciprocal_Memory_Transfer(fParabolic->Mesh(), fHyperbolic->Mesh());

}

/** @brief Update memory using the Transfer object */
void TRMSegregatedAnalysis::UpdateMemory(){
    
    fSimulationData->SetCurrentStateQ(false);
    fParabolic->UpdateMemory();
    
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    Hyperbolic()->UpdateMemory();
    fTransfer->Reciprocal_Memory_Transfer(fParabolic->Mesh(), fHyperbolic->Mesh());
    
}

/** @brief Update memory using the Transfer object at REAL n */
void TRMSegregatedAnalysis::UpdateFluxes_at_n(){

    fParabolic->UpdateMemory_at_n();
    fTransfer->un_To_Transport_Mesh(fParabolic->Mesh(), fHyperbolic->Mesh(),true);
    fTransfer->un_To_Transport_Mesh(fParabolic->Mesh(), fHyperbolic->Mesh(),false);
    
    
}

/** @brief update global state for the new euler step */
void TRMSegregatedAnalysis::UpdateGlobalSolution(){
    
    fParabolic->X() = fParabolic->X_n();
    fHyperbolic->X() = fHyperbolic->X_n();
    
}

/** @brief keep global last state for restart a euler step */
void TRMSegregatedAnalysis::KeepGlobalSolution(){
    
    fParabolic->X_n() = fParabolic->X();
    fHyperbolic->X_n() = fHyperbolic->X();
    
}

void TRMSegregatedAnalysis::PostProcessStep(bool draw_mixed_mapQ){
    
    if (draw_mixed_mapQ) {
        fParabolic->PostProcessStep();
    }
    
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    fHyperbolic->PostProcessStep();
    
}
