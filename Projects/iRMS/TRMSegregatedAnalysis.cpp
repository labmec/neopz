//
//  TRMSegregatedAnalysis.cpp
//  PZ
//
//  Created by Omar on 7/18/16.
//
//

#include "TRMSegregatedAnalysis.h"
#define NS

TRMSegregatedAnalysis::TRMSegregatedAnalysis() : TPZAnalysis() {
    
    /** @brief define the simulation data */
    fSimulationData = NULL;
    
    /** @brief define the transfer matrices */
    fTransfer = NULL;
    
    /** @brief define the elliptic system */
    fElliptic = NULL;
    
    /** @brief define the parabolic system */
    fParabolic = NULL;
    
    /** @brief define the hyperbolic system */
    fHyperbolic = NULL;
    
    /** @brief Residue error for displacement */
    ferror_displacement = 0.0;
    
    /** @brief Residue error for flux - pressure */
    ferror_flux_pressure = 0.0;
    
    /** @brief Residue error for saturations */
    ferror_saturation = 0.0;
    
    /** @brief Correction variation for displacement */
    fdx_norm_displacement = 0.0;
    
    /** @brief Correction variation for flux - pressure */
    fdx_norm_flux_pressure = 0.0;
    
    /** @brief Correction variation for saturations */
    fdx_norm_saturation = 0.0;
    
    /** @brief number of segregated corrections */
    fk_iterations = 0;
    
}

TRMSegregatedAnalysis::~TRMSegregatedAnalysis(){
    
}

/** @brief Copy constructor $ */
TRMSegregatedAnalysis::TRMSegregatedAnalysis(const TRMSegregatedAnalysis &copy)
{
    fSimulationData         = copy.fSimulationData;
    fTransfer               = copy.fTransfer;
    fElliptic               = copy.fElliptic;
    fParabolic              = copy.fParabolic;
    fHyperbolic             = copy.fHyperbolic;
    ferror_displacement     = copy.ferror_displacement;
    ferror_flux_pressure    = copy.ferror_flux_pressure;
    ferror_saturation       = copy.ferror_saturation;
    fdx_norm_displacement   = copy.fdx_norm_displacement;
    fdx_norm_flux_pressure  = copy.fdx_norm_flux_pressure;
    fdx_norm_saturation     = copy.fdx_norm_saturation;
    fk_iterations           = copy.fk_iterations;
}

/** @brief Copy assignemnt operator $ */
TRMSegregatedAnalysis & TRMSegregatedAnalysis::operator=(const TRMSegregatedAnalysis &other)
{
    if (this != & other) {  // prevent self-assignment
        
        fSimulationData         = other.fSimulationData;
        fTransfer               = other.fTransfer;
        fElliptic               = other.fElliptic;
        fParabolic              = other.fParabolic;
        fHyperbolic             = other.fHyperbolic;
        ferror_displacement     = other.ferror_displacement;
        ferror_flux_pressure    = other.ferror_flux_pressure;
        ferror_saturation       = other.ferror_saturation;
        fdx_norm_displacement   = other.fdx_norm_displacement;
        fdx_norm_flux_pressure  = other.fdx_norm_flux_pressure;
        fdx_norm_saturation     = other.fdx_norm_saturation;
        fk_iterations           = other.fk_iterations;
    }
    return *this;
}

/** @brief Resize and fill residue and solution vectors */
void TRMSegregatedAnalysis::AdjustVectors(){
    
    fElliptic->AdjustVectors();
    fParabolic->AdjustVectors();
    fHyperbolic->AdjustVectors();
}

void TRMSegregatedAnalysis::SegregatedIteration(){

    this->UpdateMemory_at_n();
    
    fParabolic->ExcecuteOneStep();

    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    this->UpdateFluxes_at_n();
    this->UpdateMemory_at_n();

    fHyperbolic->ExcecuteOneStep();
    this->UpdateMemory_at_n();    
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
    bool IsConverged_iQ = true;
    bool MustRestartQ = false;
    
    this->UpdateMemory(); // last average values
    
    for (int k = 1; k <= n; k++) {

        fk_iterations = k;
        this->SegregatedIteration();
        
        ferror_flux_pressure = fParabolic->error_norm();
        ferror_saturation = fHyperbolic->error_norm();
        
        fdx_norm_flux_pressure = fParabolic->dx_norm();
        fdx_norm_saturation = fHyperbolic->dx_norm();
        
        IsConverged_eQ = (ferror_flux_pressure < epsilon_res) &&  (ferror_saturation < epsilon_res);
        IsConverged_dQ = (fdx_norm_flux_pressure < epsilon_cor) &&  (fdx_norm_saturation < epsilon_cor);
        
        if (!fSimulationData->IsOnePhaseQ()) {
            IsConverged_iQ = (fParabolic->k_ietrarions() <= 10) &&  (fHyperbolic->k_ietrarions() <= 10);
        }

//        MustRestartQ = MustRestartStep();
        MustRestartQ = (fHyperbolic->k_ietrarions() == n) ;
        
        if((k == n || MustRestartQ)  && dt > dt_min && dt_down < 1.0){
            dt *= dt_down;
            if(dt_min > dt ){
                fSimulationData->Setdt(dt_min);
            }
            else{
                fSimulationData->Setdt(dt);
            }
            std::cout << "Segregated:: Decreasing time step to " << fSimulationData->dt()/86400.0 << "; (day): " << std::endl;
            if (MustRestartQ) {
                std::cout << "Segregated:: Force restarting current time step correction " << std::endl;
            }
            std::cout << "Segregated:: Restarting current time step correction " << std::endl;
            
            this->KeepGlobalSolution();
            k = 0;
            continue;
        }
        
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
    
    }
    
    // update Time value
    REAL current_time = fSimulationData->t() + fSimulationData->dt();
    fSimulationData->SetTime(current_time);
    
    std::cout << "Segregated:: Exit max iterations with min dt:  " << fSimulationData->dt()/86400.0 << "; (day) " << "; error: " << ferror_flux_pressure + ferror_saturation <<  "; dx: " << fdx_norm_flux_pressure + fdx_norm_saturation << std::endl;

}

/** @brief Execute a segregated iteration with fixed stress  */
void TRMSegregatedAnalysis::SegregatedIteration_Fixed_Stress(){
    

    this->UpdateMemory_at_n();
    
    // Undarined response
    if (fSimulationData->IsInitialStateQ()) {

        if (fSimulationData->IsGeomechanicQ()) {

//            fElliptic->ExcecuteOneStep();
//            fTransfer->elliptic_To_parabolic(fElliptic->Mesh(), fParabolic->Mesh());
//            fParabolic->ExcecuteOneStep();
//            fTransfer->parabolic_To_elliptic(fParabolic->Mesh(), fElliptic->Mesh());
//
//            if (!fSimulationData->IsOnePhaseQ()) {
//                fTransfer->parabolic_To_hyperbolic_volumetric(fParabolic->Mesh(), fHyperbolic->Mesh());
//                fTransfer->elliptic_To_hyperbolic(fElliptic->Mesh(), fHyperbolic->Mesh());
//            }
            
            // Fixed stress Iteration 1
            if (fSimulationData->IsOnePhaseQ()) {
                fParabolic->ExcecuteOneStep();
            }
            else{
                this->Segregated_p_h_Iteration();
            }
            if (fSimulationData->IsGeomechanicQ()) {
                fTransfer->parabolic_To_elliptic(fParabolic->Mesh(), fElliptic->Mesh());
                fElliptic->ExcecuteOneStep();
                fTransfer->elliptic_To_parabolic(fElliptic->Mesh(), fParabolic->Mesh());
                if (fSimulationData->IsTwoPhaseQ()) {
                    fTransfer->elliptic_To_hyperbolic(fElliptic->Mesh(), fHyperbolic->Mesh());
                }
            }

            
        }else{
            fParabolic->ExcecuteOneStep();
        }
        
        return;
        
    }
    
    // Fixed stress Iteration 1
    if (fSimulationData->IsOnePhaseQ()) {
        fParabolic->ExcecuteOneStep();
    }
    else{
        this->Segregated_p_h_Iteration();
    }
    if (fSimulationData->IsGeomechanicQ()) {
        fTransfer->parabolic_To_elliptic(fParabolic->Mesh(), fElliptic->Mesh());
        fElliptic->ExcecuteOneStep();
        fTransfer->elliptic_To_parabolic(fElliptic->Mesh(), fParabolic->Mesh());
        if (fSimulationData->IsTwoPhaseQ()) {
            fTransfer->elliptic_To_hyperbolic(fElliptic->Mesh(), fHyperbolic->Mesh());
        }
    }
    
    if (!fSimulationData->IsGeomechanicQ()) {
        return;
    }

    
}

/** @brief Execute a segregated iteration between parabolic and hyperbolic operators  */
void TRMSegregatedAnalysis::Segregated_p_h_Iteration(){

    fParabolic->ExcecuteOneStep();
    this->UpdateFluxes_at_n();
    this->UpdateMemory_at_n();
    fHyperbolic->ExcecuteOneStep();
    this->UpdateMemory_at_n();
    
//    fParabolic->ExcecuteOneStep();
//    this->UpdateFluxes_at_n();
//    this->UpdateMemory_at_n();
//    fHyperbolic->ExcecuteOneStep();
//    this->UpdateMemory_at_n();
    
}

void TRMSegregatedAnalysis::ExcecuteOneStep_Fixed_Stress(){
    
    
    REAL dt_min    = fSimulationData->dt_min();
    REAL dt_max    = fSimulationData->dt_max();
    REAL dt_up     = fSimulationData->dt_up();
    REAL dt_down   = fSimulationData->dt_down();
    REAL dt        = fSimulationData->dt();
    
    REAL epsilon_res = this->SimulationData()->epsilon_res();
    REAL epsilon_cor = this->SimulationData()->epsilon_cor();
    int n  =   this->SimulationData()->n_corrections();
    
    ferror_displacement = 1.0;
    ferror_flux_pressure = 1.0;
    ferror_saturation = 1.0;
    fdx_norm_displacement = 1.0;
    fdx_norm_flux_pressure = 1.0;
    fdx_norm_saturation = 1.0;
    
    bool IsConverged_eQ = false;
    bool IsConverged_dQ = false;
    bool IsConverged_iQ = true;
    bool MustRestartQ = false;
    
    this->UpdateMemory(); // last average values
    
    for (int k = 1; k <= n; k++) {
        
        this->SegregatedIteration_Fixed_Stress();
        fk_iterations = k;
        
        ferror_displacement  = fElliptic->error_norm();
        ferror_flux_pressure = fParabolic->error_norm();
        ferror_saturation    = fHyperbolic->error_norm();
        
        fdx_norm_displacement   = fElliptic->dx_norm();
        fdx_norm_flux_pressure  = fParabolic->dx_norm();
        fdx_norm_saturation     = fHyperbolic->dx_norm();
                
        IsConverged_eQ = (ferror_flux_pressure < epsilon_res) &&  (ferror_saturation < epsilon_res);
        IsConverged_dQ = (fdx_norm_flux_pressure < epsilon_cor) &&  (fdx_norm_saturation < epsilon_cor);
        
        if (!fSimulationData->IsOnePhaseQ()) {
            IsConverged_iQ = (fParabolic->k_ietrarions() <= 10) &&  (fHyperbolic->k_ietrarions() <= 10);
        }
        
        MustRestartQ = MustRestartStep();
        MustRestartQ = (fHyperbolic->k_ietrarions() == n) ;
        
        if((k == n || MustRestartQ)  && dt > dt_min && dt_down < 1.0){
            dt *= dt_down;
            if(dt_min > dt ){
                fSimulationData->Setdt(dt_min);
            }
            else{
                fSimulationData->Setdt(dt);
            }
            std::cout << "Segregated:: Decreasing time step to " << fSimulationData->dt()/86400.0 << "; (day): " << std::endl;
            if (MustRestartQ) {
                std::cout << "Segregated:: Force restarting current time step correction " << std::endl;
            }
            std::cout << "Segregated:: Restarting current time step correction " << std::endl;
            
            this->KeepGlobalSolution();
            k = 0;
            continue;
        }
        
        if((IsConverged_eQ && IsConverged_dQ) &&  IsConverged_iQ)
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
        
    }
    
    // update Time value
    REAL current_time = fSimulationData->t() + fSimulationData->dt();
    fSimulationData->SetTime(current_time);
    
    std::cout << "Segregated:: Exit max iterations with min dt:  " << fSimulationData->dt()/86400.0 << "; (day) " << "; error: " << ferror_flux_pressure + ferror_saturation <<  "; dx: " << fdx_norm_flux_pressure + fdx_norm_saturation << std::endl;
    
}

/** @brief Update memory using the Transfer object at REAL n */
void TRMSegregatedAnalysis::UpdateMemory_at_n(){
    
    fSimulationData->SetCurrentStateQ(true);
    
#ifdef NS
    
    if (fSimulationData->IsGeomechanicQ()) {
        fElliptic->UpdateMemory_at_n();
        
        // values Ae - Ap operators
        fTransfer->elliptic_To_parabolic(fElliptic->Mesh(), fParabolic->Mesh());
        fTransfer->parabolic_To_elliptic(fParabolic->Mesh(), fElliptic->Mesh());
    }
    
    fParabolic->UpdateMemory_at_n();
    
    
    if (!fSimulationData->IsOnePhaseQ()) {
        fHyperbolic->UpdateMemory_at_n();
        
        // average values Ah - Ap operators
        fTransfer->parabolic_To_hyperbolic_volumetric(fParabolic->Mesh(), fHyperbolic->Mesh());
        fTransfer->hyperbolic_To_parabolic_volumetric(fHyperbolic->Mesh(), fParabolic->Mesh());
        
        if (fSimulationData->IsGeomechanicQ()) {
            fTransfer->hyperbolic_To_elliptic(fHyperbolic->Mesh(), fElliptic->Mesh());
            fTransfer->elliptic_To_hyperbolic(fElliptic->Mesh(), fHyperbolic->Mesh());
        }
        
    }
    

    
#else
    
    fParabolic->UpdateMemory_at_n();
    
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    Hyperbolic()->UpdateMemory_at_n();
    
    if (fSimulationData->TransporResolution().first) {
        fTransfer->Reciprocal_Memory_TransferII(fParabolic->Mesh(), fHyperbolic->Mesh());
    }
    else{
        fTransfer->Reciprocal_Memory_Transfer(fParabolic->Mesh(), fHyperbolic->Mesh());
    }
    
#endif

}

/** @brief Update memory using the Transfer object */
void TRMSegregatedAnalysis::UpdateMemory(){
    
    fSimulationData->SetCurrentStateQ(false);
    
#ifdef NS
    
    if (fSimulationData->IsGeomechanicQ()) {
        fElliptic->UpdateMemory();
        
        // values Ae - Ap operators
        fTransfer->elliptic_To_parabolic(fElliptic->Mesh(), fParabolic->Mesh());
        fTransfer->parabolic_To_elliptic(fParabolic->Mesh(), fElliptic->Mesh());
    }
    
    fParabolic->UpdateMemory();
    
    if (!fSimulationData->IsOnePhaseQ()) {
        fHyperbolic->UpdateMemory();
        
        // average values Ah - Ap operators
        fTransfer->parabolic_To_hyperbolic_volumetric(fParabolic->Mesh(), fHyperbolic->Mesh());
        fTransfer->hyperbolic_To_parabolic_volumetric(fHyperbolic->Mesh(), fParabolic->Mesh());
        
        // average values Ah - Ae operators
        if (fSimulationData->IsGeomechanicQ()) {
            fTransfer->hyperbolic_To_elliptic(fHyperbolic->Mesh(), fElliptic->Mesh());
            fTransfer->elliptic_To_hyperbolic(fElliptic->Mesh(), fHyperbolic->Mesh());
        }
        
    }
    
#else
    
    fParabolic->UpdateMemory();
    
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    Hyperbolic()->UpdateMemory();
    
    if (fSimulationData->TransporResolution().first) {
        fTransfer->Reciprocal_Memory_TransferII(fParabolic->Mesh(), fHyperbolic->Mesh());
    }
    else{
        fTransfer->Reciprocal_Memory_Transfer(fParabolic->Mesh(), fHyperbolic->Mesh());
    }
    
#endif
    
}


/** @brief Update memory using the Transfer object at REAL n */
void TRMSegregatedAnalysis::UpdateFluxes_at_n(){

    fParabolic->UpdateMemory_at_n();
    
    fTransfer->parabolic_To_hyperbolic_interfaces(fParabolic->Mesh(), fHyperbolic->Mesh(),true);
    fTransfer->parabolic_To_hyperbolic_interfaces(fParabolic->Mesh(), fHyperbolic->Mesh(),false);

}

/** @brief update global state for the new euler step */
void TRMSegregatedAnalysis::UpdateGlobalSolution(){
    
#ifdef NS

    fElliptic->X()  = fElliptic->X_n();
    fParabolic->X() = fParabolic->X_n();
    fHyperbolic->X() = fHyperbolic->X_n();
    
#else
    
    fParabolic->X() = fParabolic->X_n();
    fHyperbolic->X() = fHyperbolic->X_n();
    
#endif
    
}

/** @brief keep global last state for restart a euler step */
void TRMSegregatedAnalysis::KeepGlobalSolution(){
    
#ifdef NS
    
    fElliptic->X()  = fElliptic->X();
    fParabolic->X() = fParabolic->X();
    fHyperbolic->X() = fHyperbolic->X();    
    
#else
    
    fParabolic->X() = fParabolic->X();
    fHyperbolic->X() = fHyperbolic->X();
    
#endif
    
}

/** @brief keep global last state for restart a euler step */
bool TRMSegregatedAnalysis::MustRestartStep(){
    
    int n_data = fHyperbolic->X_n().Rows();
    REAL epsilon = 1.0e-3;

    for (long i = 0; i < n_data; i++) {
        if ( (1.0 - fHyperbolic->X_n()(i,0)) < + epsilon ) {
//            fHyperbolic->X_n().Print("sw = ");
            return true;
        }
        
        if ( (fHyperbolic->X_n()(i,0)) < - epsilon ) {
//            fHyperbolic->X_n().Print("sw = ");
            return true;
        }
    }
    
    return false;
    
}


void TRMSegregatedAnalysis::PostProcessStep(bool draw_mixed_mapQ){
    
    if (draw_mixed_mapQ) {
        
        if (fSimulationData->IsGeomechanicQ()) {
            fElliptic->PostProcessStep();
        }
    
        fParabolic->PostProcessStep();
    }
    
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    fHyperbolic->PostProcessStep();
    
}
