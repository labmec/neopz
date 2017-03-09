//
//  TPZSegregatedSolver.cpp
//  PZ
//
//  Created by omar on 3/8/17.
//
//

#include "TPZSegregatedSolver.h"

TPZSegregatedSolver::TPZSegregatedSolver() : TPZAnalysis() {
    
    /** @brief define the simulation data */
    fSimulationData = NULL;
    
    /** @brief define the elliptic component */
    felliptic = new TPZElasticAnalysis;
    
    /** @brief define the parabolic component */
    fparabolic = new TPZFLuxPressureAnalysis;
    
    /** @brief define the tranfer object data */
    ftransfer = NULL;
    
    /** @brief Residue error */
    ferror    = 1.0;
    
    /** @brief Correction variation */
    fdx_norm  = 1.0;
    
    /** @brief number of newton corrections */
    fk_iterations = 0;
    
}

TPZSegregatedSolver::~TPZSegregatedSolver(){
    
}

/** @brief Copy constructor $ */
TPZSegregatedSolver::TPZSegregatedSolver(const TPZSegregatedSolver &copy)
{
    fSimulationData = copy.fSimulationData;
    felliptic       = copy.felliptic;
    fparabolic       = copy.fparabolic;
    ferror          = copy.ferror;
    fdx_norm        = copy.fdx_norm;
    
}

/** @brief Copy assignemnt operator $ */
TPZSegregatedSolver & TPZSegregatedSolver::operator=(const TPZSegregatedSolver &other)
{
    if (this != & other) {  // prevent self-assignment
        
        fSimulationData = other.fSimulationData;
        felliptic       = other.felliptic;
        fparabolic      = other.fparabolic;
        ferror          = other.ferror;
        fdx_norm        = other.fdx_norm;
    }
    return *this;
}


/** @brief execute the evolutionary problem */
void TPZSegregatedSolver::Run_Evolution(std::string elliptic, std::string parabolic){
    
    int interval = 20;
    int n = fSimulationData->n_steps();
    REAL time = 0.0;
    REAL dt = this->SimulationData()->dt();
    
    this->PostProcessStep(elliptic,parabolic);    // Initial condition
    
    for (int i = 1; i <= n; i++) {
        
        time = i * dt;
        this->SimulationData()->SetTime(time);
        this->ExcecuteOneStep();
        
        if (i%interval == 0) {
            std::cout<< "Segregated solver:: Reporting at time (s) = " << time << std::endl;
            this->PostProcessStep(elliptic,parabolic);
        }
        
        std::cout<< "Segregated solver:: Current time (s) = " << time << std::endl;
    }
    
}

void TPZSegregatedSolver::ExcecuteOneStep(){
    
    STATE epsilon_res = this->SimulationData()->epsilon_res();
    STATE epsilon_cor = this->SimulationData()->epsilon_cor();
    int n  =   this->SimulationData()->n_corrections();
    
    bool IsConverged_eQ = false;
    bool IsConverged_dQ = false;
    
    UpdateState();// Required
    Update_at_n_State();// Required
    
    for (int k = 1; k <= n; k++) {
        
        this->SegregatedIteration();
        
        REAL error_e = felliptic->error_norm();
        REAL error_p = fparabolic->error_norm();
        
        REAL dx_norm_e = felliptic->dx_norm();
        REAL dx_norm_p = fparabolic->dx_norm();
        
        IsConverged_eQ = (error_e < epsilon_res) &&  (error_p < epsilon_res);
        IsConverged_dQ = (dx_norm_e < epsilon_cor) &&  (dx_norm_p < epsilon_cor);
        
        if(/* IsConverged_eQ || */ IsConverged_dQ)
        {
            std::cout << "Geomechanic Coupling:: Converged with iterations:  " << k << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
            this->UpdateGlobalSolution();
            this->Update_at_n_State();
            return;
        }
        
    }
    
    std::cout << "Geomechanic Coupling:: Exit max iterations with min dt:  " << fSimulationData->dt() << "; (secs) " << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
    
}

/** @brief Execute a segregated iteration  */
void TPZSegregatedSolver::SegregatedIteration(){
    
    // Fixed-Stress split
    fparabolic->ExcecuteOneStep();
    ftransfer->parabolic_To_elliptic(fparabolic->Mesh(), felliptic->Mesh());
    
    felliptic->ExcecuteOneStep();
    ftransfer->elliptic_To_parabolic(felliptic->Mesh(), fparabolic->Mesh());
    
    Update_at_n_State();
    
}

/** @brief update last state solution */
void TPZSegregatedSolver::UpdateState(){
    
    fSimulationData->SetCurrentStateQ(false);
    felliptic->UpdateState();
    fparabolic->UpdateState();
    
    ftransfer->elliptic_To_parabolic(felliptic->Mesh(), fparabolic->Mesh());
    ftransfer->parabolic_To_elliptic(fparabolic->Mesh(), felliptic->Mesh());
    
}

/** @brief update current state solution */
void TPZSegregatedSolver::Update_at_n_State(){
    
    fSimulationData->SetCurrentStateQ(true);
    felliptic->Update_at_n_State();
    fparabolic->Update_at_n_State();
    
    ftransfer->elliptic_To_parabolic(felliptic->Mesh(), fparabolic->Mesh());
    ftransfer->parabolic_To_elliptic(fparabolic->Mesh(), felliptic->Mesh());
}

void TPZSegregatedSolver::PostProcessStep(std::string elliptic, std::string parabolic){
    
    felliptic->PostProcessStep(elliptic);
    fparabolic->PostProcessStep(parabolic);
}


/** @brief update global state for the new euler step */
void TPZSegregatedSolver::UpdateGlobalSolution(){
    
    felliptic->X() = felliptic->X_n();
    fparabolic->X() = fparabolic->X_n();
    
}

/** @brief keep global last state for restart a euler step */
void TPZSegregatedSolver::KeepGlobalSolution(){
    
    felliptic->X_n() = felliptic->X();
    fparabolic->X_n() = fparabolic->X();
    
}