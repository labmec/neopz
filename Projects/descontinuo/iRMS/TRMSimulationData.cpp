//
//  TRMSimulationData.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMSimulationData.h"


/** @brief Initialize the raw data */
TRMSimulationData::TRMSimulationData(){
    
    /** @brief initial state */
    fIsInitialStateQ = false;
    
    /** @brief current time state */
    fIsCurrentStateQ =  false;
    
    /** @brief one-phase flow */
    fIsOnePhaseQ   = false;
    
    /** @brief two-phase flow */
    fIsTwoPhaseQ   = false;
    
    /** @brief three-phase flow */
    fIsThreePhaseQ = false;
    
    /** @brief Definition of the flow system one - two or three phase */
    fSystemType.Resize(0);
    
    /** @brief Definition gravity field */
    fg.Resize(0);
    
    /** @brief Material identifier for interfaces */
    fInterface_mat_Id = 99999;
    
    /** @brief Store time values to be reported */
    fReportingTimes.Resize(0);
    
    /** @brief ntime steps */
    fn_steps = 0;
    
    /** @brief Initial time */
    ftime_0 = 0.0;
    
    /** @brief Final time */
    ftime_n = 0.0;
    
    /** @brief Time step */
    fdt = 0.0;
    
    /** @brief Time step */
    ftime = 0.0;
    
    /** @brief Min time step */
    fdt_min = 0.0;
    
    /** @brief Max time step */
    fdt_max = 0.0;
    
    /** @brief Increment dt factor */
    fdt_up = 0.0;
    
    /** @brief Decrement dt factor */
    fdt_down = 0.0;
    
    /** @brief Use of quasi newton method */
    fIsQuasiNewtonQ = false;
    
    /** @brief Autopointer of all the petrophysics data */
    fPetroPhysics = NULL;
    
    /** @brief phase alpha */
    fPhase_alpha = NULL;
    
    /** @brief phase beta */
    fPhase_beta = NULL;
    
    /** @brief phase gamma */
    fPhase_gamma = NULL;

    /** @brief Stores the spatial information given in maps */
    fMap = NULL;
    
    /** @brief L2 projection material id for gradient reconstruction */
    fl2_projection_material_id = 2001;
    
    /** @brief Skeleton dfault material id for MHM substructuring */
    fSkeleton_material_id = 5001;
    
    /** @brief Define the use of linear gradient reconstruction */
    fUseGradientRQ = false;
    
}

/** @brief destructor */
TRMSimulationData::~TRMSimulationData(){
    
}

/** @brief Set autopointer of the RawData */
void TRMSimulationData::SetRawData(TPZAutoPointer<TRMRawData> &RawData){
    fRawData = RawData;
    SetMap(RawData->fMap);
    SetGravity(RawData->fg);
    SetSystemType(RawData->fSystemType,RawData->fPhases);
    SetTimeControls(RawData->fn_steps, RawData->fdt, RawData->fdt_up, RawData->fdt_down, RawData->fdt_max, RawData->fdt_min, RawData->fReportingTimes);
    SetNumericControls(RawData->fn_corrections, RawData->fepsilon_res, RawData->fepsilon_cor, RawData->fIsQuasiNewtonQ);
}

/** @brief Setup reporting times and time step size */
void TRMSimulationData::SetTimeControls(int n_times, STATE dt, STATE dt_in, STATE dt_de, STATE dt_max, STATE dt_min, TPZStack< std::pair< STATE , bool> , 500 > ReportingTimes){
    fdt = dt;
    fn_steps = n_times;

    for (int it = 0; it < ReportingTimes.size(); it++) {
        fReportingTimes.Push(ReportingTimes[it].first);
        fReportingTimesMixedQ.Push(ReportingTimes[it].second);
    }
    fdt_max     = dt_max;
    fdt_min     = dt_min;
    fdt_up      = dt_in;
    fdt_down    = dt_de;
    ftime_0     = fReportingTimes[0];
    ftime_n     = fReportingTimes[ReportingTimes.size()-1];
}

/** @brief Setup reporting times and time step size */
void TRMSimulationData::SetNumericControls(int n_corrections, STATE epsilon_res, STATE epsilon_cor, bool IsQuasiNewtonQ){
    fn_corrections  = n_corrections;
    fepsilon_res    = epsilon_res;
    fepsilon_cor    = epsilon_cor;
    fIsQuasiNewtonQ = IsQuasiNewtonQ;
}

/** @brief Set phase alpha */
void TRMSimulationData::SetPhaseAlpha(TPZAutoPointer<TRMPhaseProperties> &alpha)
{
    fPhase_alpha = alpha;
}

/** @brief Get phase alpha */
TPZAutoPointer<TRMPhaseProperties> & TRMSimulationData::AlphaProp()
{
    return fPhase_alpha;
}

/** @brief Set phase beta */
void TRMSimulationData::SetPhaseBeta(TPZAutoPointer<TRMPhaseProperties> &beta)
{
    fPhase_beta = beta;
}

/** @brief Get phase beta */
TPZAutoPointer<TRMPhaseProperties> & TRMSimulationData::BetaProp()
{
    return fPhase_beta;
}

/** @brief Set phase gamma */
void TRMSimulationData::SetPhaseGamma(TPZAutoPointer<TRMPhaseProperties> &gamma)
{
    fPhase_gamma = gamma;
}

/** @brief Get phase gamma */
TPZAutoPointer<TRMPhaseProperties> & TRMSimulationData::GammaProp()
{
    return fPhase_gamma;
}



