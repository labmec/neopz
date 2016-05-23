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
    
    /** @brief Min time step */
    fdt_min = 0.0;
    
    /** @brief Max time step */
    fdt_max = 0.0;
    
    /** @brief Increment dt factor */
    fdt_up = 0.0;
    
    /** @brief Decrement dt factor */
    fdt_down = 0.0;
    
    /** @brief phase alpha */
    fPhase_alpha = NULL;
    
    /** @brief phase beta */
    fPhase_Beta = NULL;
    
    /** @brief phase gamma */
    fPhase_gamma = NULL;
    
}

/** @brief destructor */
TRMSimulationData::~TRMSimulationData(){
    
}

/** @brief Set autopointer of the RawData */
void TRMSimulationData::SetRawData(TPZAutoPointer<TRMRawData> &RawData){
    fRawData = RawData;
    SetSystemType(RawData->fSystemType);
    SetTimeControls(RawData->fn_steps, RawData->fdt, RawData->fdt_up, RawData->fdt_up);
    SetNumericControls(RawData->fn_corrections, RawData->fepsilon_res, RawData->fepsilon_cor);
}

/** @brief Setup reporting times and time step size */
void TRMSimulationData::SetTimeControls(int n_times, STATE dt, STATE dt_in, STATE dt_de){
    fn_steps = n_times;
    fReportingTimes.Resize(fn_steps, 0.0);
    for (int it = 0 ; it < fn_steps; it++) {
        fReportingTimes[it] = REAL(it+1)*dt;
    }
    ftime_0    = fReportingTimes[0];
    ftime_n    = fReportingTimes[fn_steps-1];
}

/** @brief Setup reporting times and time step size */
void TRMSimulationData::SetNumericControls(int n_corrections, STATE epsilon_res, STATE epsilon_cor){
    fn_corrections  = n_corrections;
    fepsilon_res    = epsilon_res;
    fepsilon_cor    = epsilon_cor;
}



