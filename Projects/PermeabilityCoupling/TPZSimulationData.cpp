//
//  TPZSimulationData.cpp
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#include "TPZSimulationData.h"


/** @brief Initialize the raw data */
TPZSimulationData::TPZSimulationData(){
    
    /** @brief initial state */
    fIsInitialStateQ = false;
    
    /** @brief current time state */
    fIsCurrentStateQ =  false;
    
    /** @brief Definition gravity field */
    fg.Resize(0);
    
    /** @brief Store time values to be reported */
    fReportingTimes.Resize(0);
    
    /** @brief ntime steps */
    fn_steps = 0;
    
    /** @brief Time step */
    fdt = 0.0;
    
    /** @brief Time step */
    ftime = 0.0;
    
}

/** @brief destructor */
TPZSimulationData::~TPZSimulationData(){
    
}

/** @brief Setup reporting times and time step size */
void TPZSimulationData::SetTimeControls(int n_times, REAL dt){
    
    fn_steps    = n_times;
    fdt         = dt;
    fReportingTimes.Resize(n_times, 0.0);
    for (int it = 0; it < n_times; it++) {
        fReportingTimes[it] = it*dt;
    }
    
}

/** @brief Setup reporting times and time step size */
void TPZSimulationData::SetNumericControls(int n_corrections, REAL epsilon_res, REAL epsilon_cor){
    
    fn_corrections  =   n_corrections;
    fepsilon_res    =   epsilon_res;
    fepsilon_cor    =   epsilon_cor;
    
}