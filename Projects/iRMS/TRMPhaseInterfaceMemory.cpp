//
//  TRMPhaseMemory.cpp
//  PZ
//
//  Created by Philippe Devloo on 7/5/15.
//
//

#include "TRMPhaseInterfaceMemory.h"

/** @brief Default constructor */
TRMPhaseInterfaceMemory::TRMPhaseInterfaceMemory(){
    
    /** @brief contains the normal flux per surface area */
    fun = 0.0;
    
    /** @brief contains the volumetric left average pressure at last time step */
    fp_avg_n_l = 0.0;
    
    /** @brief contains the volumetric left saturation of alhpa at last time step */
    fsa_n_l = 0.0;
    
    /** @brief contains the volumetric left saturation of alhpa at last time step */
    fsb_n_l = 0.0;
    
    /** @brief contains the volumetric right average pressure at last time step */
    fp_avg_n_r = 0.0;
    
    /** @brief contains the volumetric right saturation of alhpa at last time step */
    fsa_n_r = 0.0;
    
    /** @brief contains the volumetric right saturation of alhpa at last time step */
    fsb_n_r = 0.0;
}

/** @brief Default destructor */
TRMPhaseInterfaceMemory::~TRMPhaseInterfaceMemory(){
    
}