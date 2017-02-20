//
//  TRMPhaseMemory.cpp
//  PZ
//
//  Created by Philippe Devloo on 7/5/15.
//
//

#include "TRMPhaseMemory.h"


/** @brief Default constructor */
TRMPhaseMemory::TRMPhaseMemory(){

    fun = 0.0;
    fp_avg      = 0.0;
    fp_avg_n    = 0.0;
    fsa         = 0.0;
    fsa_n       = 0.0;
    fsb         = 0.0;
    fsb_n       = 0.0;
    
    fporosity           = 0.0;
    
    fK.Resize(3,3);
    fK.Zero();
    fK.Identity();
    
    fKinv.Resize(3,3);
    fKinv.Zero();
    fKinv.Identity();
    
    fx.Resize(3,0.0);
    
}

/** @brief Default destructor */
TRMPhaseMemory::~TRMPhaseMemory(){
    
}

