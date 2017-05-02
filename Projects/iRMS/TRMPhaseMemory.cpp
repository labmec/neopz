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

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Memory (\partial Gamma and Omega) :: Spatial memory items
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief first lame undrained parameter */
    f_lambda = 0.0;
    
    /** @brief second lame undrained parameter */
    f_mu = 0.0;
    
    /** @brief Storage constrained modulus */
    f_S_e = 0.0;
    
    /** @brief Maurice Biot coefficient */
    f_alpha = 0.0;
    
    fun = 0.0;
    fp_avg      = 0.0;
    fp_avg_n    = 0.0;
    fsa_0       = 0.0;
    fsb_0       = 0.0;
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

