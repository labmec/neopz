//
//  TRMMemory.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/25/15.
//
//

#include "TRMMemory.h"


/** @brief Default constructor */
TRMMemory::TRMMemory(){
    
    fu.Resize(3, 0.0);
    fu_n.Resize(3, 0.0);
    fdivu               = 0.0;
    fdivu_n             = 0.0;
    
    // Required
    fp              = 0.0;
    fp_n            = 0.0;
    fp_avg          = 0.0;
    fp_avg_n        = 0.0;
    fsa             = 0.0;
    fsa_n           = 0.0;
    fsb             = 0.0;
    fsb_n           = 0.0;

    
    
    fporosity           = 0.0;
    
    fK.Resize(3,3);
    fK.Zero();
    fK.Identity();
    
    fx.Resize(3,0.0);
    fw = 0.0;
    fdet = 0.0;
    frhs = 0.0;
    
}

/** @brief Default destructor */
TRMMemory::~TRMMemory(){
    
}