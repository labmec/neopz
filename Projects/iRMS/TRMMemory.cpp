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
    
    fu.Resize(0, 0.0);
    fu_n.Resize(0, 0.0);
    fdivu               = 0.0;
    fdivu_n             = 0.0;
    
    fPressure           = 0.0;
    fPressure_n         = 0.0;
    fSw                 = 0.0;
    fSw_n               = 0.0;
    fporosity           = 0.0;
    
    K.Resize(3,3);
    K.Zero();
    K.Identity();
    
    fx.Resize(0,0.0);
    fw = 0.0;
    fdet = 0.0;
    frhs = 0.0;
    
}

/** @brief Default destructor */
TRMMemory::~TRMMemory(){
    
}