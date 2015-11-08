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
    
    fPressure           = 0.0;
    fPressure_n         = 0.0;
    fSw                 = 0.0;
    fSw_n               = 0.0;
    fporosity           = 0.0;
    
    K.Resize(3,3);
    K.Zero();
    K.Identity();

}

/** @brief Default destructor */
TRMMemory::~TRMMemory(){
    
}