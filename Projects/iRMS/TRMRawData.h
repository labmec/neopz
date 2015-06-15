//
//  TRMRawData.h
//  PZ
//
//  Created by Philippe Devloo on 5/25/15.
//
//

#ifndef __PZ__TRMRawData__
#define __PZ__TRMRawData__

#include <stdio.h>
#include "pzreal.h"
class TRMRawData {
    
public:
    
    REAL fLw;
    bool fHasLiner;
    bool fHasCasing;
    
    REAL fReservoirWidth;
    REAL fReservoirLength;
    REAL fReservoirHeight;
    REAL fProdVertPosition;
    
    TRMRawData()
    {
        fLw = 0.;
        fHasLiner = true;
        fHasCasing = true;
        
        fReservoirWidth = 0.;
        fReservoirLength = 0.;
        fReservoirHeight = 0.;
        fProdVertPosition = 0.;
    }
    
    ~TRMRawData(){}
};

#endif /* defined(__PZ__TRMRawData__) */
