//
//  TRMRawData.h
//  PZ
//
//  Created by Philippe Devloo on 5/25/15.
//
// This class store the computational information required for iRMS

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
    
    TRMRawData();
    ~TRMRawData();
};

#endif /* defined(__PZ__TRMRawData__) */
