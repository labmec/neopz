//
//  TRMSimulationData.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMSimulationData__
#define __PZ__TRMSimulationData__

#include <stdio.h>
#include "pzreal.h"
#include "TRMBlackOilProperties.h"
#include "TRMRockProperties.h"
#include "TRMPetrophysicsProperties.h"
#include "TRMSpatialPropertiesMap.h"

class TRMRawData;

class TRMSimulationData {
    
protected:
    
    TRMPetrophysicsProperties fPetroProp;
    TRMBlackOilProperties fBlackOilProp;
    TRMRockProperties fRockProp;
    TRMSpatialPropertiesMap fSpatialProp;
    
    
public:
    
    enum EMaterialIds {}; // All this materials ids are related with physical boundary conditions.
    
    // Initializes the simulation data based on a raw data
    TRMSimulationData(TRMRawData &rawData);
    
    // Set get methods
    
    void SetPetroProp(TRMPetrophysicsProperties PetroProp)
    {
        fPetroProp = PetroProp;
    }
    
    TRMPetrophysicsProperties & GetPetroProp()
    {
        return fPetroProp;
    }

    void SetBlackOilProp(TRMBlackOilProperties BlackOilProp)
    {
        fBlackOilProp = BlackOilProp;
    }
    
    TRMBlackOilProperties & GetBlackOilProp()
    {
        return fBlackOilProp;
    }
  
    void SetRockProp(TRMRockProperties RockProp)
    {
        fRockProp = RockProp;
    }
    
    TRMRockProperties & GetRockProp()
    {
        return fRockProp;
    }
    
    void SetSpatialProp(TRMSpatialPropertiesMap SpatialProp)
    {
        fSpatialProp = SpatialProp;
    }
    
    TRMSpatialPropertiesMap & GetSpatialProp()
    {
        return fSpatialProp;
    }
    
    static REAL Gravity();
    

};

#endif /* defined(__PZ__TRMSimulationData__) */
