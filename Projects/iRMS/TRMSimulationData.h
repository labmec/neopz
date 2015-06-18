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
#include "TRMRawData.h"


class TRMSimulationData {
    
protected:
    
    /** @brief Stores all the petrophysics data */
    TRMPetrophysicsProperties fPetroProp;
    
    /** @brief Stores all the blackoil fluid properties */
    TRMBlackOilProperties fBlackOilProp;
    
    /** @brief Stores all the rock and geomechanic properties */
    TRMRockProperties fRockProp;
    
    /** @brief Stores the spatial information given in maps */
    TRMSpatialPropertiesMap fSpatialProp;
    
    
public:
    
    
    /** @brief Initialize the raw data */
    TRMSimulationData(TRMRawData &RawData);
    
    /** @brief destructor */
    ~TRMSimulationData();
    
    /**
     * @defgroup Water Properties of water from few field parameters
     * @brief    Implements several set/get attributes for the simulation data:
     *           TRMPetrophysicsProperties
     *           TRMBlackOilProperties
     *           TRMSpatialPropertiesMap
     *
     * @{
     */
    
    /** @brief Stores all the petrophysics data */
    void SetPetroProp(TRMPetrophysicsProperties PetroProp)
    {
        fPetroProp = PetroProp;
    }
    
    TRMPetrophysicsProperties & GetPetroProp()
    {
        return fPetroProp;
    }

    /** @brief Stores all the blackoil fluid properties */
    void SetBlackOilProp(TRMBlackOilProperties BlackOilProp)
    {
        fBlackOilProp = BlackOilProp;
    }
    
    TRMBlackOilProperties & GetBlackOilProp()
    {
        return fBlackOilProp;
    }
  
    /** @brief Stores all the rock and geomechanic properties */
    void SetRockProp(TRMRockProperties RockProp)
    {
        fRockProp = RockProp;
    }
    
    TRMRockProperties & GetRockProp()
    {
        return fRockProp;
    }
    
    /** @brief Stores the spatial information given in maps */
    void SetSpatialProp(TRMSpatialPropertiesMap SpatialProp)
    {
        fSpatialProp = SpatialProp;
    }
    
    TRMSpatialPropertiesMap & GetSpatialProp()
    {
        return fSpatialProp;
    }
    
    static REAL Gravity();
    
    // @}
    

};

#endif /* defined(__PZ__TRMSimulationData__) */
