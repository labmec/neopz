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
    
    /** @brief Autopointer of the RawData */
    TPZAutoPointer<TRMRawData> fRawData;
    
    /** @brief Stores all the petrophysics data */
    TRMPetrophysicsProperties fPetroProp;
    
    /** @brief Stores all the blackoil fluid properties */
    TRMBlackOilProperties fBlackOilProp;
    
    /** @brief Stores all the rock and geomechanic properties */
    TRMRockProperties fRockProp;
    
    /** @brief Stores the spatial information given in maps */
    TRMSpatialPropertiesMap fSpatialProp;
    
    
public:
    
    
    /** @brief default constructor */
    TRMSimulationData();
    
    /** @brief destructor */
    ~TRMSimulationData();
    
    /**
     * @defgroup Access methods
     * @brief    Implements several set/get attributes for the simulation data:
     *
     * @{
     */
    
    /** @brief Set autopointer of the RawData */
    void SetRawData(TPZAutoPointer<TRMRawData> RawData){
        fRawData = RawData;
    }
    
    /** @brief Get autopointer of the RawData */
    TPZAutoPointer<TRMRawData>  RawData(){
        return fRawData;
    }
    
    /** @brief Stores all the petrophysics data */
    void SetPetroProp(TRMPetrophysicsProperties PetroProp)
    {
        fPetroProp = PetroProp;
    }
    
    TRMPetrophysicsProperties & PetroProp()
    {
        return fPetroProp;
    }

    /** @brief Stores all the blackoil fluid properties */
    void SetBlackOilProp(TRMBlackOilProperties BlackOilProp)
    {
        fBlackOilProp = BlackOilProp;
    }
    
    TRMBlackOilProperties & BlackOilProp()
    {
        return fBlackOilProp;
    }
  
    /** @brief Stores all the rock and geomechanic properties */
    void SetRockProp(TRMRockProperties RockProp)
    {
        fRockProp = RockProp;
    }
    
    TRMRockProperties & RockProp()
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
