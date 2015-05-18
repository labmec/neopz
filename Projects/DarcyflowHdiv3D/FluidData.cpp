/*
 *  ReservoirData.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "FluidData.h"

FluidData::FluidData()
{
    /** @brief Temperature @ reservoir conditions  - F */
    fReservoirTemperature = 180.0;
    
    /** @brief Characteristic Pressure - Pa */
    fPref = 1.0;
    
    /** @brief Characteristic Density - kg/m3 */
    fRhoref = 1.0;
    
    /** @brief Characteristic viscosity - Pa s */
    fMuref = 1.0;
    
}

FluidData::~FluidData()
{
    
}



