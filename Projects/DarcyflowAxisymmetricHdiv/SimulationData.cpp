/*
 *  problemdata.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "SimulationData.h"

SimulationData::SimulationData()
{
    /** @brief Delta t - s */
    fDeltaT = 0.0;
    
    /** @brief Time - s */
    fTime = 0.0;
    
    /** @brief Maximum number of newton iterations */
    fMaxiterations = 0;
    
    /** @brief DeltaX tolerance for newton iterations */
    ftoleranceDeltaX = 0.0;
    
    /** @brief Residual tolerance for newton iterations */
    ftoleranceResiual = 0.0;
    
    /** @brief Broyden iterations */
    fIsBroyden = false;
	
}

SimulationData::~SimulationData()
{
	
}