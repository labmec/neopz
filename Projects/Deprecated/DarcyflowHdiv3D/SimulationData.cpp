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
    
    /** @brief Max Time - s */
    fMaxTime = 0.0;    
    
    /** @brief Time - s */
    fTime = 0.0;
    
    /** @brief Maximum number of newton iterations */
    fMaxiterations = 0;
    
    /** @brief Maximum number of newton iterations */
    fFixedJacobianIterations = 0;
    
    /** @brief DeltaX tolerance for newton iterations */
    ftoleranceDeltaX = 0.0;
    
    /** @brief Residual tolerance for newton iterations */
    ftoleranceResiual = 0.0;
    
    /** @brief Number of uniform mesh refinement */
    fHref = 0;
    
    /** @brief Number of uniform mesh refinement in postprocessing */
    fHrefpost = 0;
    
    /** @brief Approximation Order for velocity */
    fqorder = 1;
    
    /** @brief Approximation Order for pressure */
    fporder = 1;
    
    /** @brief Used a direct Solver */
    fIsDirect = false;
    
    /** @brief Use Conjugated Gradient method */
    fIsCG = false;
    
    /** @brief Broyden iterations */
    fIsBroyden = false;
    
    /** @brief Computes a H1 approximation */
    fIsH1approx = false;
    
    /** @brief State: n or n+1 temporal state */
    fnStep = true;
    
}

SimulationData::~SimulationData()
{
    
}

SimulationData::SimulationData(const SimulationData &copy){
    
}