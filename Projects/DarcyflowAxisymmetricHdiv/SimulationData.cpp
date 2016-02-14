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
    
    /** @brief Number of threads for assembly */
    fnthreads = 0;
    
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
    
    /** @brief Approximation Order for saturations */
    fsorder = 0;
    
    /** @brief Used a direct Solver */
    fIsDirect = false;
    
    /** @brief Optimize band width */
    fOptband = false;
    
    /** @brief Use Conjugated Gradient method */
    fIsCG = false;
    
    /** @brief Broyden iterations */
    fIsPicard = false;
    
    /** @brief Define the use of static condensation */
    fCondenseElements = false;
    
    /** @brief Define the use dimensionless formulation */
    fIsDimensionless = false;
    
    /** @brief Define the use of a linear approxiamtion of S using a gradient reconstruction procedure */
    fUseGR = false;
    
    /** @brief void material being used for GR */
    fMatL2 = 1000;
    
    /** @brief State: n or n+1 temporal state */
    fnStep = true;
    
    /** @brief Gravity */
    fGravity.Resize(2,1);
    
    /** @brief Definition of the flow system one - two and  ... three phase */
    fSystemType.Resize(0);
    
    /** @brief vector that contains materials id to integrate the production */
    fMaterialsToIntegrate.Resize(0);
    
    /** @brief Is one-phase flow? */
    fIsOnePhaseQ = false;
    
    /** @brief Is two-phase flow? */
    fIsTwoPhaseQ = false;
    
    /** @brief Is three-phase flow? */
    fIsThreePhaseQ = false;
    
    /** @brief Using the Linear grativational segregation model */
    fIsLinearSegregationsQ = true;
    
    /** @brief Is axisymmetric analysis */
    fIsAxisymmetricQ = false;
    
    /** @brief Is Impes analysis */
    fIsImpesQ = false;
    
    /** @brief Is hydrostatic boundary condition */
    fIsHydrostaticBCQ =  false;
    
    /** @brief Counterclockwise rotation angle */
    fAngle = 0.0;
    
    /** @brief Is a mesh with geometric progression */
    fIsMeshwithPGQ = false;
    
    /** @brief Compute a variable K distribution */
    fIsHeterogeneousQ = false;
    
    /** @brief geometric progression ratio >= 1 */
    fpg_ratio = 1.0;
    
    /** @brief Store time values to be reported */
    fTimesToPrint.Resize(0);
    
    /** @brief Ouput directory name */
    foutdirectory = "Dump";
    
    /** @brief GID mesh file */
    fGIDfile = "Batata.dump";

    /** @brief Time scale for dimensionless calculations */
    ftime_scale = 1.0;
    
    /** @brief Length scale for dimensionless calculations */
    flength_scale = 1.0;
    
    /** @brief Velocity scale for dimensionless calculations */
    fvelocity_scale = 1.0;
    
}

SimulationData::~SimulationData()
{
    
}