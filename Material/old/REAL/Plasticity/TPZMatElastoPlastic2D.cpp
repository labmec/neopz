///*
// *  pzelastoplastic2D.cpp
// *  ElastoPlasticModels
// *
// *  Created by Diogo Cecilio on 10/25/10.
// *  Copyright 2010 __MyCompanyName__. All rights reserved.
// *
// */

#include "TPZMatElastoPlastic2D_impl.h"

#include "pzbndcond.h"
#include "TPZLadeKim.h"
#include "TPZSandlerDimaggio.h"
#include "TPZYCDruckerPrager.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "TPZElasticCriterion.h"
#include "TPZYCMohrCoulomb.h"
#include "TPZMohrCoulomb.h"
#include "TPZDruckerPrager.h"
#include "TPZYCWillamWarnke.h"
#include "TPZWillamWarnke.h"
#include "TPZVonMises.h"
#include "TPZYCVonMises.h"
#include "TPZYCModifiedMohrCoulomb.h"
#include "TPZYCCamClayPV.h"
#include "TPZMatElastoPlastic2DTranslator.h"
#include "TPZSandlerDimaggioTranslator.h"
#include "TPZPlasticStepPVTranslator.h"
#include "TPZYCMohrCoulombPVTranslator.h"
#include "TPZSandlerExtendedTranslator.h"
#include "TPZYCCamClayPVTranslator.h"


template class TPZMatElastoPlastic2D<TPZPlasticStep<TPZYCModifiedMohrCoulomb, TPZThermoForceA, TPZElasticResponse>, TPZElastoPlasticMem>;
//template class TPZMatElastoPlastic2D<TPZModifiedMohrCoulomb>;

template class TPZMatElastoPlastic2D<TPZPlasticStep<TPZYCWillamWarnke, TPZThermoForceA, TPZElasticResponse> , TPZElastoPlasticMem>;
template class TPZMatElastoPlastic2D<TPZWillamWarnke>;

template class TPZMatElastoPlastic2D<TPZLadeKim, TPZElastoPlasticMem>;
template class TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>, TPZElastoPlasticMem>;
template class TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>, TPZElastoPlasticMem>;


template class TPZMatElastoPlastic2D<TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> , TPZElastoPlasticMem>;
template class TPZMatElastoPlastic2D<TPZDruckerPrager>;


template class TPZMatElastoPlastic2D<TPZPlasticStep<TPZYCMohrCoulomb, TPZThermoForceA, TPZElasticResponse>, TPZElastoPlasticMem>;
template class TPZMatElastoPlastic2D<TPZMohrCoulomb>;

template class TPZMatElastoPlastic2D<TPZPlasticStep<TPZYCVonMises, TPZThermoForceA, TPZElasticResponse>, TPZElastoPlasticMem>;
template class TPZMatElastoPlastic2D<TPZVonMises>;


template class TPZMatElastoPlastic2D<TPZLadeKim, TPZPoroElastoPlasticMem>;
template class TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>, TPZPoroElastoPlasticMem>;
template class TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>, TPZPoroElastoPlasticMem>;
//template class TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>, TPZPoroElastoPlasticMem>;
template class TPZMatElastoPlastic2D<TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> , TPZPoroElastoPlasticMem>;

template class TPZRestoreClassWithTranslator<TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>, TPZElastoPlasticMem>, TPZMatElastoPlastic2DTranslator<TPZSandlerDimaggioTranslator<SANDLERDIMAGGIOSTEP1TRANSLATOR>, TPZElastoPlasticMemTranslator>>;
template class TPZRestoreClassWithTranslator<TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>, TPZElastoPlasticMem>, TPZMatElastoPlastic2DTranslator<TPZSandlerDimaggioTranslator<SANDLERDIMAGGIOSTEP2TRANSLATOR>, TPZElastoPlasticMemTranslator>>;


template class TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem>;
template class TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem>;
template class TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCCamClayPV,TPZElasticResponse> , TPZElastoPlasticMem>;


template class TPZRestoreClassWithTranslator<TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem>, TPZMatElastoPlastic2DTranslator<TPZPlasticStepPVTranslator<TPZYCMohrCoulombPVTranslator,TPZElasticResponseTranslator> , TPZElastoPlasticMemTranslator>>;
template class TPZRestoreClassWithTranslator<TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse>, TPZElastoPlasticMem>, TPZMatElastoPlastic2DTranslator<TPZPlasticStepPVTranslator<TPZSandlerExtendedTranslator,TPZElasticResponseTranslator>, TPZElastoPlasticMemTranslator>>;
template class TPZRestoreClassWithTranslator<TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCCamClayPV,TPZElasticResponse>, TPZElastoPlasticMem>, TPZMatElastoPlastic2DTranslator<TPZPlasticStepPVTranslator<TPZYCCamClayPVTranslator,TPZElasticResponseTranslator>, TPZElastoPlasticMemTranslator>>;

template class TPZMatElastoPlastic2D<TPZElasticCriterion , TPZElastoPlasticMem>;
template class TPZMatElastoPlastic2D<TPZElasticCriterion , TPZPoroElastoPlasticMem>;
