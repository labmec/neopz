//$Id: pzelastoplastic.cpp,v 1.33 2010-10-18 15:37:59 diogo Exp $

#include "TPZMatElastoPlastic_impl.h"

#include "TPZLadeKim.h"  
#include "TPZSandlerDimaggio.h"
#include "TPZYCDruckerPrager.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "TPZYCMohrCoulomb.h"
#include "TPZMohrCoulomb.h"
#include "TPZDruckerPrager.h"
#include "TPZYCWillamWarnke.h"
#include "TPZWillamWarnke.h"
#include "TPZVonMises.h"
#include "TPZYCVonMises.h"
#include "TPZYCModifiedMohrCoulomb.h"
#include "TPZSandlerExtended.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZElasticCriterion.h"
#include "TPZYCCamClayPV.h"

template class TPZMatElastoPlastic<TPZPlasticStep<TPZYCModifiedMohrCoulomb, TPZThermoForceA, TPZElasticResponse>, TPZElastoPlasticMem>;
//template class TPZMatElastoPlastic<TPZModifiedMohrCoulomb>;

template class TPZMatElastoPlastic<TPZPlasticStep<TPZYCWillamWarnke, TPZThermoForceA, TPZElasticResponse> , TPZElastoPlasticMem>;
template class TPZMatElastoPlastic<TPZWillamWarnke>;


template class TPZMatElastoPlastic<TPZLadeKim, TPZElastoPlasticMem>;
template class TPZMatElastoPlastic<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>, TPZElastoPlasticMem>;
template class TPZMatElastoPlastic<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>, TPZElastoPlasticMem>;


template class TPZMatElastoPlastic<TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> , TPZElastoPlasticMem>;
template class TPZMatElastoPlastic<TPZDruckerPrager>;


template class TPZMatElastoPlastic<TPZPlasticStep<TPZYCMohrCoulomb, TPZThermoForceA, TPZElasticResponse>, TPZElastoPlasticMem>;
template class TPZMatElastoPlastic<TPZMohrCoulomb>;

template class TPZMatElastoPlastic<TPZPlasticStep<TPZYCVonMises, TPZThermoForceA, TPZElasticResponse>, TPZElastoPlasticMem>;
template class TPZMatElastoPlastic<TPZVonMises>;

template class TPZMatElastoPlastic<TPZLadeKim, TPZPoroElastoPlasticMem>;
template class TPZMatElastoPlastic<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>, TPZPoroElastoPlasticMem>;
template class TPZMatElastoPlastic<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>, TPZPoroElastoPlasticMem>;
template class TPZMatElastoPlastic<TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> , TPZPoroElastoPlasticMem>;

template class TPZMatElastoPlastic<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem>;
template class TPZMatElastoPlastic<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem>;
template class TPZMatElastoPlastic<TPZPlasticStepPV<TPZYCCamClayPV,TPZElasticResponse> , TPZElastoPlasticMem>;

template class TPZMatElastoPlastic<TPZElasticCriterion , TPZElastoPlasticMem>;
template class TPZMatElastoPlastic<TPZElasticCriterion , TPZPoroElastoPlasticMem>;


