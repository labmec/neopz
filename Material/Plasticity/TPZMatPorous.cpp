/**
 * @file
 */


#include "TPZMatPorous_impl.h"

#include "TPZLadeKim.h"  
#include "TPZSandlerDimaggio.h"
#include "TPZMatElastoPlastic.h"
#include "TPZYCDruckerPrager.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"


/*
template class TPZMatPorous< TPZMatElastoPlastic<TPZLadeKim, TPZPoroElastoPlasticMem> >;
template class TPZMatPorous< TPZMatElastoPlastic<TPZSandlerDimaggio, TPZPoroElastoPlasticMem> >;
*/

template class TPZMatPorous< TPZLadeKim, TPZPoroElastoPlasticMem >;
template class TPZMatPorous< TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>, TPZPoroElastoPlasticMem >;
template class TPZMatPorous<TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> , TPZPoroElastoPlasticMem>;
