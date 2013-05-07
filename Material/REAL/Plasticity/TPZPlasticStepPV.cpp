/**
 * @file
 */


#include "TPZPlasticStepPV.h"
#include "TPZElasticResponse.h"


#include "TPZYCSandlerDimaggio.h"

#include "pzlog.h"


template class TPZPlasticStepPV<TPZYCSandlerDimaggio, TPZElasticResponse>;



