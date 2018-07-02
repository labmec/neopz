/* 
 * File:   TPZPlasticStepTranslator.cpp
 * Author: thiago
 * 
 * Created on 12 de Março de 2018, 19:47
 */

#include "TPZPlasticStepTranslator.h"
#include "TPZYCSandlerDimaggioTranslator.h"
#include "TPZSandlerDimaggioThermoForceATranslator.h"
#include "TPZElasticResponseTranslator.h"
#include "TPZYCSandlerDimaggioLTranslator.h"
#include "TPZYCSandlerDimaggioL2Translator.h"


template class TPZPlasticStepTranslator<TPZYCSandlerDimaggioTranslator, TPZSandlerDimaggioThermoForceATranslator, TPZElasticResponseTranslator>;

template class TPZPlasticStepTranslator<TPZYCSandlerDimaggioLTranslator, TPZSandlerDimaggioThermoForceATranslator, TPZElasticResponseTranslator>;

template class TPZPlasticStepTranslator<TPZYCSandlerDimaggioL2Translator, TPZSandlerDimaggioThermoForceATranslator, TPZElasticResponseTranslator>;