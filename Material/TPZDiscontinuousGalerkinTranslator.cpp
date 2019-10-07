//
//  TPZDiscontinuousGalerkinTranslator.cpp
//  pz
//
//  Created by Omar Durán on 9/26/19.
//

#include "TPZChunkInTranslation.h"
#include "TPZDiscontinuousGalerkinTranslator.h"

TPZDiscontinuousGalerkinTranslator::TPZDiscontinuousGalerkinTranslator() {
}

TPZDiscontinuousGalerkinTranslator::TPZDiscontinuousGalerkinTranslator(const TPZDiscontinuousGalerkinTranslator & other) {
}

void TPZDiscontinuousGalerkinTranslator::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    parentTranslator.UpdateStream(chunk,toVersion);
}

TPZDiscontinuousGalerkinTranslator::~TPZDiscontinuousGalerkinTranslator() {
}
