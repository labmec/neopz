//
//  TPZPorousElasticCriterionTranslator.cpp
//  pz
//
//  Created by Omar Durán on 9/5/19.
//

#include "TPZPorousElasticCriterionTranslator.h"


TPZPorousElasticCriterionTranslator::TPZPorousElasticCriterionTranslator() {
    
}

TPZPorousElasticCriterionTranslator::TPZPorousElasticCriterionTranslator(const TPZPorousElasticCriterionTranslator & other) {
    
}

void TPZPorousElasticCriterionTranslator::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    tpzPlasticStateTranslatorSTATE.UpdateStream(chunk, toVersion);
    tpzElasticResponseTranslator.UpdateStream(chunk, toVersion);
    tpzPorousElasticResponseTranslator.UpdateStream(chunk, toVersion);
}

TPZPorousElasticCriterionTranslator::~TPZPorousElasticCriterionTranslator() {
    
}
