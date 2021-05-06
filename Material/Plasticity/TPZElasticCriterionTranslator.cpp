/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZElasticCriterionTranslator.cpp
 * Author: thiago
 * 
 * Created on 12 de Março de 2018, 15:22
 */

#include "TPZElasticCriterionTranslator.h"

TPZElasticCriterionTranslator::TPZElasticCriterionTranslator() {
}

TPZElasticCriterionTranslator::TPZElasticCriterionTranslator(const TPZElasticCriterionTranslator& orig) {
}

void TPZElasticCriterionTranslator::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    tpzPlasticStateTranslatorSTATE.UpdateStream(chunk, toVersion);
    tpzElasticResponseTranslator.UpdateStream(chunk, toVersion);
}

TPZElasticCriterionTranslator::~TPZElasticCriterionTranslator() {
}

