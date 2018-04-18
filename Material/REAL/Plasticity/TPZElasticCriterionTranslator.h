/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZElasticCriterionTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 15:22
 */

#ifndef TPZELASTICCRITERIONTRANSLATOR_H
#define TPZELASTICCRITERIONTRANSLATOR_H

#include "TPZChunkTranslator.h"
#include "TPZPlasticStateTranslator.h"
#include "TPZElasticResponseTranslator.h"

class TPZElasticCriterionTranslator : public TPZChunkTranslator {
public:
    TPZElasticCriterionTranslator();
    TPZElasticCriterionTranslator(const TPZElasticCriterionTranslator& orig);
    
    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    
    virtual ~TPZElasticCriterionTranslator();
private:
    TPZPlasticStateTranslator<STATE> tpzPlasticStateTranslatorSTATE;
    TPZElasticResponseTranslator tpzElasticResponseTranslator;
};

#endif /* TPZELASTICCRITERIONTRANSLATOR_H */

