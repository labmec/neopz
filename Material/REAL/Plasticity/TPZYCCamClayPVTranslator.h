/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZYCCamClayPVTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 15:30
 */

#ifndef TPZYCCAMCLAYPVTRANSLATOR_H
#define TPZYCCAMCLAYPVTRANSLATOR_H

#include "TPZChunkTranslator.h"
#include "TPZElasticResponseTranslator.h"

class TPZYCCamClayPVTranslator : public TPZChunkTranslator {
public:
    TPZYCCamClayPVTranslator();
    TPZYCCamClayPVTranslator(const TPZYCCamClayPVTranslator& orig);
    
    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    
    virtual ~TPZYCCamClayPVTranslator();
private:
    TPZElasticResponseTranslator tpzElasticResponseTranslator;
};

#endif /* TPZYCCAMCLAYPVTRANSLATOR_H */

