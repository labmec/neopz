/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZSandlerExtendedTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 15:39
 */

#ifndef TPZSANDLEREXTENDEDTRANSLATOR_H
#define TPZSANDLEREXTENDEDTRANSLATOR_H

#include "TPZChunkTranslator.h"
#include "TPZElasticResponseTranslator.h"

class TPZSandlerExtendedTranslator : public TPZChunkTranslator {
public:
    TPZSandlerExtendedTranslator();
    TPZSandlerExtendedTranslator(const TPZSandlerExtendedTranslator& orig);
    
    void UpdateStream(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    
    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    
    virtual ~TPZSandlerExtendedTranslator();
    
private:
    
    TPZElasticResponseTranslator tpzElasticResponseTranslator;
    
    void UpdateAttributesV1(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    void UpdateFromV2(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    void UpdateAttributesV2(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    
};

#endif /* TPZSANDLEREXTENDEDTRANSLATOR_H */

