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
    
    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    
    virtual ~TPZSandlerExtendedTranslator();
private:
    TPZElasticResponseTranslator tpzElasticResponseTranslator;
};

#endif /* TPZSANDLEREXTENDEDTRANSLATOR_H */

