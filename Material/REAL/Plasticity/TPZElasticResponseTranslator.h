/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZElasticResponseTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 15:26
 */

#ifndef TPZELASTICRESPONSETRANSLATOR_H
#define TPZELASTICRESPONSETRANSLATOR_H

#include "TPZChunkTranslator.h"

class TPZElasticResponseTranslator : public TPZChunkTranslator {
    
public:
    
    TPZElasticResponseTranslator();
    
    TPZElasticResponseTranslator(const TPZElasticResponseTranslator& orig);
    
    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    
    virtual ~TPZElasticResponseTranslator();

};

#endif /* TPZELASTICRESPONSETRANSLATOR_H */

