/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZYCMohrCoulombPVTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 15:35
 */

#ifndef TPZYCMOHRCOULOMBPVTRANSLATOR_H
#define TPZYCMOHRCOULOMBPVTRANSLATOR_H

#include "TPZChunkTranslator.h"
#include "TPZElasticResponseTranslator.h"

class TPZYCMohrCoulombPVTranslator : public TPZChunkTranslator {
public:
    TPZYCMohrCoulombPVTranslator();
    TPZYCMohrCoulombPVTranslator(const TPZYCMohrCoulombPVTranslator& orig);

    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    
    virtual ~TPZYCMohrCoulombPVTranslator();
private:
    TPZElasticResponseTranslator tpzElasticResponseTranslator;
};

#endif /* TPZYCMOHRCOULOMBPVTRANSLATOR_H */

