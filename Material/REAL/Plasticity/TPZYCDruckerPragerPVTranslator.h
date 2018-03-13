/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZYCDruckerPragerPVTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 15:33
 */

#ifndef TPZYCDRUCKERPRAGERPVTRANSLATOR_H
#define TPZYCDRUCKERPRAGERPVTRANSLATOR_H

#include "TPZChunkTranslator.h"
#include "TPZYCCamClayPVTranslator.h"

class TPZYCDruckerPragerPVTranslator : public TPZChunkTranslator {
public:
    TPZYCDruckerPragerPVTranslator();
    TPZYCDruckerPragerPVTranslator(const TPZYCDruckerPragerPVTranslator& orig);

    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    
    virtual ~TPZYCDruckerPragerPVTranslator();
private:
    TPZYCCamClayPVTranslator tpzYCCamClayPVTranslator ;
};

#endif /* TPZYCDRUCKERPRAGERPVTRANSLATOR_H */

