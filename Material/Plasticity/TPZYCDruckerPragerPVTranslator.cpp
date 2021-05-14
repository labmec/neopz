/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZYCDruckerPragerPVTranslator.cpp
 * Author: thiago
 * 
 * Created on 12 de Março de 2018, 15:33
 */

#include "TPZYCDruckerPragerPVTranslator.h"

TPZYCDruckerPragerPVTranslator::TPZYCDruckerPragerPVTranslator() {
}

TPZYCDruckerPragerPVTranslator::TPZYCDruckerPragerPVTranslator(const TPZYCDruckerPragerPVTranslator& orig) {
}

void TPZYCDruckerPragerPVTranslator::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    tpzYCCamClayPVTranslator.UpdateStream(chunk, toVersion);
}

TPZYCDruckerPragerPVTranslator::~TPZYCDruckerPragerPVTranslator() {
}

