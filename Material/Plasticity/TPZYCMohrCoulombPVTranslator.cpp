/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZYCMohrCoulombPVTranslator.cpp
 * Author: thiago
 * 
 * Created on 12 de Março de 2018, 15:35
 */

#include "TPZYCMohrCoulombPVTranslator.h"
#include "TPZChunkInTranslation.h"

TPZYCMohrCoulombPVTranslator::TPZYCMohrCoulombPVTranslator() {
}

TPZYCMohrCoulombPVTranslator::TPZYCMohrCoulombPVTranslator(const TPZYCMohrCoulombPVTranslator& orig) {
}

void TPZYCMohrCoulombPVTranslator::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    REAL fPhi;
    chunk.mOldStream.Read(&fPhi);
    chunk.mNewStream.Write(&fPhi);
    REAL fPsi;
    chunk.mOldStream.Read(&fPsi);
    chunk.mNewStream.Write(&fPsi);
    REAL fc;
    chunk.mOldStream.Read(&fc);
    chunk.mNewStream.Write(&fc);
    REAL fEpsPlasticBar;
    chunk.mOldStream.Read(&fEpsPlasticBar);
    chunk.mNewStream.Write(&fEpsPlasticBar);
    tpzElasticResponseTranslator.UpdateStream(chunk, toVersion);
}

TPZYCMohrCoulombPVTranslator::~TPZYCMohrCoulombPVTranslator() {
}

