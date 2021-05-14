/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZYCCamClayPVTranslator.cpp
 * Author: thiago
 * 
 * Created on 12 de Março de 2018, 15:30
 */

#include "TPZYCCamClayPVTranslator.h"
#include "TPZChunkInTranslation.h"

TPZYCCamClayPVTranslator::TPZYCCamClayPVTranslator() {
}

TPZYCCamClayPVTranslator::TPZYCCamClayPVTranslator(const TPZYCCamClayPVTranslator& orig) {
}

void TPZYCCamClayPVTranslator::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    REAL fGamma;
    chunk.mOldStream.Read(&fGamma);
    chunk.mNewStream.Write(&fGamma);
    REAL fM;
    chunk.mOldStream.Read(&fM);
    chunk.mNewStream.Write(&fM);
    REAL fPt;
    chunk.mOldStream.Read(&fPt);
    chunk.mNewStream.Write(&fPt);
    REAL fLogHardening;
    chunk.mOldStream.Read(&fLogHardening);
    chunk.mNewStream.Write(&fLogHardening);
    REAL fLogBulkModulus; 
    chunk.mOldStream.Read(&fLogBulkModulus);
    chunk.mNewStream.Write(&fLogBulkModulus);
    REAL fA0; 
    chunk.mOldStream.Read(&fA0);
    chunk.mNewStream.Write(&fA0);
    REAL fE0;     
    chunk.mOldStream.Read(&fE0);
    chunk.mNewStream.Write(&fE0);
    tpzElasticResponseTranslator.UpdateStream(chunk, toVersion);
}

TPZYCCamClayPVTranslator::~TPZYCCamClayPVTranslator() {
}

