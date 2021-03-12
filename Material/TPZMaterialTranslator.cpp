/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZMaterialTranslator.cpp
 * Author: thiago
 * 
 * Created on 12 de Março de 2018, 18:51
 */

#include "TPZMaterialTranslator.h"
#include "TPZChunkInTranslation.h"

TPZMaterialTranslator::TPZMaterialTranslator() {
}

TPZMaterialTranslator::TPZMaterialTranslator(const TPZMaterialTranslator& orig) {
}

void TPZMaterialTranslator::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    int fId;
    chunk.mOldStream.Read(&fId, 1);
    chunk.mNewStream.Write(&fId, 1);
    REAL gBigNumber;
    chunk.mOldStream.Read(&gBigNumber, 1);
    chunk.mNewStream.Write(&gBigNumber, 1);
    TPZPersistenceManager::TranslateNextPointer(chunk, toVersion); // fForcingFunction 
    TPZPersistenceManager::TranslateNextPointer(chunk, toVersion); // fExactSol 
    TPZPersistenceManager::TranslateNextPointer(chunk, toVersion); // fTimeDependentForcingFunction 
    TPZPersistenceManager::TranslateNextPointer(chunk, toVersion); // fTimedependentFunctionExact 
    TPZPersistenceManager::TranslateNextPointer(chunk, toVersion); // fBCForcingFunction 
    TPZPersistenceManager::TranslateNextPointer(chunk, toVersion); // fTimedependentBCForcingFunction 
    bool fLinearContext;
    chunk.mOldStream.Read(fLinearContext);
    chunk.mNewStream.Write(fLinearContext);
    int fNumLoadCases;
    chunk.mOldStream.Read(&fNumLoadCases);
    chunk.mNewStream.Write(&fNumLoadCases);
    int fPostProcIndex;
    chunk.mOldStream.Read(&fPostProcIndex);
    chunk.mNewStream.Write(&fPostProcIndex);
}

TPZMaterialTranslator::~TPZMaterialTranslator() {
}

