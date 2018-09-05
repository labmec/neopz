/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZElasticResponseTranslator.cpp
 * Author: thiago
 * 
 * Created on 12 de Março de 2018, 15:26
 */

#include "TPZElasticResponseTranslator.h"
#include "pzreal.h"
#include "TPZChunkInTranslation.h"

TPZElasticResponseTranslator::TPZElasticResponseTranslator() {
}

TPZElasticResponseTranslator::TPZElasticResponseTranslator(const TPZElasticResponseTranslator& orig) {
}

void TPZElasticResponseTranslator::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    REAL fLambda;
    chunk.mOldStream.Read(&fLambda);
    chunk.mNewStream.Write(&fLambda);
    REAL fMu;
    chunk.mOldStream.Read(&fMu);
    chunk.mNewStream.Write(&fMu);
}


TPZElasticResponseTranslator::~TPZElasticResponseTranslator() {
}

