//
//  TPZNullMaterialTranslator.cpp
//  pz
//
//  Created by Omar Durán on 9/26/19.
//

#include "TPZNullMaterialTranslator.h"
#include "TPZChunkInTranslation.h"

TPZNullMaterialTranslator::TPZNullMaterialTranslator() {
}

TPZNullMaterialTranslator::TPZNullMaterialTranslator(const TPZNullMaterialTranslator & other) {
    
}

void TPZNullMaterialTranslator::UpdateAttributes(TPZChunkInTranslation & chunk, const std::map<std::string, uint64_t>& toVersion) {
    
    parentTranslator.UpdateAttributes(chunk, toVersion);
    int fDim;
    chunk.mOldStream.Read(&fDim);
    chunk.mNewStream.Write(&fDim);
    
    int fNState;
    chunk.mOldStream.Read(&fNState);
    chunk.mNewStream.Write(&fNState);
}


TPZNullMaterialTranslator::~TPZNullMaterialTranslator() {
    
}

