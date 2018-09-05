/* 
 * File:   TPZYCSandlerDimaggioTranslator.cpp
 * Author: thiago
 * 
 * Created on 12 de Março de 2018, 19:41
 */

#include "TPZYCSandlerDimaggioTranslator.h"
#include "TPZChunkInTranslation.h"

TPZYCSandlerDimaggioTranslator::TPZYCSandlerDimaggioTranslator() {
}

TPZYCSandlerDimaggioTranslator::TPZYCSandlerDimaggioTranslator(const TPZYCSandlerDimaggioTranslator& orig) {
}

void TPZYCSandlerDimaggioTranslator::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    REAL fA;
    chunk.mOldStream.Read(&fA);
    chunk.mNewStream.Write(&fA);
    REAL fB;
    chunk.mOldStream.Read(&fB);
    chunk.mNewStream.Write(&fB);
    REAL fC;
    chunk.mOldStream.Read(&fC);
    chunk.mNewStream.Write(&fC);
    REAL fD;
    chunk.mOldStream.Read(&fD);
    chunk.mNewStream.Write(&fD);
    REAL fW;
    chunk.mOldStream.Read(&fW);
    chunk.mNewStream.Write(&fW);
    REAL fR;
    chunk.mOldStream.Read(&fR);
    chunk.mNewStream.Write(&fR);
}

TPZYCSandlerDimaggioTranslator::~TPZYCSandlerDimaggioTranslator() {
}

