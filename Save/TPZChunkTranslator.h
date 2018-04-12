/* 
 * File:   TPZChunkTranslator.h
 * Author: thiago
 *
 * Created on 7 de Março de 2018, 17:14
 */

#ifndef TPZCHUNKTRANSLATOR_H
#define TPZCHUNKTRANSLATOR_H

#include "TPZChunkInTranslation.h"


class TPZChunkTranslator {
    
public:
    virtual void UpdateStream(TPZChunkInTranslation &chunk, const std::map<std::string, uint64_t> &toVersion) {
        UpdateAttributes(chunk, toVersion);
    }
    
    virtual void UpdateAttributes(TPZChunkInTranslation &chunk, const std::map<std::string, uint64_t> &toVersion)=0;
    
};

#endif /* TPZCHUNKTRANSLATOR_H */

