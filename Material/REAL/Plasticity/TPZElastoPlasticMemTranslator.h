/* 
 * File:   TPZElastoPlasticMemTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 18:58
 */

#ifndef TPZELASTOPLASTICMEMTRANSLATOR_H
#define TPZELASTOPLASTICMEMTRANSLATOR_H

#include "TPZChunkTranslator.h"
#include "TPZPlasticStateTranslator.h"

class TPZElastoPlasticMemTranslator : public TPZChunkTranslator {
public:
    TPZElastoPlasticMemTranslator();
    TPZElastoPlasticMemTranslator(const TPZElastoPlasticMemTranslator& orig);
  

    virtual void UpdateStream(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);

    virtual ~TPZElastoPlasticMemTranslator();
private:
    
    void UpdateFromV1(TPZChunkInTranslation &chunk, const std::map<std::string, uint64_t> &toVersion);
    TPZPlasticStateTranslator<REAL> tpzPlasticStateTranslatorREAL;
    TPZTensorTranslator<REAL> tpzTensorTranslatorREAL;
};

#endif /* TPZELASTOPLASTICMEMTRANSLATOR_H */

