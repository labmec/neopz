/* 
 * File:   TPZMatElastoPlasticTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 18:35
 */

#ifndef TPZMATELASTOPLASTICTRANSLATOR_H
#define TPZMATELASTOPLASTICTRANSLATOR_H

#include "TPZMatWithMemTranslator.h"
#include "TPZElastoPlasticMemTranslator.h"

template <class TTranslator, class TMEMTranslator = TPZElastoPlasticMemTranslator>
class TPZMatElastoPlasticTranslator : public TPZMatWithMemTranslator<TMEMTranslator>{
public:
    TPZMatElastoPlasticTranslator();
    TPZMatElastoPlasticTranslator(const TPZMatElastoPlasticTranslator<TTranslator, TMEMTranslator>& orig);
    
    virtual void UpdateStream(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) override;

    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) override;

    virtual ~TPZMatElastoPlasticTranslator();
private:
    TTranslator tTranslator;
    TPZMatWithMemTranslator<TMEMTranslator> parentTranslator;
};


template <class TTranslator, class TMEMTranslator>
TPZMatElastoPlasticTranslator<TTranslator, TMEMTranslator>::TPZMatElastoPlasticTranslator() {
}

template <class TTranslator, class TMEMTranslator>
TPZMatElastoPlasticTranslator<TTranslator, TMEMTranslator>::TPZMatElastoPlasticTranslator(const TPZMatElastoPlasticTranslator<TTranslator, TMEMTranslator>& orig) {
}

template <class TTranslator, class TMEMTranslator>
void TPZMatElastoPlasticTranslator<TTranslator, TMEMTranslator>::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion){
    parentTranslator.UpdateStream(chunk, toVersion);
    TPZManVector<REAL, 3> fForce(3);    
    chunk.mOldStream.Read(&fForce[0], 3);
    chunk.mNewStream.Write(&fForce[0], 3);
    TPZManVector<REAL,3> fPostProcessDirection(3);
    chunk.mOldStream.Read(&fPostProcessDirection[0], 3);
    chunk.mNewStream.Write(&fPostProcessDirection[0], 3);
    tTranslator.UpdateStream(chunk, toVersion);
    REAL fTol;
    chunk.mOldStream.Read(&fTol);
    chunk.mNewStream.Write(&fTol);
}


template <class TTranslator, class TMEMTranslator>
TPZMatElastoPlasticTranslator<TTranslator, TMEMTranslator>::~TPZMatElastoPlasticTranslator() {
}

template <class TTranslator, class TMEMTranslator>
void TPZMatElastoPlasticTranslator<TTranslator, TMEMTranslator>::UpdateStream(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    this->UpdateAttributes(chunk, toVersion);
}


#endif /* TPZMATELASTOPLASTICTRANSLATOR_H */

