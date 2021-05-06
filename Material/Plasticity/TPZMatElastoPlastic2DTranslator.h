/* 
 * File:   TPZMatElastoPlastic2DTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 17:02
 */

#ifndef TPZMATELASTOPLASTIC2DTRANSLATOR_H
#define TPZMATELASTOPLASTIC2DTRANSLATOR_H

#include "TPZElastoPlasticMemTranslator.h"
#include "TPZMatElastoPlasticTranslator.h"

template <typename TTranslator, typename TMEMTranslator = TPZElastoPlasticMemTranslator>
class TPZMatElastoPlastic2DTranslator : public TPZMatElastoPlasticTranslator<TTranslator,TMEMTranslator> {
public:
    TPZMatElastoPlastic2DTranslator();
    TPZMatElastoPlastic2DTranslator(const TPZMatElastoPlastic2DTranslator<TTranslator, TMEMTranslator>& orig);
    
    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);

    virtual ~TPZMatElastoPlastic2DTranslator();
private:
TPZMatElastoPlasticTranslator<TTranslator,TMEMTranslator> parentTranslator;
};


template <typename TTranslator, typename TMEMTranslator>
TPZMatElastoPlastic2DTranslator<TTranslator, TMEMTranslator>::TPZMatElastoPlastic2DTranslator() {
}

template <typename TTranslator, typename TMEMTranslator>
TPZMatElastoPlastic2DTranslator<TTranslator, TMEMTranslator>::TPZMatElastoPlastic2DTranslator(const TPZMatElastoPlastic2DTranslator<TTranslator, TMEMTranslator>& orig) {
}

template <typename TTranslator, typename TMEMTranslator>
void TPZMatElastoPlastic2DTranslator<TTranslator, TMEMTranslator>::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion){
    parentTranslator.UpdateStream(chunk, toVersion);
    int classid;
    chunk.mOldStream.Read(&classid);
    chunk.mNewStream.Write(&classid);
}

template <typename TTranslator, typename TMEMTranslator>
TPZMatElastoPlastic2DTranslator<TTranslator, TMEMTranslator>::~TPZMatElastoPlastic2DTranslator() {
}



#endif /* TPZMATELASTOPLASTIC2DTRANSLATOR_H */

