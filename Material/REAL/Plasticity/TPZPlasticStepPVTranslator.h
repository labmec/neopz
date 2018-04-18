/* 
 * File:   TPZPlasticStepPVTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 15:46
 */

#ifndef TPZPLASTICSTEPPVTRANSLATOR_H
#define TPZPLASTICSTEPPVTRANSLATOR_H

#include "TPZChunkTranslator.h"
#include "TPZPlasticStateTranslator.h"

template <typename YC_Translator, typename ER_Translator>
class TPZPlasticStepPVTranslator : public TPZChunkTranslator {
public:
    TPZPlasticStepPVTranslator();
    TPZPlasticStepPVTranslator(const TPZPlasticStepPVTranslator<YC_Translator, ER_Translator>& orig);
    
    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
   
    virtual ~TPZPlasticStepPVTranslator();
private:
    YC_Translator yc_Translator;
    ER_Translator er_Translator;
    TPZPlasticStateTranslator<STATE> tpzPlasticStateTranslatorSTATE;
};

template <typename YC_Translator, typename ER_Translator>
TPZPlasticStepPVTranslator<YC_Translator, ER_Translator>::TPZPlasticStepPVTranslator(){
    
}

template <typename YC_Translator, typename ER_Translator>
TPZPlasticStepPVTranslator<YC_Translator, ER_Translator>::TPZPlasticStepPVTranslator(const TPZPlasticStepPVTranslator<YC_Translator, ER_Translator>& orig){
    
}

template <typename YC_Translator, typename ER_Translator>
void TPZPlasticStepPVTranslator<YC_Translator, ER_Translator>::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion){
    yc_Translator.UpdateStream(chunk, toVersion);
    er_Translator.UpdateStream(chunk, toVersion);
    REAL fResTol;
    chunk.mOldStream.Read(&fResTol);
    chunk.mNewStream.Write(&fResTol);
    int fMaxNewton;
    chunk.mOldStream.Read(&fMaxNewton);
    chunk.mNewStream.Write(&fMaxNewton);
    tpzPlasticStateTranslatorSTATE.UpdateStream(chunk, toVersion);
}

template <typename YC_Translator, typename ER_Translator>
TPZPlasticStepPVTranslator<YC_Translator, ER_Translator>::~TPZPlasticStepPVTranslator(){
}

#endif /* TPZPLASTICSTEPPVTRANSLATOR_H */

