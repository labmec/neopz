/* 
 * File:   TPZPlasticStepTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 19:47
 */

#ifndef TPZPLASTICSTEPTRANSLATOR_H
#define TPZPLASTICSTEPTRANSLATOR_H

#include "TPZChunkTranslator.h"
#include "TPZPlasticStateTranslator.h"

template <typename YC_t, typename TF_t, typename ER_t>
class TPZPlasticStepTranslator : public TPZChunkTranslator {
public:
    TPZPlasticStepTranslator();
    TPZPlasticStepTranslator(const TPZPlasticStepTranslator<YC_t, TF_t, ER_t>& orig);

    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);

    virtual ~TPZPlasticStepTranslator();
protected:
    YC_t fYCTranslator;
    TF_t fTFTranslator;
    ER_t fERTranslator;
    TPZPlasticStateTranslator<REAL> tpzPlasticStateTranslatorREAL;
};

template <typename YC_t, typename TF_t, typename ER_t>
TPZPlasticStepTranslator<YC_t, TF_t, ER_t>::TPZPlasticStepTranslator() {
}

template <typename YC_t, typename TF_t, typename ER_t>
TPZPlasticStepTranslator<YC_t, TF_t, ER_t>::TPZPlasticStepTranslator(const TPZPlasticStepTranslator<YC_t, TF_t, ER_t>& orig) {
}

template <typename YC_t, typename TF_t, typename ER_t>
void TPZPlasticStepTranslator<YC_t, TF_t, ER_t>::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    fYCTranslator.UpdateStream(chunk, toVersion);
    fTFTranslator.UpdateStream(chunk, toVersion);
    fERTranslator.UpdateStream(chunk, toVersion);
    REAL fResTol;
    chunk.mOldStream.Read(&fResTol);
    chunk.mNewStream.Write(&fResTol);
    REAL fIntegrTol;
    chunk.mOldStream.Read(&fIntegrTol);
    chunk.mNewStream.Write(&fIntegrTol);
    int fMaxNewton;
    chunk.mOldStream.Read(&fMaxNewton);
    chunk.mNewStream.Write(&fMaxNewton);
    REAL fMinLambda;
    chunk.mOldStream.Read(&fMinLambda);
    chunk.mNewStream.Write(&fMinLambda);
    REAL fMinStepSize;
    chunk.mOldStream.Read(&fMinStepSize);
    chunk.mNewStream.Write(&fMinStepSize);
    tpzPlasticStateTranslatorREAL.UpdateStream(chunk, toVersion); //fN
    int fMaterialTensionSign;
    chunk.mOldStream.Read(&fMaterialTensionSign);
    chunk.mNewStream.Write(&fMaterialTensionSign);
    int fInterfaceTensionSign;
    chunk.mOldStream.Read(&fInterfaceTensionSign);
    chunk.mNewStream.Write(&fInterfaceTensionSign);

}

template <typename YC_t, typename TF_t, typename ER_t>
TPZPlasticStepTranslator<YC_t, TF_t, ER_t>::~TPZPlasticStepTranslator() {
}



#endif /* TPZPLASTICSTEPTRANSLATOR_H */

