/* 
 * File:   TPZSandlerDimaggioTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 19:25
 */

#ifndef TPZSANDLERDIMAGGIOTRANSLATOR_H
#define TPZSANDLERDIMAGGIOTRANSLATOR_H

#include "TPZChunkInTranslation.h"

template<typename SANDLERDIMAGGIOPARENTTRANSLATOR>
class TPZSandlerDimaggioTranslator : public SANDLERDIMAGGIOPARENTTRANSLATOR {
public:
    TPZSandlerDimaggioTranslator();

    virtual void UpdateStream(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);

    TPZSandlerDimaggioTranslator(const TPZSandlerDimaggioTranslator<SANDLERDIMAGGIOPARENTTRANSLATOR>& orig);
    virtual ~TPZSandlerDimaggioTranslator();
private:
    void UpdateFromV1(TPZChunkInTranslation &chunk, const std::map<std::string, uint64_t> &toVersion);
    SANDLERDIMAGGIOPARENTTRANSLATOR parentTranslator;
};

template<typename SANDLERDIMAGGIOPARENTTRANSLATOR>
TPZSandlerDimaggioTranslator<SANDLERDIMAGGIOPARENTTRANSLATOR>::TPZSandlerDimaggioTranslator() {
}

template<typename SANDLERDIMAGGIOPARENTTRANSLATOR>
TPZSandlerDimaggioTranslator<SANDLERDIMAGGIOPARENTTRANSLATOR>::TPZSandlerDimaggioTranslator(const TPZSandlerDimaggioTranslator<SANDLERDIMAGGIOPARENTTRANSLATOR>& orig) {
}

template<typename SANDLERDIMAGGIOPARENTTRANSLATOR>
void TPZSandlerDimaggioTranslator<SANDLERDIMAGGIOPARENTTRANSLATOR>::UpdateStream(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    auto old_version = chunk.mOldVersion["NeoPZ"];
    switch (old_version) {
        case 1:
            UpdateFromV1(chunk, toVersion);
            break;
        default:
            UpdateAttributes(chunk, toVersion);
            break;
    }
}

template<typename SANDLERDIMAGGIOPARENTTRANSLATOR>
void TPZSandlerDimaggioTranslator<SANDLERDIMAGGIOPARENTTRANSLATOR>::UpdateFromV1(TPZChunkInTranslation &chunk, const std::map<std::string, uint64_t> &toVersion) {
    parentTranslator.UpdateStream(chunk, toVersion);
    REAL fA, fB, fC, fD, fR, fW;
    
    chunk.mOldStream.Read(&fA);
    chunk.mOldStream.Read(&fB);
    chunk.mOldStream.Read(&fC);
    chunk.mOldStream.Read(&fD);
    chunk.mOldStream.Read(&fR);
    chunk.mOldStream.Read(&fW);	
    
    REAL fLambda;
    chunk.mOldStream.Read(&fLambda);	
    REAL fMu;
    chunk.mOldStream.Read(&fMu);	
    
    REAL fResTol;
    chunk.mOldStream.Read(&fResTol);	
    REAL fIntegrTol;
    chunk.mOldStream.Read(&fIntegrTol);	
    int fMaxNewton;
    chunk.mOldStream.Read(&fMaxNewton);	
    REAL fMinLambda;
    chunk.mOldStream.Read(&fMinLambda);

    TPZManVector<REAL, 6> fEpsT_fData(6);
    chunk.mOldStream.Read(fEpsT_fData);
    TPZManVector<REAL, 6> fEpsP_fData(6);
    chunk.mOldStream.Read(fEpsP_fData);
    REAL m_hardening;
    chunk.mOldStream.Read(&m_hardening);
}
 
template<typename SANDLERDIMAGGIOPARENTTRANSLATOR>
void TPZSandlerDimaggioTranslator<SANDLERDIMAGGIOPARENTTRANSLATOR>::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    parentTranslator.UpdateStream(chunk, toVersion);
}

template<typename SANDLERDIMAGGIOPARENTTRANSLATOR>
TPZSandlerDimaggioTranslator<SANDLERDIMAGGIOPARENTTRANSLATOR>::~TPZSandlerDimaggioTranslator() {
}



#endif /* TPZSANDLERDIMAGGIOTRANSLATOR_H */

