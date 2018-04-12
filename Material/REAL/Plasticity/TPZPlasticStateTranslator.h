/* 
 * File:   TPZPlasticStateTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 14:54
 */

#ifndef TPZPLASTICSTATETRANSLATOR_H
#define TPZPLASTICSTATETRANSLATOR_H

#include "TPZChunkTranslator.h"
#include "TPZTensorTranslator.h"

template <typename T>
class TPZPlasticStateTranslator : public TPZChunkTranslator {
public:
    TPZPlasticStateTranslator();
    TPZPlasticStateTranslator(const TPZPlasticStateTranslator<T>& orig);
    
    virtual void UpdateStream(TPZChunkInTranslation &chunk, const std::map<std::string, uint64_t> &toVersion);
    void UpdateAttributes(TPZChunkInTranslation &chunk, const std::map<std::string, uint64_t> &toVersion);
    
    virtual ~TPZPlasticStateTranslator();
private:
    void UpdateFromV1(TPZChunkInTranslation &chunk, const std::map<std::string, uint64_t> &toVersion);
    TPZTensorTranslator<T> tpzTensorTranslatorT;
};


template <typename T>
TPZPlasticStateTranslator<T>::TPZPlasticStateTranslator() {
}

template <typename T>
TPZPlasticStateTranslator<T>::TPZPlasticStateTranslator(const TPZPlasticStateTranslator<T>& orig) {
}

template <typename T>
void TPZPlasticStateTranslator<T>::UpdateStream(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
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

template <typename T>
void TPZPlasticStateTranslator<T>::UpdateFromV1(TPZChunkInTranslation &chunk, const std::map<std::string, uint64_t> &toVersion) {
    tpzTensorTranslatorT.UpdateStream(chunk, toVersion); // fEpsT
    tpzTensorTranslatorT.UpdateStream(chunk, toVersion); // fEpsP
    T fAlpha;
    chunk.mOldStream.Read(&fAlpha);
    chunk.mNewStream.Write(&fAlpha);
    int fMType = 1;
    chunk.mNewStream.Write(&fMType);
}

template <typename T>
void TPZPlasticStateTranslator<T>::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    tpzTensorTranslatorT.UpdateStream(chunk, toVersion); // fEpsT
    tpzTensorTranslatorT.UpdateStream(chunk, toVersion); // fEpsP
    T fAlpha;
    chunk.mOldStream.Read(&fAlpha);
    chunk.mNewStream.Write(&fAlpha);
    int fMType;
    chunk.mOldStream.Read(&fMType);
    chunk.mNewStream.Write(&fMType);
}

template <typename T>
TPZPlasticStateTranslator<T>::~TPZPlasticStateTranslator() {
}

#endif /* TPZPLASTICSTATETRANSLATOR_H */

