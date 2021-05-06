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
    tpzTensorTranslatorT.UpdateStream(chunk, toVersion); // m_eps_t
    tpzTensorTranslatorT.UpdateStream(chunk, toVersion); // m_eps_p
    T m_hardening;
    chunk.mOldStream.Read(&m_hardening);
    chunk.mNewStream.Write(&m_hardening);
    int m_m_type = 1;
    chunk.mNewStream.Write(&m_m_type);
}

template <typename T>
void TPZPlasticStateTranslator<T>::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    tpzTensorTranslatorT.UpdateStream(chunk, toVersion); // m_eps_t
    tpzTensorTranslatorT.UpdateStream(chunk, toVersion); // m_eps_p
    T m_hardening;
    chunk.mOldStream.Read(&m_hardening);
    chunk.mNewStream.Write(&m_hardening);
    int m_m_type;
    chunk.mOldStream.Read(&m_m_type);
    chunk.mNewStream.Write(&m_m_type);
}

template <typename T>
TPZPlasticStateTranslator<T>::~TPZPlasticStateTranslator() {
}

#endif /* TPZPLASTICSTATETRANSLATOR_H */

