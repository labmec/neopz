/* 
 * File:   TPZElastoPlasticMemTranslator.cpp
 * Author: thiago
 * 
 * Created on 12 de Março de 2018, 18:58
 */

#include "TPZElastoPlasticMemTranslator.h"
#include "TPZTensor.h"

TPZElastoPlasticMemTranslator::TPZElastoPlasticMemTranslator() {
}

TPZElastoPlasticMemTranslator::TPZElastoPlasticMemTranslator(const TPZElastoPlasticMemTranslator& orig) {
}

void TPZElastoPlasticMemTranslator::UpdateStream(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
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

void TPZElastoPlasticMemTranslator::UpdateFromV1(TPZChunkInTranslation &chunk, const std::map<std::string, uint64_t> &toVersion) {
    TPZManVector<REAL, 6> fSigma(6);
    chunk.mOldStream.Read(&fSigma[0],6);
    chunk.mNewStream.Write(fSigma);
    TPZManVector<REAL, 6> m_eps_t(6);
    chunk.mOldStream.Read(&m_eps_t[0],6);
    chunk.mNewStream.Write(m_eps_t);
    TPZManVector<REAL, 6> m_eps_p(6);
    chunk.mOldStream.Read(&m_eps_p[0],6);
    chunk.mNewStream.Write(m_eps_p);
    REAL m_hardening;
    chunk.mOldStream.Read(&m_hardening);
    chunk.mNewStream.Write(&m_hardening);
    int m_m_type = 1;
    chunk.mNewStream.Write(&m_m_type);
    int fPlasticSteps;
    chunk.mOldStream.Read(&fPlasticSteps,1);
    chunk.mNewStream.Write(&fPlasticSteps,1);
    TPZManVector<REAL,3> fDisplacement(3);
    chunk.mOldStream.Read(&fDisplacement[0],3);
    chunk.mNewStream.Write(fDisplacement);
}


void TPZElastoPlasticMemTranslator::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    tpzTensorTranslatorREAL.UpdateStream(chunk, toVersion);  // fSigma
    tpzPlasticStateTranslatorREAL.UpdateStream(chunk, toVersion); // fPlasticState
    int fPlasticSteps;
    chunk.mOldStream.Read(&fPlasticSteps,1);
    chunk.mNewStream.Write(&fPlasticSteps,1);
    TPZManVector<REAL,3> fDisplacement;
    chunk.mOldStream.Read(fDisplacement);
    chunk.mNewStream.Write(fDisplacement);
}

TPZElastoPlasticMemTranslator::~TPZElastoPlasticMemTranslator() {
}

int TPZElastoPlasticMemTranslator::GetClassId() const {
    return classid;
}

void TPZElastoPlasticMemTranslator::SetClassId(int classid) {
    TPZElastoPlasticMemTranslator::classid = classid;
}

int TPZElastoPlasticMemTranslator::classid = -1;
