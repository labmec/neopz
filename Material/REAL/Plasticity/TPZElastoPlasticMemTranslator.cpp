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
    TPZManVector<REAL, 6> m_sigma(6);
    chunk.mOldStream.Read(&m_sigma[0],6);
    chunk.mNewStream.Write(m_sigma);
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
    int m_plastic_steps;
    chunk.mOldStream.Read(&m_plastic_steps,1);
    chunk.mNewStream.Write(&m_plastic_steps,1);
    TPZManVector<REAL,3> m_u(3);
    chunk.mOldStream.Read(&m_u[0],3);
    chunk.mNewStream.Write(m_u);
}


void TPZElastoPlasticMemTranslator::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    tpzTensorTranslatorREAL.UpdateStream(chunk, toVersion);  // m_sigma
    tpzPlasticStateTranslatorREAL.UpdateStream(chunk, toVersion); // m_elastoplastic_state
    int m_plastic_steps;
    chunk.mOldStream.Read(&m_plastic_steps,1);
    chunk.mNewStream.Write(&m_plastic_steps,1);
    TPZManVector<REAL,3> m_u;
    chunk.mOldStream.Read(m_u);
    chunk.mNewStream.Write(m_u);
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
