//
//  TPZPorousElasticResponseTranslator.cpp
//  pz
//
//  Created by Omar Durán on 9/5/19.
//

#include "pzreal.h"
#include "TPZPorousElasticResponseTranslator.h"
#include "TPZChunkInTranslation.h"

TPZPorousElasticResponseTranslator::TPZPorousElasticResponseTranslator() {
}

TPZPorousElasticResponseTranslator::TPZPorousElasticResponseTranslator(const TPZPorousElasticResponseTranslator & other) {
    
}

void TPZPorousElasticResponseTranslator::UpdateAttributes(TPZChunkInTranslation & chunk, const std::map<std::string, uint64_t>& toVersion) {
    
    /// Logarithmic bulk modulus
    STATE m_kappa;
    chunk.mOldStream.Read(&m_kappa);
    chunk.mNewStream.Write(&m_kappa);
    
    /// Elastic tensile strengh
    STATE m_pt_el;
    chunk.mOldStream.Read(&m_pt_el);
    chunk.mNewStream.Write(&m_pt_el);
    
    /// Initial void ratio
    STATE m_e_0;
    chunk.mOldStream.Read(&m_e_0);
    chunk.mNewStream.Write(&m_e_0);
    
    /// Initial equivalent pressure stress
    STATE m_p_0;
    chunk.mOldStream.Read(&m_p_0);
    chunk.mNewStream.Write(&m_p_0);
    
    /// Poisson ratio
    STATE m_nu;
    chunk.mOldStream.Read(&m_nu);
    chunk.mNewStream.Write(&m_nu);
    
    /// Second lamé parameter
    STATE m_mu;
    chunk.mOldStream.Read(&m_mu);
    chunk.mNewStream.Write(&m_mu);
    
    /// Directive for define constant shear modulus calculations (false means constant Poisson ratio)
    bool m_is_G_constant_Q;
    chunk.mOldStream.Read(m_is_G_constant_Q);
    chunk.mNewStream.Write(m_is_G_constant_Q);
    
    /// Directive for define Plain stress state or plane strain state
    bool m_plane_stress_Q;
    chunk.mOldStream.Read(m_plane_stress_Q);
    chunk.mNewStream.Write(m_plane_stress_Q);
    
}


TPZPorousElasticResponseTranslator::~TPZPorousElasticResponseTranslator() {
    
}
