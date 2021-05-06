//
//  TPZElasticResponse.h
//  pz
//
//  Created by Erick Slis Raggio Santos on 04/07/2009.
//

#include "TPZElasticResponse.h"

int TPZElasticResponse::ClassId() const{
    return Hash("TPZElasticResponse");
}

TPZElasticResponse::TPZElasticResponse() : m_lambda(0.), m_mu(0.) {
    m_epsilon_star.Zero();
    m_sigma_star.Zero();
}

TPZElasticResponse::TPZElasticResponse(const TPZElasticResponse & other) {
    m_lambda        = other.m_lambda;
    m_mu            = other.m_mu;
    m_epsilon_star      = other.m_epsilon_star;
    m_sigma_star    = other.m_sigma_star;
}

TPZElasticResponse & TPZElasticResponse::operator=(const TPZElasticResponse & other) {
    m_lambda        = other.m_lambda;
    m_mu            = other.m_mu;
    m_epsilon_star      = other.m_epsilon_star;
    return *this;
}


void TPZElasticResponse::Write(TPZStream& buf, int withclassid) const { //ok
    buf.Write(&m_lambda);
    buf.Write(&m_mu);
    m_epsilon_star.Write(buf,withclassid);
    m_sigma_star.Write(buf,withclassid);
}

void TPZElasticResponse::Read(TPZStream& buf, void* context) { //ok
    buf.Read(&m_lambda);
    buf.Read(&m_mu);
    m_epsilon_star.Read(buf, context);
    m_sigma_star.Read(buf, context);
}


const char * TPZElasticResponse::Name() const {
    return "TPZElasticResponse";
}

void TPZElasticResponse::Print(std::ostream & out) const {
    out << this->Name();
    out << "\n Young = " << E();
    out << "\n Poisson = " << Poisson();
    out << "\n m_lambda = " << m_lambda;
    out << "\n m_mu = " << m_mu;
    m_epsilon_star.Print(out);
    m_sigma_star.Print(out);
}

void TPZElasticResponse::De(TPZFMatrix<REAL> & De) {
    REAL Mu2 = 2 * m_mu;
    
    De.Zero();
    
    De(_XX_, _XX_) += m_lambda;
    De(_XX_, _YY_) += m_lambda;
    De(_XX_, _ZZ_) += m_lambda;
    De(_YY_, _XX_) += m_lambda;
    De(_YY_, _YY_) += m_lambda;
    De(_YY_, _ZZ_) += m_lambda;
    De(_ZZ_, _XX_) += m_lambda;
    De(_ZZ_, _YY_) += m_lambda;
    De(_ZZ_, _ZZ_) += m_lambda;
    
    int i;
    for (i = 0; i < 6; i++)De(i, i) += Mu2;
}

void TPZElasticResponse::SetEngineeringData(REAL Eyoung, REAL Poisson) {
    
    m_lambda = Poisson * Eyoung / ((1. + Poisson)*(1. - 2. * Poisson));
    m_mu = Eyoung / (2. * (1. + Poisson));
}

void TPZElasticResponse::SetLameData(REAL lambda, REAL mu) {
    
    m_lambda = lambda;
    m_mu = mu;
}

REAL TPZElasticResponse::Lambda() const {
    return m_lambda;
}

REAL TPZElasticResponse::K() const {
    return m_lambda + 2. * m_mu / 3.;
}

REAL TPZElasticResponse::Mu() const {
    return m_mu;
}

REAL TPZElasticResponse::G() const {
    return Mu();
}

REAL TPZElasticResponse::E() const {
    REAL E = m_mu * (3. * m_lambda + 2. * m_mu) / (m_lambda + m_mu);
    return E;
}

REAL TPZElasticResponse::Poisson() const {
    REAL poisson = m_lambda / (2. * (m_lambda + m_mu));
    return poisson;
}

void TPZElasticResponse::SetReferenceStrainData(TPZTensor<REAL> & eps_star){
    m_epsilon_star = eps_star;
}

TPZTensor<REAL> & TPZElasticResponse::ReferenceStrainData(){
    return m_epsilon_star;
}

void TPZElasticResponse::SetReferenceStressData(TPZTensor<REAL> & sigma_star){
    m_sigma_star = sigma_star;
}

TPZTensor<REAL> & TPZElasticResponse::ReferenceStressData(){
    return m_sigma_star;
}

