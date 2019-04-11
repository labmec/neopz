//
// Created by natalia on 28/03/19.
//

#include "TRKSolution.h"

TRKSolution::TRKSolution() {

}

TRKSolution::TRKSolution(const TRKSolution &  other) {
    m_material = other.m_material;
    m_re = other.m_re;
    m_rw = other.m_rw;
    m_sigma = other.m_sigma;
    m_pw = other.m_pw;
    m_sigma0 = other.m_sigma0;
    m_theta = other.m_theta;
}

TRKSolution & TRKSolution::operator=(const TRKSolution &  other) {
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    m_material = other.m_material;
    m_re = other.m_re;
    m_rw = other.m_rw;
    m_sigma = other.m_sigma;
    m_pw = other.m_pw;
    m_sigma0 = other.m_sigma0;
    m_theta = other.m_theta;

    return *this;
}

TRKSolution::~TRKSolution() {

}

void TRKSolution::SetMaterial(TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > *material) {
    m_material = material;
}

TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > * TRKSolution::Material() {
    return m_material;
}

void TRKSolution::SetExternRadius(REAL re) {
    m_re = re;
}

REAL TRKSolution::ExternRadius() {
    return m_re;
}

void TRKSolution::SetWellboreRadius(REAL rw) {
    m_rw = rw;
}

REAL TRKSolution::WellboreRadius() {
    return m_rw;
}

void TRKSolution::SetStressXYZ(TPZTensor<REAL> &sigma, REAL theta) {
    m_sigma = sigma;
    m_theta = theta;
}

TPZTensor<REAL> TRKSolution::StressXYZ() {
    return m_sigma;
}

REAL TRKSolution::Theta() {
    return m_theta;
}

void TRKSolution::SetWellborePressure(REAL pw) {
    m_pw = pw;
}

REAL TRKSolution::WellborePressure() {
    return m_pw;
}

void TRKSolution::SetInitialStress(REAL sigma0) {
    m_sigma0 = sigma0;
}

REAL TRKSolution::InitialStress() {
    return m_sigma0;
}

void TRKSolution::F (REAL r, REAL ur, REAL sigma_r, REAL &d_ur, REAL &d_sigmar) {
    REAL lambda = m_material->GetDefaultMemory().m_ER.Lambda();
    REAL G = m_material->GetDefaultMemory().m_ER.G();

    d_ur = (r*sigma_r-lambda*ur)/(r*lambda+2*G*r);
    d_sigmar = (-sigma_r + (2*G*ur/r + lambda*(ur/r + (r*sigma_r-lambda*ur)/(r*(lambda + 2*G)))))/r;
}

void TRKSolution::ParametersAtRe(TPZFNMatrix<3,REAL> &sigma, REAL &u_re) {
    REAL nu = m_material->GetDefaultMemory().m_ER.Poisson();
    REAL G = m_material->GetDefaultMemory().m_ER.G();

    sigma.Resize(1,3);

    // Stress at re
    sigma(0,0) = (m_sigma.XX() + m_sigma.YY()) / 2 * (1 - pow(m_rw / m_re, 2)) + (1 - 4 * pow(m_rw / m_re, 2) + 3 * pow(m_rw / m_re, 4)) * (m_sigma.XX() - m_sigma.YY()) / 2 * cos(2 * m_theta) +
                        m_sigma.XY() * (1 - 4 * pow(m_rw / m_re, 2) + 3 * pow(m_rw / m_re, 4)) * sin(2 * m_theta) + m_pw * pow(m_rw / m_re, 2) - m_sigma0;
    sigma(0,1) = (m_sigma.XX() + m_sigma.YY()) / 2 * (1 + pow(m_rw / m_re, 2)) - (1 + 3 * pow(m_rw / m_re, 4)) * (m_sigma.XX() - m_sigma.YY()) / 2 * cos(2 * m_theta) -
                        m_sigma.XY() * (1 + 3 * pow(m_rw / m_re, 4)) * sin(2 * m_theta) - m_pw * pow(m_rw / m_re, 2) - m_sigma0;
    sigma(0,2) = m_sigma.ZZ() - nu * (2 * (m_sigma.XX() - m_sigma.YY()) * (m_rw / m_re) * (m_rw / m_re) * cos(2 * m_theta) + 4 * m_sigma.XY() * (m_rw / m_re) * (m_rw / m_re) * sin(2 * m_theta)) - m_sigma0;

    // Displacement at re
    u_re = -(m_re * sigma(0,2) - m_re * sigma(0,1)) / (2 * G);
}

void TRKSolution::RKProcess(int np, std::ostream &out, bool euler) {
    REAL h = (m_rw - m_re) / np;

    REAL lambda = m_material->GetDefaultMemory().m_ER.Lambda();
    REAL G = m_material->GetDefaultMemory().m_ER.G();

    TPZVec<REAL> r(np+1);
    TPZVec<REAL> u(np+1);
    TPZFNMatrix<3,REAL> sigma(np+1,3,0.);

    // Displacement and stress at re
    TPZFNMatrix<3,REAL> sigma_re;
    REAL u_re;
    ParametersAtRe(sigma_re, u_re);

    r[0] = m_re;
    u[0] = u_re;
    sigma(0,0) = sigma_re(0,0);
    sigma(0,1) = sigma_re(0,1);
    sigma(0,2) = sigma_re(0,2);

    for (int i = 0; i < np; i++) {
        REAL du_k1;
        REAL dsigma_rr_k1;

        //k1
        F(r[i], u[i], sigma(i,0), du_k1, dsigma_rr_k1);

        if (euler == false) {
            REAL du_k2, du_k3, du_k4;
            REAL dsigma_rr_k2, dsigma_rr_k3, dsigma_rr_k4;
            //k2
            F(r[i] + h / 2., u[i] + h * du_k1 / 2., sigma(i, 0) + h * dsigma_rr_k1 / 2., du_k2, dsigma_rr_k2);

            //k3
            F(r[i] + h / 2., u[i] + h * du_k2 / 2., sigma(i, 0) + h * dsigma_rr_k2 / 2., du_k3, dsigma_rr_k3);

            //k4
            F(r[i] + h, u[i] + h * du_k3, sigma(i, 0) + h * dsigma_rr_k3, du_k4, dsigma_rr_k4);

            //u_ip1, sigma_ip1
            u[i + 1] = u[i] + 1. / 6. * h * (du_k1 + 2. * du_k2 + 2. * du_k3 + du_k4);
            sigma(i+1,0) = sigma(i,0) + 1. / 6. * h * (dsigma_rr_k1 + 2. * dsigma_rr_k2 + 2. * dsigma_rr_k3 + dsigma_rr_k4);

        } else if (euler == true) {
            //u_ip1, sigma_ip1
            u[i + 1] = u[i] + h * du_k1;
            sigma(i+1,0) = sigma(i,0) + h * dsigma_rr_k1;
        }

        r[i + 1] = r[i] + h;
        sigma(i+1,1) = 2 * G * u[i + 1] / r[i + 1] + lambda * (u[i + 1] / r[i + 1] + (r[i + 1] * sigma(i+1,0) - lambda * u[i + 1]) / (r[i + 1] * (lambda + 2 * G)));
        sigma(i+1,2) = (-lambda*(lambda* (sigma(i+1,0) - 2* sigma(i+1,0) - sigma(i+1,1)) - 2 *G *(sigma(i+1,0) + sigma(i+1,1))))/(2*(lambda + G)*(lambda + 2*G));
    }

    out << "radius" << "  " << "u" << "   " << "sigma_rr" << "   " << "sigma_tt" << "   " << "sigma_zz" << std::endl;
    for (int i = 0; i < np; i++) {
        out << r[i] << "  " << u[i] << "   " << sigma(i,0) << "   " << sigma(i,1) << "   " << sigma(i,2) << std::endl;
    }
}
