//
// Created by natalia on 28/03/19.
//

#include "TRKSolution.h"

TRKSolution::TRKSolution() {

}

TRKSolution::TRKSolution(const TRKSolution &  other) : m_initial_state_memory() {
    m_elastoplastic_model = other.m_elastoplastic_model;
    m_re = other.m_re;
    m_rw = other.m_rw;
    m_ure = other.m_ure;
    m_sigma_re = other.m_sigma_re;
    m_n_points = other.m_n_points;
    m_memory_vector.resize(0);
    m_initial_state_memory = other.m_initial_state_memory;

}

TRKSolution & TRKSolution::operator=(const TRKSolution &  other) {
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    m_elastoplastic_model = other.m_elastoplastic_model;
    m_re = other.m_re;
    m_rw = other.m_rw;
    m_ure = other.m_ure;
    m_sigma_re = other.m_sigma_re;
    m_n_points = other.m_n_points;
    m_memory_vector = other.m_memory_vector;
    m_initial_state_memory = other.m_initial_state_memory;
    return *this;
}

TRKSolution::~TRKSolution() {

}

void TRKSolution::SetElastoPlasticModel(TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> & model){
    m_elastoplastic_model = model;
}

void TRKSolution::SetExternalRadius(REAL re) {
    m_re = re;
}

REAL TRKSolution::ExternalRadius() {
    return m_re;
}

void TRKSolution::SetWellboreRadius(REAL rw) {
    m_rw = rw;
}

REAL TRKSolution::WellboreRadius() {
    return m_rw;
}


void TRKSolution::SetRadialDisplacement(REAL ure) {
    m_ure = ure;
}

REAL TRKSolution::RadialRadialDisplacement() {
    return m_ure;
}

void TRKSolution::SetRadialStress(REAL sigma_re) {
    m_sigma_re = sigma_re;
}

REAL TRKSolution::RadialStress() {
    return m_sigma_re;
}

void TRKSolution::SetInitialStateMemory(TPZElastoPlasticMem memory){
    m_initial_state_memory = memory;
}

void TRKSolution::SetNumberOfPoints(int n_points) {
    m_n_points = n_points;
}

int TRKSolution::GetNumberOfPoints() {
    return m_n_points;
}

void TRKSolution::FillPointsMemory(){
    m_memory_vector.Resize(m_n_points+1);
    for (auto & item : m_memory_vector) {
        item = m_initial_state_memory;
    }
}

void TRKSolution::RKProcess(std::ostream &out, bool euler) {
    
    if (m_n_points <= 0) {
        DebugStop();
    }
    
    REAL h = (m_rw - m_re) / m_n_points;
    TPZVec<REAL> r(m_n_points+1);
    TPZVec<REAL> u(m_n_points+1);
    TPZFNMatrix<10,REAL> sigma(m_n_points+1,3,0.);
    
    // Displacement and stress at re
    r[0] = m_re;
    u[0] = m_ure;
    
    /// NVB the RK approximation is based on Hydristatic assumption.
    sigma(0,0) = m_sigma_re;
    sigma(0,1) = m_sigma_re;
    sigma(0,2) = m_sigma_re;
    
    m_memory_vector[0].m_sigma.XX() = m_sigma_re;
    m_memory_vector[0].m_sigma.YY() = m_sigma_re;
    m_memory_vector[0].m_sigma.ZZ() = m_sigma_re;
    
    REAL lambda = m_memory_vector[0].m_ER.Lambda();
    REAL G = m_memory_vector[0].m_ER.G();

    for (int i = 0; i < m_n_points; i++) {
        REAL du_k1;
        REAL dsigma_rr_k1;
        
        /// Assuming that Lamé parameters suffer small change between two RK points
        /// http://www.ecs.umass.edu/~arwade/courses/str-mech/polar.pdf
        
        //k1
        F(r[i], u[i], sigma(i,0), du_k1, dsigma_rr_k1, lambda, G);
        
        if (euler == false) {
            REAL du_k2, du_k3, du_k4;
            REAL dsigma_rr_k2, dsigma_rr_k3, dsigma_rr_k4;
            //k2
            F(r[i] + h / 2., u[i] + h * du_k1 / 2., sigma(i, 0) + h * dsigma_rr_k1 / 2., du_k2, dsigma_rr_k2, lambda, G);

            //k3
            F(r[i] + h / 2., u[i] + h * du_k2 / 2., sigma(i, 0) + h * dsigma_rr_k2 / 2., du_k3, dsigma_rr_k3, lambda, G);

            //k4
            F(r[i] + h, u[i] + h * du_k3, sigma(i, 0) + h * dsigma_rr_k3, du_k4, dsigma_rr_k4, lambda, G);

            //u_ip1, sigma_ip1
            r[i + 1] = r[i] + h;
            u[i + 1] = u[i] + 1. / 6. * h * (du_k1 + 2. * du_k2 + 2. * du_k3 + du_k4);
            sigma(i+1,0) = sigma(i,0) + 1. / 6. * h * (dsigma_rr_k1 + 2. * dsigma_rr_k2 + 2. * dsigma_rr_k3 + dsigma_rr_k4);

        } else if (euler == true) {
        
            //u_ip1, sigma_ip1
            r[i + 1] = r[i] + h;
            u[i + 1] = u[i] + h * du_k1;
            sigma(i+1,0) = sigma(i,0) + h * dsigma_rr_k1;
        }
        
        /// update elastoplastic state
        REAL sigma_r, sigma_t, sigma_z;
        {
            REAL r_pone = r[i + 1];
            REAL ur_pone = u[i + 1];
            REAL last_sigma_r = sigma(i,0);
            sigma_r = last_sigma_r + h * dsigma_rr_k1;
            sigma_t = (lambda*r_pone*sigma_r + 4*G*(G + lambda)*ur_pone)/((2*G + lambda)*r_pone);
            REAL nu = lambda / (2.0*(lambda+G));
            REAL Ey = G * (3.0*lambda+2.0*G) / (lambda+G);
            REAL eps_r = (1+nu)*((1-nu)*sigma_r - nu*sigma_t) / Ey;
            REAL eps_t = (1+nu)*((1-nu)*sigma_t - nu*sigma_r) / Ey;
            sigma_z = lambda*(eps_r + ur_pone/r_pone);
            
            TPZTensor<REAL> eps_total, eps_e, eps_p ,sigma;
            eps_total.Zero();
            eps_total.XX() = eps_r;
            eps_total.YY() = eps_t;
        
            sigma.Zero();
            sigma.XX() = sigma_r;
            sigma.YY() = sigma_t;
            sigma.ZZ() = sigma_z;
            
            TPZPlasticState<REAL> state = m_memory_vector[i+1].m_elastoplastic_state;
            TPZFMatrix<REAL> Dep(6,6,0.0);
            m_elastoplastic_model.SetState(state);
            m_elastoplastic_model.ApplyStrainComputeSigma(eps_total, sigma, &Dep);
            lambda = Dep(0,5);
            G = Dep(4,4)/2.0;
            
            state = m_elastoplastic_model.GetState();
            m_memory_vector[i+1].m_elastoplastic_state = state;
            m_memory_vector[i+1].m_sigma = sigma;

//            std::cout << "lambda = " << lambda << std::endl;
//            std::cout << "G = " << G << std::endl;
            
            
        }
        
        sigma(i+1,0) = sigma_r;
        sigma(i+1,1) = sigma_t;
        sigma(i+1,2) = sigma_z;
        
    }
    
    out << "radius" << "  " << "u" << "   " << "sigma_rr" << "   " << "sigma_tt" << "   " << "sigma_zz" << "  " << "eps_rr" << "   " << "eps_tt" << "   " << "eps_zz" << "  " << "eps_p_rr" << "   " << "eps_p_tt" << "   " << "eps_p_zz" << std::endl;
    for (int i = 0; i < m_n_points; i++) {
        TPZTensor<REAL> eps_t = m_memory_vector[i].m_elastoplastic_state.m_eps_t;
        TPZTensor<REAL> eps_p = m_memory_vector[i].m_elastoplastic_state.m_eps_p;
        TPZTensor<REAL> sigma = m_memory_vector[i].m_sigma;
        
        out << r[i] << "  " << u[i] << "   " << sigma.XX() << "   " << sigma.YY() << "   " << sigma.ZZ() << "   " << eps_t.XX() << "   " << eps_t.YY() << "   " << eps_t.ZZ() << "   " << eps_p.XX() << "   " << eps_p.YY() << "   " << eps_p.ZZ() << std::endl;
    }
}

void TRKSolution::F (REAL r, REAL ur, REAL sigma_r, REAL &d_ur, REAL &d_sigmar, REAL & lambda, REAL & G){
    
    REAL sigma_t = (lambda*r*sigma_r + 4*G*(G + lambda)*ur)/((2*G + lambda)*r);
    d_ur = (r*sigma_r-lambda*ur)/(r*lambda+2*G*r);
    d_sigmar = (-sigma_r + sigma_t)/r;
    
}
