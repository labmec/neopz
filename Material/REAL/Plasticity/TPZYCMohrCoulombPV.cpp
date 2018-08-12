#include "TPZYCMohrCoulombPV.h"

typedef TFad<3, REAL> fadtype;

TPZYCMohrCoulombPV::TPZYCMohrCoulombPV() : fPhi(0.), fPsi(0.), fc(0.), fER(), fEpsPlasticBar(0.) {

}

TPZYCMohrCoulombPV::TPZYCMohrCoulombPV(REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER) : fPhi(Phi), fPsi(Psi), fc(c), fER(ER), fEpsPlasticBar(0.) {

}

TPZYCMohrCoulombPV::TPZYCMohrCoulombPV(const TPZYCMohrCoulombPV &cp)
{
    TPZYCMohrCoulombPV::operator=(cp);
}

TPZYCMohrCoulombPV & TPZYCMohrCoulombPV::operator=(const TPZYCMohrCoulombPV &cp) {
    fPhi = cp.fPhi;
    fPsi = cp.fPsi;
    fc = cp.fc;
    fEpsPlasticBar = cp.fEpsPlasticBar;
    fER = cp.fER;
    return *this;
}

void TPZYCMohrCoulombPV::Print(std::ostream& out) const {
    out << "TPZYCMohrCoulombPV\n";
    out << "Phi: " << fPhi << std::endl;
    out << "Psi: " << fPsi << std::endl;
    out << "c: " << fc << std::endl;
}


int TPZYCMohrCoulombPV::ClassId() const{
    return Hash("TPZYCMohrCoulombPV");
}

void TPZYCMohrCoulombPV::Read(TPZStream& buf, void* context) { //ok
    buf.Read(&fPhi);
    buf.Read(&fPsi);
    buf.Read(&fc);
    buf.Read(&fEpsPlasticBar);
    fER.Read(buf, context);
}

void TPZYCMohrCoulombPV::Write(TPZStream& buf, int withclassid) const { //ok
    buf.Write(&fPhi);
    buf.Write(&fPsi);
    buf.Write(&fc);
    buf.Write(&fEpsPlasticBar);
    fER.Write(buf, withclassid);
}

REAL TPZYCMohrCoulombPV::InitialDamage(const TPZVec<REAL> &stress_p) const{
    
    std::cout << "TPZYCMohrCoulombPV::There is no damage variable for this model at the current time." << std::endl;
    
    TPZVec<REAL> phi(3);
    REAL alpha = 0.0;
    REAL tol = 1.0e-10;
    Phi(stress_p, alpha, phi);
    bool Is_valid_stress_on_cap_Q =  fabs(phi[0]) < tol || phi[0] < 0.0;
    
    if (!Is_valid_stress_on_cap_Q) {
        std::cerr << "TPZYCMohrCoulombPV::Invalid stress state." << std::endl;
        DebugStop();
    }

    return 0.0;
}

void TPZYCMohrCoulombPV::Phi(TPZVec<STATE> sig_vec, STATE alpha, TPZVec<STATE> &phi)const {
    phi.resize(3);
    for (int i = 0; i < 3; i++) phi[i] = 0;
    phi[0] = PhiPlane(sig_vec);
    phi[2] = PhiPlane(sig_vec); // Consistency with two surfaces models
}

template <class T>
void TPZYCMohrCoulombPV::PlasticityFunction(const T epsp, T &c, T &H) const {
    c = fc; // c(epsp)
    H = 0.; // dc(epsp)/depsp
}

template<class T>
TPZVec<T> TPZYCMohrCoulombPV::SigmaElastPV(const TPZVec<T> &deform) const {
    T trace = deform[0] + deform[1] + deform[2];
    TPZVec<T> sigma(3, 0.);
    
    sigma = trace * fER.Lambda() + 2 * fER.G() * deform[0];
    sigma = trace * fER.Lambda() + 2 * fER.G() * deform[1];
    sigma = trace * fER.Lambda() + 2 * fER.G() * deform[2];
    
    return sigma;
}

template<class T>
T TPZYCMohrCoulombPV::PhiPlane(const TPZVec<T> &sigma) const {
    const REAL sinphi = sin(fPhi);
    const REAL cosphi = cos(fPhi);
    
    T c, H;
    PlasticityFunction(T(fEpsPlasticBar), c, H);

    return sigma[0] - sigma[2] + (sigma[0] + sigma[2]) * sinphi - 2. * c*cosphi;
}

template<class T>
bool TPZYCMohrCoulombPV::ReturnMapPlane(const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew) const {
    sigma_projected = sigma_trial;
    TPZManVector<T, 3> eigenvalues = sigma_projected;
    const REAL sinphi = sin(fPhi);
    const REAL sinpsi = sin(fPsi);
    const REAL cosphi = cos(fPhi);
    const REAL sinphi2 = sinphi*sinphi;
    const REAL cosphi2 = 1. - sinphi2;
    const REAL constA = 4. * fER.G() *(1. + sinphi * sinpsi / 3.) + 4. * fER.K() * sinphi*sinpsi;
    T c, H;
    T epsbar = T(fEpsPlasticBar + memory.fGamma[0]*2. * cosphi);
    PlasticityFunction(epsbar, c, H);
    T phi = eigenvalues[0] - eigenvalues[2]+(eigenvalues[0] + eigenvalues[2]) * sinphi - 2. * c*cosphi;
    T gamma = memory.fGamma[0];
    REAL phival = TPZExtractVal::val(phi);
    int n_iterations = 30; // @TODO : Define a numeric controls manager object and use it to obtain this information
    int i;
    bool stop_criterion;
    for (i = 0; i < n_iterations; i++) {
        T jac = -constA - T(4. * cosphi2) * H;
        T delta_gamma = - phi / jac;
        gamma += delta_gamma;
        phi = eigenvalues[0] - eigenvalues[2]+(eigenvalues[0] + eigenvalues[2]) * sinphi - 2. * c * cosphi - constA * gamma;
        phival = TPZExtractVal::val(phi);
        stop_criterion = IsZero(phival);
        if (stop_criterion) {
            break;
        }
    }

#ifdef PZDEBUG
    if (i == n_iterations) {
        DebugStop();
    }
#endif
    
    epsbar = T(fEpsPlasticBar) + gamma * T(2. * cosphi);
    memory.fGamma[0] = TPZExtractVal::val(gamma);
    eigenvalues[0] -= T(2. * fER.G()*(1 + sinpsi / 3.) + 2. * fER.K() * sinpsi) * gamma;
    eigenvalues[1] += T((4. * fER.G() / 3. - fER.K()*2.) * sinpsi) * gamma;
    eigenvalues[2] += T(2. * fER.G()*(1 - sinpsi / 3.) - 2. * fER.K() * sinpsi) * gamma;
    sigma_projected = eigenvalues;
    epsbarnew = TPZExtractVal::val(epsbar);
    
    bool check_validity_Q = (TPZExtractVal::val(eigenvalues[0]) > TPZExtractVal::val(eigenvalues[1]) || IsZero(eigenvalues[0]-eigenvalues[1])) && (TPZExtractVal::val(eigenvalues[1]) > TPZExtractVal::val(eigenvalues[2]) || IsZero(eigenvalues[1]-eigenvalues[2]));
    return (check_validity_Q);
}

void TPZYCMohrCoulombPV::ComputePlaneTangent(TPZMatrix<REAL> &tang, REAL &epsbarp) const {
    
    const REAL sin_phi = sin(fPhi);
    const REAL sin_psi = sin(fPsi);
    const REAL G = fER.G(), K = fER.K();
    const REAL denominator = 6.0 * G + 2.0 * (G + 3.0 * K) * sin_phi * sin_psi;
    
    REAL epsbar = epsbarp;
    REAL c, H;
    PlasticityFunction(epsbar, c, H);
    
    tang.Redim(3, 3);
    
    // First column
    tang(0, 0) = (sin_phi - 1.0) * (-3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;
    tang(1, 0) = (2.0 * G - 3.0 * K) * (sin_phi + 1.0) * sin_psi / denominator;
    tang(2, 0) = -(sin_phi + 1.0) * (-3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;
    
    // Second column
    tang(0, 1) = 0.0;
    tang(1, 1) = 1.0;
    tang(2, 1) = 0.0;
    
    // Third column
    tang(0, 2) = -(sin_phi - 1.0) * (3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;
    tang(1, 2) = (2.0 * G - 3.0 * K) * (sin_phi - 1.0) * sin_psi / denominator;
    tang(2, 2) = (sin_phi + 1.0) * (3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;
    
}

template<class T>
bool TPZYCMohrCoulombPV::ReturnMapLeftEdge(const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew) const {
    
    sigma_projected = sigma_trial;
    TPZManVector<T, 3> eigenvalues = sigma_projected;
    const REAL sinphi = sin(fPhi);
    const REAL sinpsi = sin(fPsi);
    const REAL cosphi = cos(fPhi);
    const REAL sinphi2 = sinphi*sinphi;
    const REAL cosphi2 = 1. - sinphi2;
    TPZManVector<T, 2> gamma(2, 0.), phi(2, 0.), sigma_bar(2, 0.), ab(2, 0.);
    gamma[0] = memory.fGamma[0];
    gamma[1] = memory.fGamma[1];
    
    TPZManVector<REAL, 2> phival(2, 0.);
    TPZManVector<TPZManVector<T, 2>, 2> jac(2), jac_inv(2);
    for (int i = 0; i < 2; i++) {
        jac[i].Resize(2, 0.);
        jac_inv[i].Resize(2, 0.);
    }
    
    sigma_bar[0] = eigenvalues[0] - eigenvalues[2]+(eigenvalues[0] + eigenvalues[2]) * T(sinphi);
    sigma_bar[1] = eigenvalues[1] - eigenvalues[2]+(eigenvalues[1] + eigenvalues[2]) * T(sinphi);
    T c, H;
    ab[0] = T(4. * fER.G()*(1 + sinphi * sinpsi / 3.) + 4. * fER.K() * sinphi * sinpsi);
    ab[1] = T(2. * fER.G()*(1. - sinphi - sinpsi - sinphi * sinpsi / 3.) + 4. * fER.K() * sinphi * sinpsi);
    T epsbar = T(fEpsPlasticBar) + (gamma[0] + gamma[1]) * T(2. * cosphi);
    PlasticityFunction(epsbar, c, H);
    
    phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - T(2. * cosphi) * c;
    phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - T(2. * cosphi) * c;
    T res = (fabs(phival[0]) + fabs(phival[1]));
    int n_iterations = 30; // @TODO : Define a numeric controls manager object and use it to obtain this information
    int i;
    bool stop_criterion;
    for (i = 0; i < n_iterations; i++) {
        
        jac[0][0] = -ab[0] - T(4. * cosphi2) * H;
        jac[1][0] = -ab[1] - T(4. * cosphi2) * H;
        jac[0][1] = -ab[1] - T(4. * cosphi2) * H;
        jac[1][1] = -ab[0] - T(4. * cosphi2) * H;
        
        T det_jac = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        
#ifdef PZDEBUG
        if(IsZero(det_jac)){
            std::cerr << "TPZYCMohrCoulombPV:: Singular jacobian." << std::endl;
            DebugStop();
        }
#endif

        jac_inv[0][0] = jac[1][1] / det_jac;
        jac_inv[1][0] = -jac[1][0] / det_jac;
        jac_inv[0][1] = -jac[0][1] / det_jac;
        jac_inv[1][1] = jac[0][0] / det_jac;
        
        gamma[0] -= (jac_inv[0][0] * phi[0] + jac_inv[0][1] * phi[1]);
        gamma[1] -= (jac_inv[1][0] * phi[0] + jac_inv[1][1] * phi[1]);
        
        epsbar = T(fEpsPlasticBar)+(gamma[0] + gamma[1]) * T(2. * cosphi);
        PlasticityFunction(epsbar, c, H);
        
        phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - T(2. * cosphi) * c;
        phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - T(2. * cosphi) * c;
        phival[0] = TPZExtractVal::val(phi[0]);
        phival[1] = TPZExtractVal::val(phi[1]);
        res = (fabs(phival[0]) + fabs(phival[1]));
        
        stop_criterion = IsZero(res);
        if (stop_criterion) {
            break;
        }
    }
    
#ifdef PZDEBUG
    if (i == n_iterations) {
        DebugStop();
    }
#endif

    memory.fGamma[0] = TPZExtractVal::val(gamma[0]);
    memory.fGamma[1] = TPZExtractVal::val(gamma[1]);
    eigenvalues[0] += -T(2. * fER.G()*(1 + sinpsi / 3.) + 2. * fER.K() * sinpsi) * gamma[0] + T((4. * fER.G() / 3. - 2. * fER.K()) * sinpsi) * gamma[1];
    eigenvalues[1] += T((4. * fER.G() / 3. - fER.K()*2.) * sinpsi) * gamma[0] - T(2. * fER.G()*(1. + sinpsi / 3.) + 2. * fER.K() * sinpsi) * gamma[1];
    eigenvalues[2] += T(2. * fER.G()*(1 - sinpsi / 3.) - 2. * fER.K() * sinpsi)*(gamma[0] + gamma[1]);
    sigma_projected = eigenvalues;
    epsbarnew = TPZExtractVal::val(epsbar);

    bool check_validity_Q = (TPZExtractVal::val(eigenvalues[0]) > TPZExtractVal::val(eigenvalues[1]) || IsZero(eigenvalues[0]-eigenvalues[1])) && (TPZExtractVal::val(eigenvalues[1]) > TPZExtractVal::val(eigenvalues[2]) || IsZero(eigenvalues[1]-eigenvalues[2]));
    return (check_validity_Q);
}

/**
 * @brief Computes dsigmapr/dsigmatr for the ReturnMapLeftEdge
 */
void TPZYCMohrCoulombPV::ComputeLeftEdgeTangent(TPZMatrix<REAL> &tang, REAL &epsbarp) const {
    
    const REAL sin_phi = sin(fPhi);
    const REAL sin_psi = sin(fPsi);
    const REAL G = fER.G(), K = fER.K();
    const REAL a = 4.0 * G * (1.0 + (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;
    const REAL b = 2.0 * G * (1.0 - sin_phi - sin_psi - (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;
    
    REAL epsbar = epsbarp;
    REAL c, H;
    PlasticityFunction(epsbar, c, H);
    
    tang.Redim(3, 3);
    
    // First column
    tang(0, 0) = (-3*b*b + 3*a*(a - 2*G*(1 + sin_phi)) -
                  2*(a*G + 2*b*G + 3*a*K - 3*b*K)*(1 + sin_phi)*sin_psi)/(3.*(a - b)*(a + b));
    tang(1, 0) = (2*(1 + sin_phi)*(a*(2*G - 3*K)*sin_psi + 3*b*K*sin_psi + b*G*(3 + sin_psi)))/
    (3.*(a*a - b*b));
    tang(2, 0) = (-2*(1 + sin_phi)*(G*(-3 + sin_psi) + 3*K*sin_psi))/(3.*(a + b));
    
    // Second column
    tang(0, 1) = (2*(1 + sin_phi)*(a*(2*G - 3*K)*sin_psi + 3*b*K*sin_psi + b*G*(3 + sin_psi)))/
    (3.*(a*a - b*b));
    tang(1, 1) = (-3*b*b + 3*a*(a - 2*G*(1 + sin_phi)) -
                  2*(a*G + 2*b*G + 3*a*K - 3*b*K)*(1 + sin_phi)*sin_psi)/(3.*(a - b)*(a + b));
    tang(2, 1) = (-2*(1 + sin_phi)*(G*(-3 + sin_psi) + 3*K*sin_psi))/(3.*(a + b));
    
    // Third column
    tang(0, 2) = (2*(-1 + sin_phi)*(G*(-3 + sin_psi) - 6*K*sin_psi))/(3.*(a + b));
    tang(1, 2) = (2*(-1 + sin_phi)*(G*(-3 + sin_psi) - 6*K*sin_psi))/(3.*(a + b));
    tang(2, 2) = (3*a + 3*b - 4*(-1 + sin_phi)*(G*(-3 + sin_psi) + 3*K*sin_psi))/(3.*(a + b));
    
}

/**
 * @brief Implements the return map in the right edge of the surface
 */
template<class T>
bool TPZYCMohrCoulombPV::ReturnMapRightEdge(const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew) const {
    
    sigma_projected = sigma_trial;
    TPZManVector<T, 3> eigenvalues = sigma_projected;
    const REAL sinphi = sin(fPhi);
    const REAL sinpsi = sin(fPsi);
    const REAL cosphi = cos(fPhi);
    const REAL sinphi2 = sinphi*sinphi;
    const REAL cosphi2 = 1. - sinphi2;
    const REAL KV = fER.K();
    const REAL GV = fER.G();
    
    TPZManVector<T, 2> gamma(2, 0.), phi(2, 0.), sigma_bar(2, 0.), ab(2, 0.);
    gamma[0] = memory.fGamma[0];
    gamma[1] = memory.fGamma[1];
    TPZManVector<REAL, 2> phival(2, 0.);
    TPZManVector<TPZManVector<T, 2>, 2> jac(2), jac_inv(2);
    for (int i = 0; i < 2; i++) {
        jac[i].Resize(2, 0.);
        jac_inv[i].Resize(2, 0.);
    }
    
    sigma_bar[0] = eigenvalues[0] - eigenvalues[2]+(eigenvalues[0] + eigenvalues[2]) * T(sinphi);
    sigma_bar[1] = eigenvalues[0] - eigenvalues[1]+(eigenvalues[0] + eigenvalues[1]) * T(sinphi);
    T c, H;
    ab[0] = T(4. * GV * (1 + sinphi * sinpsi / 3.) + 4. * KV * sinphi * sinpsi);
    ab[1] = T(2. * GV * (1. + sinphi + sinpsi - sinphi * sinpsi / 3.) + 4. * KV * sinphi * sinpsi);
    T epsbar = T(fEpsPlasticBar)+(gamma[0] + gamma[1]) * T(2. * cosphi);
    PlasticityFunction(epsbar, c, H);
    
    phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - T(2. * cosphi) * c;
    phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - T(2. * cosphi) * c;

    T res = (fabs(phival[0]) + fabs(phival[1]));
    int n_iterations = 30; // @TODO : Define a numeric controls manager object and use it to obtain this information
    int i;
    bool stop_criterion;
    for (i = 0; i < n_iterations; i++) {

        jac[0][0] = -ab[0] - T(4. * cosphi2) * H;
        jac[1][0] = -ab[1] - T(4. * cosphi2) * H;
        jac[0][1] = -ab[1] - T(4. * cosphi2) * H;
        jac[1][1] = -ab[0] - T(4. * cosphi2) * H;
        
        T det_jac = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        
#ifdef PZDEBUG
        if(IsZero(det_jac)){
            std::cerr << "TPZYCMohrCoulombPV:: Singular jacobian." << std::endl;
            DebugStop();
        }
#endif
        
        jac_inv[0][0] = jac[1][1] / det_jac;
        jac_inv[1][0] = -jac[1][0] / det_jac;
        jac_inv[0][1] = -jac[0][1] / det_jac;
        jac_inv[1][1] = jac[0][0] / det_jac;
        
        gamma[0] -= (jac_inv[0][0] * phi[0] + jac_inv[0][1] * phi[1]);
        gamma[1] -= (jac_inv[1][0] * phi[0] + jac_inv[1][1] * phi[1]);
        
        epsbar = T(fEpsPlasticBar)+(gamma[0] + gamma[1]) * T(2. * cosphi);
        PlasticityFunction(epsbar, c, H);

        phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - T(2. * cosphi) * c;
        phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - T(2. * cosphi) * c;
        phival[0] = TPZExtractVal::val(phi[0]);
        phival[1] = TPZExtractVal::val(phi[1]);
        res = (fabs(phival[0]) + fabs(phival[1]));
        
        stop_criterion = IsZero(res);
        if (stop_criterion) {
            break;
        }
    }
    
#ifdef PZDEBUG
    if (i == n_iterations) {
        DebugStop();
    }
#endif

    memory.fGamma[0] = TPZExtractVal::val(gamma[0]);
    memory.fGamma[1] = TPZExtractVal::val(gamma[1]);

    eigenvalues[0] -= T(2. * GV * (1 + sinpsi / 3.) + 2. * KV * sinpsi)*(gamma[0] + gamma[1]);
    eigenvalues[1] += T((4. * GV / 3. - KV * 2.) * sinpsi) * gamma[0] + T(2. * GV * (1. - sinpsi / 3.) - 2. * KV * sinpsi) * gamma[1];
    eigenvalues[2] += T(2. * GV * (1 - sinpsi / 3.) - 2. * KV * sinpsi) * gamma[0] + T((4. * GV / 3. - 2. * KV) * sinpsi) * gamma[1];
    sigma_projected = eigenvalues;
    epsbarnew = TPZExtractVal::val(epsbar);

    bool check_validity_Q = (TPZExtractVal::val(eigenvalues[0]) > TPZExtractVal::val(eigenvalues[1]) || IsZero(eigenvalues[0]-eigenvalues[1])) && (TPZExtractVal::val(eigenvalues[1]) > TPZExtractVal::val(eigenvalues[2]) || IsZero(eigenvalues[1]-eigenvalues[2]));
    return (check_validity_Q);
}

/**
 * @brief Computes dsigmapr/dsigmatr for the ReturnMapRightEdge
 */
void TPZYCMohrCoulombPV::ComputeRightEdgeTangent(TPZMatrix<REAL> &tang, REAL &epsbarp) const {
    
    const REAL sin_phi = sin(fPhi);
    const REAL sin_psi = sin(fPsi);
    const REAL G = fER.G(), K = fER.K();
    const REAL a = 4.0 * G * (1.0 + (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;
    const REAL b = 2.0 * G * (1.0 + sin_phi + sin_psi - (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;

    REAL epsbar = epsbarp;
    REAL c, H;
    PlasticityFunction(epsbar, c, H);
    
    tang.Redim(3, 3);
    
    // First column
    tang(0, 0) = (3.0*a + 3.0*b - 4.0*(1.0 + sin_phi)*(3.0*K*sin_psi + G*(3.0 + sin_psi)))/(3.*(a + b));
    tang(1, 0) = (2.0*(1.0 + sin_phi)*(-6.0*K*sin_psi + G*(3.0 + sin_psi)))/(3.*(a + b));
    tang(2, 0) = (2.0*(1.0 + sin_phi)*(-6.0*K*sin_psi + G*(3.0 + sin_psi)))/(3.*(a + b));
    
    // Second column
    tang(0, 1) = (-2*(-1 + sin_phi)*(3*K*sin_psi + G*(3 + sin_psi)))/(3.*(a + b));
    tang(1, 1) = (-3*b*b + 3*a*(a + 2*G*(-1 + sin_phi)) -
                  2*(a*G + 2*b*G + 3*a*K - 3*b*K)*(-1 + sin_phi)*sin_psi)/(3.*(a - b)*(a + b));
    tang(2, 1) = (2*(-1 + sin_phi)*(b*G*(-3 + sin_psi) + a*(2*G - 3*K)*sin_psi + 3*b*K*sin_psi))/
    (3.*(a*a - b*b));
    
    // Third column
    tang(0, 2) = (-2*(-1 + sin_phi)*(3*K*sin_psi + G*(3 + sin_psi)))/(3.*(a + b));
    tang(1, 2) = (2*(-1 + sin_phi)*(b*G*(-3 + sin_psi) + a*(2*G - 3*K)*sin_psi + 3*b*K*sin_psi))/
    (3.*(a*a - b*b));
    tang(2, 2) = (-3*b*b + 3*a*(a + 2*G*(-1 + sin_phi)) -
                  2*(a*G + 2*b*G + 3*a*K - 3*b*K)*(-1 + sin_phi)*sin_psi)/(3.*(a - b)*(a + b));
    
}

template<class T>
bool TPZYCMohrCoulombPV::ReturnMapApex(const TPZVec<T> &sigmatrial, TPZVec<T> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew) const {
    
    const REAL K = fER.K();
    const REAL sinpsi = sin(fPsi);
    const REAL cosphi = cos(fPhi);
    const REAL cotphi = 1. / tan(fPhi);
    T ptrnp1 = 0.;
    for (int i = 0; i < 3; i++) {
        ptrnp1 += T(sigmatrial[i]);
    }
    ptrnp1 /= 3.;
    T DEpsPV = 0.;
    T epsbarnp1 = T(fEpsPlasticBar);
    T c, H;
    PlasticityFunction(epsbarnp1, c, H);

    T alpha = cos(fPhi) / sin(fPsi);
    REAL tol = 1.e-8;

    T res = c * cotphi - ptrnp1;
    T pnp1;
    
    int n_iterations = 30; // @TODO : Define a numeric controls manager object and use it to obtain this information
    int i;
    bool stop_criterion;
    for (i = 0; i < n_iterations; i++) {
        const T jac = H * T(cosphi * cotphi) / T(sinpsi) + T(K);
        DEpsPV -= res / jac;

        epsbarnp1 = T(fEpsPlasticBar) + T(alpha) * DEpsPV;
        pnp1 = ptrnp1 - T(K) * DEpsPV;
        PlasticityFunction(epsbarnp1, c, H);
        res = c * cotphi - pnp1;
        
        stop_criterion = IsZero(res);
        if (stop_criterion) {
            break;
        }
    }
    
#ifdef PZDEBUG
    if (i == n_iterations) {
        DebugStop();
    }
#endif
    
    epsbarnew = TPZExtractVal::val(epsbarnp1);
    for (int i = 0; i < 3; i++) {
        sigma_projected[i] = pnp1;
    }
    return true;
}

void TPZYCMohrCoulombPV::ComputeApexGradient(TPZMatrix<REAL> & gradient, REAL & eps_bar_p) const {
    
    REAL c, H;
    const REAL cosphi = cos(fPhi);
    const REAL sinpsi = sin(fPsi);
    const REAL cotphi = 1. / tan(fPhi);
    const REAL K = fER.K();
    const REAL alpha = cosphi / sinpsi;
    PlasticityFunction(eps_bar_p, c, H);
    const REAL num = H * alpha * cotphi / K;
    const REAL denom = 1. + num;
    const REAL dpdptr = num / denom;
    const REAL dsigdsigtr = dpdptr / 3.;
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            gradient(i, j) = dsigdsigtr;
        }
    }
    
}

void TPZYCMohrCoulombPV::ProjectSigma(const TPZVec<STATE> & sigma_trial, STATE k_prev, TPZVec<STATE> & sigma, STATE &k_proj, int & m_type, TPZFMatrix<REAL> * gradient) {
    
    bool require_gradient_Q = true;
    if (!gradient) {
        require_gradient_Q = false;
    }
    
#ifdef PZDEBUG
    if (require_gradient_Q) {
        // Check for required dimensions of tangent
        if (!(gradient->Rows() == 3 && gradient->Cols() == 3)) {
            std::cerr << "Unable to compute the gradient operator. Required gradient array dimensions are 3x3." << std::endl;
            DebugStop();
        }
    }
#endif
    
    TComputeSequence memory;
    this->SetEpsBar(k_prev);
    REAL epsbartemp = k_prev; // it will be defined by the correct returnmap
    
    bool check_validity_Q;
#ifdef PZDEBUG
    // Check if we are in the correct sextant
    check_validity_Q = (TPZExtractVal::val(sigma_trial[0]) > TPZExtractVal::val(sigma_trial[1]) || IsZero(sigma_trial[0]-sigma_trial[1])) && (TPZExtractVal::val(sigma_trial[1]) > TPZExtractVal::val(sigma_trial[2]) || IsZero(sigma_trial[1]-sigma_trial[2]));
    if (!check_validity_Q) {
        DebugStop();
    }
#endif

    REAL phi = PhiPlane<REAL>(sigma_trial);
    bool elastic_update_Q = IsZero(phi) || phi < 0.0;
    if (elastic_update_Q) {
        m_type = 0; // Elastic behavior
        memory.fWhichPlane = TComputeSequence::EElastic;
        memory.fGamma.Resize(0);
        sigma = sigma_trial;
        
        if (require_gradient_Q) {
            gradient->Identity();
        }
        return;
    }
    
    m_type = 1; // failure behavior
    TPZVec<REAL> sigma_projected;
    memory.fGamma.Resize(1);
    memory.fGamma[0] = 0.;
    check_validity_Q = this->ReturnMapPlane<REAL>(sigma_trial, sigma_projected, memory, epsbartemp);
    if (check_validity_Q) {
        k_proj = epsbartemp;
        this->SetEpsBar(k_proj);
        sigma = sigma_projected;
        memory.fWhichPlane = TComputeSequence::EMainPlane;
        
        if (require_gradient_Q) {
            ComputePlaneTangent(*gradient, epsbartemp);
        }
    } else {
        memory.fGamma.Resize(2);
        memory.fGamma[0] = 0.;
        memory.fGamma[1] = 0.;
        bool IsEdge = false;

        const REAL sinpsi = sin(fPsi);
        REAL val = (1 - sinpsi) * sigma_trial[0] - 2. * sigma_trial[1] + (1 + sinpsi) * sigma_trial[2];
        if (val > 0.) {
            IsEdge = this->ReturnMapRightEdge<REAL>(sigma_trial, sigma_projected, memory, epsbartemp);
            memory.fWhichPlane = TComputeSequence::ERightEdge;
            
            if (require_gradient_Q) {
                ComputeRightEdgeTangent(*gradient, epsbartemp);
            }
        } else {
            IsEdge = this->ReturnMapLeftEdge<REAL>(sigma_trial, sigma_projected, memory, epsbartemp);
            memory.fWhichPlane = TComputeSequence::ELeftEdge;
            
            if (require_gradient_Q) {
                ComputeLeftEdgeTangent(*gradient, epsbartemp);
            }
        }
        if (!IsEdge) {
            m_type = -1; // Tensile behavior
            this->ReturnMapApex(sigma_trial, sigma_projected, memory, epsbartemp);
            memory.fWhichPlane = TComputeSequence::EApex;
            
            if (require_gradient_Q) {
                ComputeApexGradient(*gradient, epsbartemp);
            }
        }

        k_proj = epsbartemp;
        this->SetEpsBar(k_proj);
        sigma = sigma_projected;
    }
}

template void TPZYCMohrCoulombPV::PlasticityFunction<REAL>(const REAL epsp, REAL &c, REAL &H) const;
template void TPZYCMohrCoulombPV::PlasticityFunction<fadtype>(const fadtype epsp, fadtype &c, fadtype &H) const;

template TPZVec<REAL> TPZYCMohrCoulombPV::SigmaElastPV<REAL>(const TPZVec<REAL> &deform) const;
template TPZVec<fadtype> TPZYCMohrCoulombPV::SigmaElastPV<fadtype>(const TPZVec<fadtype> &deform) const;

template REAL TPZYCMohrCoulombPV::PhiPlane<REAL>(const TPZVec<REAL> &sigma) const;
template fadtype TPZYCMohrCoulombPV::PhiPlane<fadtype>(const TPZVec<fadtype> &sigma) const;

template bool TPZYCMohrCoulombPV::ReturnMapPlane<REAL>(const TPZVec<REAL> &sigma_trial, TPZVec<REAL> &sigma_projected,
        TPZYCMohrCoulombPV::TComputeSequence &memory, REAL &epsbarnew) const;
template bool TPZYCMohrCoulombPV::ReturnMapPlane<fadtype>(const TPZVec<fadtype> &sigma_trial, TPZVec<fadtype> &sigma_projected,
        TPZYCMohrCoulombPV::TComputeSequence &memory, REAL &epsbarnew) const;

template bool TPZYCMohrCoulombPV::ReturnMapLeftEdge<REAL>(const TPZVec<REAL> &sigma_trial, TPZVec<REAL> &sigma_projected,
        TPZYCMohrCoulombPV::TComputeSequence &memory, REAL &epsbarnew) const;
template bool TPZYCMohrCoulombPV::ReturnMapLeftEdge<fadtype>(const TPZVec<fadtype> &sigma_trial, TPZVec<fadtype> &sigma_projected,
        TPZYCMohrCoulombPV::TComputeSequence &memory, REAL &epsbarnew) const;

template bool TPZYCMohrCoulombPV::ReturnMapRightEdge<REAL>(const TPZVec<REAL> &sigma_trial, TPZVec<REAL> &sigma_projected,
        TPZYCMohrCoulombPV::TComputeSequence &memory, REAL &epsbarnew) const;
template bool TPZYCMohrCoulombPV::ReturnMapRightEdge<fadtype>(const TPZVec<fadtype> &sigma_trial, TPZVec<fadtype> &sigma_projected,
        TPZYCMohrCoulombPV::TComputeSequence &memory, REAL &epsbarnew) const;

template bool TPZYCMohrCoulombPV::ReturnMapApex<REAL>(const TPZVec<REAL> &sigma_trial, TPZVec<REAL> &sigma_projected,
        TPZYCMohrCoulombPV::TComputeSequence &memory, REAL &epsbarnew) const;
template bool TPZYCMohrCoulombPV::ReturnMapApex<fadtype>(const TPZVec<fadtype> &sigma_trial, TPZVec<fadtype> &sigma_projected,
        TPZYCMohrCoulombPV::TComputeSequence &memory, REAL &epsbarnew) const;
