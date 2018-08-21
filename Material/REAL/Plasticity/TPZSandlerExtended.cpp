//
//  pzsandlerextPV.cpp
//  PZ
//
//  Created by Diogo Cecilio on 9/3/13.
//
//

#include "TPZSandlerExtended.h"
#include "pzlog.h"
#include "pzreferredcompel.h"
#include "TPZHWTools.h"
#include "TPZConvergenceException.h"
#include "TPZInconsistentStateException.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.poroelastoplastic"));
#endif

#ifdef LOG4CXX
static LoggerPtr loggerConvTest(Logger::getLogger("ConvTest"));
#endif

TPZSandlerExtended::TPZSandlerExtended() : ftol(1e-5), fA(0), fB(0), fC(0), fD(0), fW(0), fK(0), fR(0), fG(0), fPhi(0), fN(0), fPsi(0), fE(0), fnu(0), fkappa_0(0) {
}

TPZSandlerExtended::TPZSandlerExtended(const TPZSandlerExtended & copy) {
    ftol = copy.ftol;
    fA = copy.fA;
    fB = copy.fB;
    fC = copy.fC;
    fD = copy.fD;
    fW = copy.fW;
    fK = copy.fK;
    fR = copy.fR;
    fG = copy.fG;
    fPhi = copy.fPhi;
    fN = copy.fN;
    fPsi = copy.fPsi;
    fE = copy.fE;
    fnu = copy.fnu;
    fkappa_0 = copy.fkappa_0;
    fElasticResponse = copy.fElasticResponse;
}

TPZSandlerExtended::TPZSandlerExtended(STATE A, STATE B, STATE C, STATE D, STATE K, STATE G, STATE W, STATE R, STATE Phi, STATE N, STATE Psi, STATE kappa_0) :
fA(A), fB(B), fC(C), fD(D), fW(W), fK(K), fR(R), fG(G), fPhi(Phi), fN(N), fPsi(Psi), fkappa_0(kappa_0) {
    fE = (9. * fK * fG) / (3. * fK + fG);
    fnu = ((3. * fK)-(2. * fG)) / (2 * (3. * fK + fG));
    TPZElasticResponse ER;
    ER.SetUp(fE, fnu);
    fElasticResponse = ER;
    ftol = 1.e-10;
}

TPZSandlerExtended::~TPZSandlerExtended() {

}

template <class T>
T TPZSandlerExtended::F(const T x) const {
    return (fA - fC * exp(x * fB) - fPhi * x);
}

template <class T>
T TPZSandlerExtended::DF(const T x) const {
    return -(exp(fB * x) * fB * fC) - fPhi;
}

STATE TPZSandlerExtended::GetF(STATE x) const {
    return F(x);
}

template<class T>
T TPZSandlerExtended::X(const T L) const {
    return (L - fR * F(L));
}

STATE TPZSandlerExtended::GetX(STATE k) {
    return X(k);
}

void TPZSandlerExtended::SetUp(STATE A, STATE B, STATE C, STATE D, STATE K, STATE G, STATE W, STATE R, STATE Phi, STATE N, STATE Psi) {
    fA = A;
    fB = B;
    fC = C;
    fD = D;
    fK = K;
    fG = G;
    fW = W;
    fR = R;
    fPhi = Phi;
    fN = N;
    fPsi = Psi;
    fE = (9. * fK * fG) / (3. * fK + fG);
    fnu = ((3. * fK)-(2. * fG)) / (2 * (3. * fK + fG));
    TPZElasticResponse ER;
    ER.SetUp(fE, fnu);
    fElasticResponse = ER;

}

void TPZSandlerExtended::SetInitialDamage(STATE kappa_0) {
    fkappa_0 = kappa_0;
}

void TPZSandlerExtended::SetElasticResponse(const TPZElasticResponse &ER) {
    fElasticResponse = ER;
    fE = ER.E();
    fnu = ER.Poisson();
    fK = ER.K();
    fG = ER.G();
}

TPZElasticResponse TPZSandlerExtended::GetElasticResponse() const {
    return fElasticResponse;
}

void TPZSandlerExtended::Read(TPZStream& buf, void* context) { //ok
    buf.Read(&ftol);
    buf.Read(&fA);
    buf.Read(&fB);
    buf.Read(&fC);
    buf.Read(&fD);
    buf.Read(&fW);
    buf.Read(&fK);
    buf.Read(&fR);
    buf.Read(&fG);
    buf.Read(&fPhi);
    buf.Read(&fN);
    buf.Read(&fPsi);
    buf.Read(&fE);
    buf.Read(&fnu);
    buf.Read(&fkappa_0);
    fElasticResponse.Read(buf, context);
}

void TPZSandlerExtended::Write(TPZStream& buf, int withclassid) const { //ok
    buf.Write(&ftol);
    buf.Write(&fA);
    buf.Write(&fB);
    buf.Write(&fC);
    buf.Write(&fD);
    buf.Write(&fW);
    buf.Write(&fK);
    buf.Write(&fR);
    buf.Write(&fG);
    buf.Write(&fPhi);
    buf.Write(&fN);
    buf.Write(&fPsi);
    buf.Write(&fE);
    buf.Write(&fnu);
    buf.Write(&fkappa_0);
    fElasticResponse.Write(buf, withclassid);
}

TPZElasticResponse TPZSandlerExtended::GetElasticResponse() {
    return fElasticResponse;
}

STATE TPZSandlerExtended::GetR() {
    return fR;
}

template<class T>
T TPZSandlerExtended::EpsEqX(T X) const {
    STATE CPer = CPerturbation();
    return (fW * (exp(fD * X) - 1 + CPer * X));
}

template<class T>
T TPZSandlerExtended::EpsEqk(const T k) const {
    return EpsEqX(X(k));
}

void TPZSandlerExtended::Firstk(STATE &epsp, STATE &k) const {
    STATE f, df, kn1, kn, resnorm, diff;
    int counter = 1;
    resnorm = 1;
    kn = epsp; //chute inicial
    kn1 = kn;
    while (resnorm > ftol && counter < 30) {

        f = EpsEqk(kn) - epsp;
        df = fD * exp(fD * (kn - (fA - fC * exp(fB * kn)) * fR))*(1 + fB * fC * exp(fB * kn) * fR) * fW;
        //df=fD*exp(fD*(kn - fR*(fA - fC*exp(fB*kn) - kn*fPhi)))*fW*(1 - fR*(-(fB*fC*exp(fB*kn)) - fPhi));
        kn1 = kn - f / df;
        diff = kn1 - kn;
        resnorm = sqrt(diff * diff);
        kn = kn1;
        counter++;

    }
    k = kn1;
}

/// Compute the normal function to the failure surface based on a reference point (I1_ref,f1(I1_ref))
STATE TPZSandlerExtended::NormalToF1(STATE I1, STATE I1_ref) const {
    
#ifdef PZDEBUG
    if (I1 < I1_ref) { // normal function is constructed for  I1 >= I1_ref
        DebugStop();
    }
#endif
    
    STATE normal_f1 = (exp(fB*I1_ref)*fA*fB*fC - exp(2*fB*I1_ref)*fB*pow(fC,2) + I1 - I1_ref)/(exp(fB*I1_ref)*fB*fC);
    return  normal_f1;
}

REAL TPZSandlerExtended::InitialDamage(const TPZVec<REAL> &stress_pv) const {
    
    TPZManVector<REAL,2> f(2);
    YieldFunction(stress_pv, 0.0, f);
    
    bool Is_valid_stress_Q = fabs(f[0]) < ftol || f[0] < 0.0;

    if (Is_valid_stress_Q) {
        
        int n_iter = 30; // @TODO:: Variable to GUI.
        REAL I1 = (stress_pv[0])+(stress_pv[1])+(stress_pv[2]);
        REAL res, jac, dk, k;
        
        k = 0.0; // initial guess
        bool stop_criterion_Q = false;
        int i;
        for (i = 0; i < n_iter; i++) {
            res = I1 - X(k);
            stop_criterion_Q = fabs(res) < ftol;
            if (stop_criterion_Q) {
                break;
            }
            jac = - 1.0 - fB * fC * fR * exp(fB*k);
            dk =  - res /jac;
            k+=dk;
        }
        
        if (!stop_criterion_Q) {
            throw TPZConvergenceException(ftol, n_iter, res, i, "TPZSandlerExtended::InitialDamage:: Newton process did not converge in hydrostatic direction.");
        }
        
        stop_criterion_Q = false;
        REAL J2 = (1.0/3.0) * (stress_pv[0]*stress_pv[0] + stress_pv[1]*stress_pv[1] + stress_pv[2]*stress_pv[2] - stress_pv[1]*stress_pv[2] - stress_pv[0]*stress_pv[2] - stress_pv[0]*stress_pv[1]);
        
        k = X(k) + k; // guess from the outer part of the cap
        
        for (int i = 0; i < n_iter; i++) {
            
            res = -1 + pow(I1,2)/(pow(fA - exp(fB*k)*fC,2)*pow(fR,2)) +
            J2/pow(fA - exp(fB*k)*fC,2) -
            (2*I1*k)/(pow(fA - exp(fB*k)*fC,2)*pow(fR,2)) +
            pow(k,2)/(pow(fA - exp(fB*k)*fC,2)*pow(fR,2));
            
            stop_criterion_Q = fabs(res) < ftol;
            if (stop_criterion_Q) {
                break;
            }
            jac = (2*(fA*(-I1 + k) + exp(fB*k)*fC*
                      (I1 + fB*pow(I1,2) + fB*pow(fR,2)*J2 - k - 2*fB*I1*k + fB*pow(k,2))))/
            (pow(fA - exp(fB*k)*fC,3)*pow(fR,2));
            dk =  - res /jac;
            k+=dk;
        }
        
        if (!stop_criterion_Q) {
            throw TPZConvergenceException(ftol, n_iter, res, i, "TPZSandlerExtended::InitialDamage:: Newton process did not converge in deviatoric direction.");
        }
        
        YieldFunction(stress_pv, k, f);
        bool Is_valid_stress_on_cap_Q =  fabs(f[1]) < ftol || f[1] < 0.0;
        
        if (!Is_valid_stress_on_cap_Q) {
            throw TPZInconsistentStateException("TPZSandlerExtended::InitialDamage: Invalid stress state over cap.");
        }
        return k;
        
    }
    else{
        throw TPZInconsistentStateException("TPZSandlerExtended::InitialDamage: Invalid stress state over failure surface.");
    }
    
    return -1;
    
}

template<class T>
T TPZSandlerExtended::ResLF2(const TPZVec<T> &trial_stress, T theta, T beta, T k, STATE kprev) const {
    T trial_I1 = (trial_stress[0])+(trial_stress[1])+(trial_stress[2]);
    T I1 = fR * F(k) * cos(theta) + k;
    T delepsp = EpsEqk(k) - EpsEqk(kprev);
    return (3. * fK * delepsp - (trial_I1 - I1));
}

template<class T>
T TPZSandlerExtended::ResLF2IJ(const TPZVec<T> &sigtrialIJ, T theta, T k, STATE kprev) const {

    T I1tr = sigtrialIJ[0];
    T I1 = fR * F(k) * cos(theta) + k;
    T delepsp = EpsEqk(k) - EpsEqk(kprev);
    return (3. * fK * delepsp - (I1tr - I1));
}

/// Compute the residual of the equation which defines the update of the damage variable

STATE TPZSandlerExtended::ResLF1(const TPZVec<STATE> &sigtrial, const TPZVec<STATE> &sigproj, const STATE k, const STATE kprev) const {
    STATE I1trial = (sigtrial[0] + sigtrial[1] + sigtrial[2]);
    STATE I1proj = (sigproj[0] + sigproj[1] + sigproj[2]);
    STATE delepsp = EpsEqk(k) - EpsEqk(kprev);
    return (3. * fK * delepsp - (I1trial - I1proj));

}

/// Compute the derivative of the equation which determines the evolution of k
// the derivative are given in terms of k

STATE TPZSandlerExtended::DResLF1(const TPZVec<STATE> &sigtrial, const TPZVec<STATE> &sigproj, const STATE k, const STATE kprev) const {
    STATE expfBk = exp(fB * k);
    STATE dreskk = 3. * exp(fD * (-((fA - expfBk * fC) * fR) + k)) * fD * fK *
            (1. + expfBk * fB * fC * fR) * fW;
    return dreskk;
}


/// Compute the derivative of the equation which determines the evolution of k

void TPZSandlerExtended::DResLF2(const TPZVec<STATE> &pt, STATE theta, STATE beta, STATE k, STATE kprev, TPZVec<STATE> &dresl) const {
    STATE expfBk = exp(fB * k);
    STATE sintheta = sin(theta);
    STATE costheta = cos(theta);

    STATE dreskk = 1. + costheta * (-(expfBk * fB * fC) - fPhi) * fR +
            3. * exp(fD * (-((fA - expfBk * fC) * fR) + k)) * fD * fK *
            (1. + expfBk * fB * fC * fR) * fW;


    STATE dresktheta = -fR * (fA - fC * expfBk - k * fPhi) * sintheta;
    dresl[0] = dresktheta;
    dresl[1] = 0;
    dresl[2] = dreskk;


}

void TPZSandlerExtended::F1Cyl(STATE xi, STATE beta, TPZVec<STATE> &f1cyl) const {
    STATE sqrt2 = M_SQRT2;
    STATE sqrt3 = sqrt(3.);
    STATE gamma = 0.5 * (1 + (1 - sin(3 * beta)) / fPsi + sin(3 * beta));
    STATE I1 = xi*sqrt3;
    STATE F1 = F(I1);
    STATE sqrtj2 = F1 / gamma;
    STATE rho = sqrt2*sqrtj2;
    f1cyl[0] = xi;
    f1cyl[1] = rho;
    f1cyl[2] = beta;

}

void TPZSandlerExtended::SurfaceParamF1(TPZVec<STATE> &sigproj, STATE &xi, STATE &beta) const {
    TPZManVector<STATE> sigHWCyl(3);
    TPZHWTools::FromPrincipalToHWCyl(sigproj, sigHWCyl);
    xi = sigHWCyl[0];
    beta = sigHWCyl[2];
#ifdef PZDEBUG
    STATE dist = DistF1(sigproj, xi, beta);
    if (fabs(dist) > ftol) {
        DebugStop();
    }
#endif
}

void TPZSandlerExtended::F2Cyl(STATE theta, STATE beta, STATE k, TPZVec<STATE> &f2cyl) const {
    const STATE M_SQRT3 = sqrt(3.);
    const STATE gamma = 0.5 * (1.0 + sin(3.0 * beta) + (1.0 - sin(3.0 * beta)) / fPsi);
    const STATE Fk = F(k);
    const STATE var = fR * Fk * cos(theta);
    const STATE sig1_star = (k - var) / M_SQRT3; // The definition for theta was corrected.
    const STATE sqrtj2 = Fk * sin(theta) / gamma;
    const STATE rho = M_SQRT2*sqrtj2;
    const STATE xi = sig1_star;
    f2cyl[0] = xi;
    f2cyl[1] = rho;
    f2cyl[2] = beta;

}

void TPZSandlerExtended::SurfaceParamF2(const TPZVec<STATE> &sigproj, const STATE k, STATE &theta, STATE &beta) const {
    TPZManVector<STATE> sigHWCyl(3);
    TPZHWTools::FromPrincipalToHWCyl(sigproj, sigHWCyl);
    //    STATE xi,rho;
    //    xi=sigHWCyl[0];
    //    rho=sigHWCyl[1];
    STATE I1 = sigHWCyl[0] * sqrt(3.);
    beta = sigHWCyl[2];
    STATE Fk = F(k);
    STATE gamma = 0.5 * (1 + (1 - sin(3 * beta)) / fPsi + sin(3 * beta));
    STATE costheta = (I1 - k) / (fR * Fk);
    STATE sqrtj2 = sigHWCyl[1] / sqrt(2.);
    STATE sintheta = sqrtj2 * gamma / (Fk - fN);
    theta = atan2(sintheta, costheta);
    //theta = acos(costheta);
    //STATE theta2 = atan((rho*sin(beta))/xi);
#ifdef PZDEBUG
    STATE err = 1. - sintheta * sintheta - costheta*costheta;
    STATE dist = DistF2(sigproj, theta, beta, k);
    if (fabs(dist) > ftol || err > ftol) {
        DebugStop();
    }
#endif
}

STATE TPZSandlerExtended::DistF1(const TPZVec<STATE> &pt, STATE xi, STATE beta) const {
    TPZManVector<STATE, 3> cyl(3);
    F1Cyl(xi, beta, cyl);
    TPZManVector<STATE, 3> cart(3);
    TPZHWTools::FromHWCylToHWCart(cyl, cart);
    TPZManVector<STATE, 3> carttrial(3);
    TPZHWTools::FromPrincipalToHWCart(pt, carttrial);
    return ((1. / (3. * fK))*(carttrial[0] - cart[0])*(carttrial[0] - cart[0]))
            +(1. / (2. * fG))*((carttrial[1] - cart[1])*(carttrial[1] - cart[1])+(carttrial[2] - cart[2])*(carttrial[2] - cart[2]));

}

STATE TPZSandlerExtended::DistF2(const TPZVec<STATE> &pt, STATE theta, STATE beta, STATE k) const {
    TPZManVector<STATE, 3> cart_trial(3); // cyl and cart trial are the same variable, it is renamed as cart_trial
    TPZManVector<STATE, 3> cart(3);
    F2Cyl(theta, beta, k, cart_trial);
    TPZHWTools::FromHWCylToHWCart(cart_trial, cart);
    TPZHWTools::FromPrincipalToHWCart(pt, cart_trial);
    return ((1. / (3. * fK))*(cart_trial[0] - cart[0])*(cart_trial[0] - cart[0]))
            +(1. / (2. * fG))*((cart_trial[1] - cart[1])*(cart_trial[1] - cart[1])+(cart_trial[2] - cart[2])*(cart_trial[2] - cart[2]));
    
}

STATE TPZSandlerExtended::DistF2IJ(const TPZVec<STATE> &sigtrialIJ, STATE theta, STATE k) const {
    STATE I1 = sigtrialIJ[0];
    STATE sqJ2 = sigtrialIJ[1];
    STATE Fk;
    Fk = F(k);
    STATE y = (sqJ2 - Fk * sin(theta));
    STATE x = 1. / (3 * fK)*(I1 - (k + Fk * fR * cos(theta)));
    STATE res = x * x / (9. * fK) + y * y / (fG);
    return res;

}

template<class T>
void TPZSandlerExtended::FromThetaKToSigIJ(const T &theta, const T &K, TPZVec<T> &sigIJ) const {
    T Fk = F(K);
    sigIJ[0] = K + Fk * fR * cos(theta);
    sigIJ[1] = Fk * sin(theta);
}

/**
 * compute the value of the equation which determines the orthogonality of the projection
 */
template<class T>
void TPZSandlerExtended::DDistF2IJ(TPZVec<T> &sigtrialIJ, T theta, T L, STATE LPrev, TPZVec<T> &ddistf2) const {
    T I1 = sigtrialIJ[0];
    T sqJ2 = sigtrialIJ[1];
    T Fk;
    Fk = F(L);
    T y = (sqJ2 - Fk * sin(theta));
    T x = (I1 - (L + Fk * fR * cos(theta)));
    ddistf2[0] = T(2.) * x * Fk * fR * sin(theta) / T(9. * fK) - T(2.) * y * Fk * cos(theta);
    ddistf2[1] = ResLF2IJ(sigtrialIJ, theta, L, LPrev);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "x = " << x << " y = " << y << " theta = " << theta << " res = " << ddistf2[0];
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

void TPZSandlerExtended::DDistFunc1(const TPZVec<STATE> &pt, STATE xi, STATE beta, TPZFMatrix<STATE> &ddistf1) const {
    STATE sigstar1, sigstar2, sigstar3, DFf, Ff, I1, sb, cb, DGamma;
    TPZManVector<STATE,3> ptcart(3);
    sb = sin(beta);
    cb = cos(beta);
    STATE sin3b = sin(3 * beta);
    STATE cos3b = cos(3 * beta);
    TPZHWTools::FromPrincipalToHWCart(pt, ptcart);
    sigstar1 = ptcart[0];
    sigstar2 = ptcart[1];
    sigstar3 = ptcart[2];
    I1 = xi * sqrt(3);
    Ff = F(I1);
    const REAL gamma = (1 + sin3b + (1 - sin3b) / fPsi) / 2.;
    const REAL gamma2 = gamma*gamma;
    const REAL gamma3 = gamma*gamma2;
    DFf = DF(I1);
    DGamma = 0.5*(3 * cos3b - (3 * cos3b) / fPsi);
    const REAL sqrt2 = M_SQRT2;
    const REAL sqrt3 = sqrt(3.);

    ddistf1.Resize(2, 1);
    ddistf1(0, 0) = (6 * sqrt3 * DFf * Ff * fK - 3 * sqrt3 * sqrt2 * DFf * fK * gamma * (cb * sigstar2 + sb * sigstar3) + 2 * fG * gamma2 * (-sigstar1 + xi)) / (3. * fG * fK * gamma2);
    ddistf1(1, 0) = -((Ff * sqrt2 * (gamma2 * (-(sb * sigstar2) + cb * sigstar3) - DGamma * gamma * (cb * sigstar2 + sb * sigstar3) + DGamma * Ff * sqrt2)) / (fG * gamma3));

}

void TPZSandlerExtended::D2DistFunc1(const TPZVec<STATE> &pt, STATE xi, STATE beta, TPZFMatrix<STATE> &tangentf1) const {
    STATE sig2, sig3, DFf, Gamma, Ff, I1, sb, cb, D2Ff, DGamma, D2Gamma, Gamma2, Gamma3, Sqrt2, Sqrt3;
    TPZVec<STATE> ptcart(3);
    sb = sin(beta);
    cb = cos(beta);
    STATE sin3b = sin(3 * beta);
    STATE cos3b = cos(3 * beta);
    TPZHWTools::FromPrincipalToHWCart(pt, ptcart);
    //    STATE sig1 = ptcart[0];
    sig2 = ptcart[1];
    sig3 = ptcart[2];
    I1 = xi * sqrt(3);
    Ff = F(I1);
    Gamma = (1. + sin3b + (1. - sin3b) / fPsi) / 2.;
    DFf = -DF(I1);
    D2Ff = -(exp(fB * I1) * pow(fB, 2.) * fC);
    DGamma = (3. * cos3b - (3. * cos3b) / fPsi) / 2.;
    D2Gamma = (-9. * sin(3. * beta) + (9. * sin(3. * beta)) / fPsi) / 2.;
    Gamma2 = Gamma*Gamma;
    Gamma3 = Gamma*Gamma2;
    //    STATE D2Gamma2=D2Gamma*D2Gamma;
    Sqrt2 = sqrt(2);
    Sqrt3 = sqrt(3);

    tangentf1.Resize(2, 2);


    tangentf1(0, 0) = (18 * (DFf * DFf + D2Ff * Ff) * fK + 2 * fG * Gamma2 -
            9 * D2Ff * fK * Gamma * (cb * sig2 + sb * sig3) * Sqrt2) / (3. * fG * fK * Gamma2);

    tangentf1(0, 1) = -((Sqrt3 * Sqrt2 * DFf * (Gamma2 * (-(sb * sig2) + cb * sig3) - DGamma * Gamma * (cb * sig2 + sb * sig3) +
            2 * DGamma * Ff * Sqrt2)) / (fG * Gamma3));

    tangentf1(1, 1) = -((Ff * (-6 * DGamma * DGamma * Ff - Gamma3 * (cb * sig2 + sb * sig3) * Sqrt2 -
            Gamma2 * (2 * DGamma * (-(sb * sig2) + cb * sig3) + D2Gamma * (cb * sig2 + sb * sig3)) *
            Sqrt2 + 2 * Gamma * (D2Gamma * Ff + DGamma * DGamma * (cb * sig2 + sb * sig3) * Sqrt2))) /
            (fG * Gamma2 * Gamma2));

    tangentf1(1, 0) = tangentf1(0, 1);

}


/// Compute the derivative of the distance function to the failure function and the result of Residue 1 (failure)
template<class T>
void TPZSandlerExtended::Res1(const TPZVec<T> &trial_stress, T i1, T beta, T k, T kprev, TPZVec<T> & residue_1) const{
    
    // In this implementation the definition for theta is given by the angle formed from -I1 axis to sqrt(J2) axis with origin on the damage variable kappa.
    
    residue_1.Resize(3, 1);
    STATE CX0   = X_0();
    STATE CPer  = CPerturbation();
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE CG    = fE/(2.0*(1.0 + fnu));
    
    TPZManVector<REAL,3> rhw_sigma(3);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_sigma);
    
    residue_1[0] = (2*(i1 - sqrt(3)*rhw_sigma[0]))/(9.*CK) + (2*sqrt(2)*exp(fB*i1)*fB*fC*fPsi*cos(beta)*
                                               (-2*sqrt(2)*(fA - exp(fB*i1)*fC)*fPsi*cos(beta) + rhw_sigma[1]*(1 + fPsi + (-1 + fPsi)*sin(3*beta)))
                                               )/(CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2)) +
    (2*sqrt(2)*exp(fB*i1)*fB*fC*fPsi*sin(beta)*
     (-2*sqrt(2)*(fA - exp(fB*i1)*fC)*fPsi*sin(beta) + rhw_sigma[2]*(1 + fPsi + (-1 + fPsi)*sin(3*beta)))
     )/(CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    
    
    residue_1[1] = (2*sqrt(2)*(fA - exp(fB*i1)*fC)*fPsi*(cos(beta)*(1 + fPsi + 8*(-1 + fPsi)*pow(sin(beta),3))*
                                          (2*sqrt(2)*(fA - exp(fB*i1)*fC)*fPsi*sin(beta) -
                                           rhw_sigma[2]*(1 + fPsi + (-1 + fPsi)*sin(3*beta))) +
                                          (2*sqrt(2)*(fA - exp(fB*i1)*fC)*fPsi*cos(beta) -
                                           rhw_sigma[1]*(1 + fPsi + (-1 + fPsi)*sin(3*beta)))*
                                          (-3*(-1 + fPsi)*cos(beta)*cos(3*beta) - sin(beta)*(1 + fPsi + (-1 + fPsi)*sin(3*beta)))))/
    (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),3));
    
    residue_1[2] = -i1 + 3*CK*fW*(-exp(-(fD*(CX0 + fA*fR - exp(fB*k)*fC*fR - k))) +
                   exp(-(fD*(CX0 + fA*fR - exp(fB*kprev)*fC*fR - kprev))) - CPer*exp(fB*k)*fC*fR +
                   CPer*exp(fB*kprev)*fC*fR + CPer*(-k + kprev)) + sqrt(3)*rhw_sigma[0];
    
}

template<class T>
void TPZSandlerExtended::Res2(const TPZVec<T> &trial_stress, T theta, T beta, T k, T kprev, TPZVec<T> &residue_2) const {

    // In this implementation the definition for theta is given by the angle formed from -I1 axis to sqrt(J2) axis with origin on the damage variable kappa.
    
    residue_2.Resize(3, 1);
    STATE CX0   = X_0();
    STATE CPer  = CPerturbation();
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE CG    = fE/(2.0*(1.0 + fnu));
    
    TPZManVector<REAL,3> rhw_sigma(3);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_sigma);
    


    residue_2[0] = (2.0*(fA - exp(fB*k)*fC)*(CG*fR*(k - sqrt(3.0)*rhw_sigma[0] - (fA - exp(fB*k)*fC)*fR*cos(theta))*
                                               pow(1.0 + fPsi + (-1.0 + fPsi)*sin(3.0*beta),2.0)*sin(theta) -
                                               9.0*sqrt(2.0)*CK*fPsi*cos(beta)*cos(theta)*
                                               ((1.0 + fPsi)*rhw_sigma[1] + (-1.0 + fPsi)*rhw_sigma[1]*sin(3.0*beta) -
                                                2.0*sqrt(2.0)*(fA - exp(fB*k)*fC)*fPsi*cos(beta)*sin(theta)) -
                                               9.0*sqrt(2.0)*CK*fPsi*cos(theta)*sin(beta)*
                                               ((1 + fPsi)*rhw_sigma[2] + (-1.0 + fPsi)*rhw_sigma[2]*sin(3.0*beta) -
                                                2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*sin(beta)*sin(theta))))/
                                                (9.0*CG*CK*pow(1.0 + fPsi + (-1.0 + fPsi)*sin(3.0*beta),2.0));
    
    residue_2[1] = (-2.0*sqrt(2.0)*(fA - exp(fB*k)*fC)*fPsi*sin(theta)*
                    ((-3.0*(-1.0 + fPsi)*cos(beta)*cos(3.0*beta) - sin(beta)*(1.0 + fPsi + (-1.0 + fPsi)*sin(3.0*beta)))*
                     (rhw_sigma[1]*(1.0 + fPsi + (-1.0 + fPsi)*sin(3.0*beta)) -
                      2.0*sqrt(2.0)*(fA - exp(fB*k)*fC)*fPsi*cos(beta)*sin(theta)) +
                     cos(beta)*(1.0 + fPsi + 8.0*(-1.0 + fPsi)*pow(sin(beta),3.0))*
                     (rhw_sigma[2]*(1.0 + fPsi + (-1.0 + fPsi)*sin(3.0*beta)) -
                      2.0*sqrt(2.0)*(fA - exp(fB*k)*fC)*fPsi*sin(beta)*sin(theta))))/
    (CG*pow(1.0 + fPsi + (-1.0 + fPsi)*sin(3.0*beta),3.0));
    
    
    residue_2[2] = -k + 3.0*CK*fW*(-exp(-(fD*(CX0 + fA*fR - exp(fB*k)*fC*fR - k))) +
                                   exp(-(fD*(CX0 + fA*fR - exp(fB*kprev)*fC*fR - kprev))) - CPer*exp(fB*k)*fC*fR +
                                   CPer*exp(fB*kprev)*fC*fR - CPer*k + CPer*kprev) + sqrt(3.0)*rhw_sigma[0] +
                                    (fA - exp(fB*k)*fC)*fR*cos(theta);
    
}

template<class T>
void TPZSandlerExtended::Res2CoVertex(const TPZVec<T> &trial_stress, T beta, T k, T kprev, TPZVec<T> & residue_covertex) const{
    
    // In this implementation the definition for theta is given by the angle formed from -I1 axis to sqrt(J2) axis with origin on the damage variable kappa.
    
    residue_covertex.Resize(2, 1);
    STATE CX0   = X_0();
    STATE CPer  = CPerturbation();
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE CG    = fE/(2.0*(1.0 + fnu));
    STATE theta = M_PI_2;
    
    TPZManVector<REAL,3> rhw_sigma(3);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_sigma);
    
    residue_covertex[0] = (2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*(cos(beta)*(1 + fPsi + 8*(-1 + fPsi)*pow(sin(beta),3))*
                                         (2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*sin(beta) - rhw_sigma[2]*(1 + fPsi + (-1 + fPsi)*sin(3*beta)))
                                         + (2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*cos(beta) -
                                            rhw_sigma[1]*(1 + fPsi + (-1 + fPsi)*sin(3*beta)))*
                                         (-3*(-1 + fPsi)*cos(beta)*cos(3*beta) - sin(beta)*(1 + fPsi + (-1 + fPsi)*sin(3*beta)))))/
    (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),3));
    
    
    residue_covertex[1] = -k + 3*CK*fW*(-exp(-(fD*(CX0 + fA*fR - exp(fB*k)*fC*fR - k))) +
                  exp(-(fD*(CX0 + fA*fR - exp(fB*kprev)*fC*fR - kprev))) - CPer*exp(fB*k)*fC*fR +
                  CPer*exp(fB*kprev)*fC*fR + CPer*(-k + kprev)) + sqrt(3)*rhw_sigma[0];

}

template<class T>
void TPZSandlerExtended::Res2Vertex(const TPZVec<T> &trial_stress, T k, T kprev, T & residue_vertex) const{
    
    STATE CX0   = X_0();
    STATE CPer  = CPerturbation();
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE theta = 0.0;
    
    TPZManVector<REAL,3> rhw_sigma(3);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_sigma);
    
    residue_vertex = -k + 3*CK*fW*(-exp(-(fD*(CX0 + fA*fR - exp(fB*k)*fC*fR - k))) +
                  exp(-(fD*(CX0 + fA*fR - exp(fB*kprev)*fC*fR - kprev))) - CPer*exp(fB*k)*fC*fR +
                  CPer*exp(fB*kprev)*fC*fR - CPer*k + CPer*kprev) + sqrt(3)*rhw_sigma[0] + (fA - exp(fB*k)*fC)*fR*cos(theta);
}

/// Compute the jacobian function of the f1 (failure) distance as a function of i1, beta and k
void TPZSandlerExtended::Jacobianf1(const TPZVec<STATE> &trial_stress, STATE i1, STATE beta, STATE k, TPZFMatrix<STATE> &jacobianf1)const{
 
    jacobianf1.Resize(3, 3);
    STATE CX0   = X_0();
    STATE CPer  = CPerturbation();
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE CG    = fE/(2.0*(1.0 + fnu));
    
    TPZManVector<REAL,3> rhw_sigma(3);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_sigma);
    
    
    jacobianf1(0,0) = (2*(36*CK*exp(fB*i1)*pow(fB,2)*fC*(-fA + 2*exp(fB*i1)*fC)*pow(fPsi,2)*pow(cos(beta),2) +
                          36*CK*exp(fB*i1)*pow(fB,2)*fC*(-fA + 2*exp(fB*i1)*fC)*pow(fPsi,2)*
                          pow(sin(beta),2) + 9*sqrt(2)*CK*exp(fB*i1)*pow(fB,2)*fC*fPsi*rhw_sigma[1]*cos(beta)*
                          (1 + fPsi + (-1 + fPsi)*sin(3*beta)) +
                          9*sqrt(2)*CK*exp(fB*i1)*pow(fB,2)*fC*fPsi*rhw_sigma[2]*sin(beta)*
                          (1 + fPsi + (-1 + fPsi)*sin(3*beta)) + CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2)))/
    (9.*CG*CK*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    
    
    jacobianf1(0,1) =    -((sqrt(2)*exp(fB*i1)*fB*fC*fPsi*(-((3 + 2*fPsi + 3*pow(fPsi,2))*rhw_sigma[2]*cos(beta)) +
                                                           5*(-1 + pow(fPsi,2))*rhw_sigma[1]*cos(2*beta) + 24*sqrt(2)*fA*fPsi*cos(3*beta) -
                                                           24*sqrt(2)*exp(fB*i1)*fC*fPsi*cos(3*beta) - 24*sqrt(2)*fA*pow(fPsi,2)*cos(3*beta) +
                                                           24*sqrt(2)*exp(fB*i1)*fC*pow(fPsi,2)*cos(3*beta) - rhw_sigma[1]*cos(4*beta) +
                                                           pow(fPsi,2)*rhw_sigma[1]*cos(4*beta) + 2*rhw_sigma[2]*cos(5*beta) - 4*fPsi*rhw_sigma[2]*cos(5*beta) +
                                                           2*pow(fPsi,2)*rhw_sigma[2]*cos(5*beta) - rhw_sigma[2]*cos(7*beta) + 2*fPsi*rhw_sigma[2]*cos(7*beta) -
                                                           pow(fPsi,2)*rhw_sigma[2]*cos(7*beta) + 3*rhw_sigma[1]*sin(beta) + 2*fPsi*rhw_sigma[1]*sin(beta) +
                                                           3*pow(fPsi,2)*rhw_sigma[1]*sin(beta) + 5*rhw_sigma[2]*sin(2*beta) - 5*pow(fPsi,2)*rhw_sigma[2]*sin(2*beta) -
                                                           rhw_sigma[2]*sin(4*beta) + pow(fPsi,2)*rhw_sigma[2]*sin(4*beta) + 2*rhw_sigma[1]*sin(5*beta) -
                                                           4*fPsi*rhw_sigma[1]*sin(5*beta) + 2*pow(fPsi,2)*rhw_sigma[1]*sin(5*beta) + rhw_sigma[1]*sin(7*beta) -
                                                           2*fPsi*rhw_sigma[1]*sin(7*beta) + pow(fPsi,2)*rhw_sigma[1]*sin(7*beta)))/
                           (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),3)));
    
    
    jacobianf1(0,2) =      0;
    
    jacobianf1(1,0) =      -((sqrt(2)*exp(fB*i1)*fB*fC*fPsi*(-((3 + 2*fPsi + 3*pow(fPsi,2))*rhw_sigma[2]*cos(beta)) +
                                                             5*(-1 + pow(fPsi,2))*rhw_sigma[1]*cos(2*beta) + 24*sqrt(2)*fA*fPsi*cos(3*beta) -
                                                             24*sqrt(2)*exp(fB*i1)*fC*fPsi*cos(3*beta) - 24*sqrt(2)*fA*pow(fPsi,2)*cos(3*beta) +
                                                             24*sqrt(2)*exp(fB*i1)*fC*pow(fPsi,2)*cos(3*beta) - rhw_sigma[1]*cos(4*beta) +
                                                             pow(fPsi,2)*rhw_sigma[1]*cos(4*beta) + 2*rhw_sigma[2]*cos(5*beta) - 4*fPsi*rhw_sigma[2]*cos(5*beta) +
                                                             2*pow(fPsi,2)*rhw_sigma[2]*cos(5*beta) - rhw_sigma[2]*cos(7*beta) + 2*fPsi*rhw_sigma[2]*cos(7*beta) -
                                                             pow(fPsi,2)*rhw_sigma[2]*cos(7*beta) + 3*rhw_sigma[1]*sin(beta) + 2*fPsi*rhw_sigma[1]*sin(beta) +
                                                             3*pow(fPsi,2)*rhw_sigma[1]*sin(beta) + 5*rhw_sigma[2]*sin(2*beta) - 5*pow(fPsi,2)*rhw_sigma[2]*sin(2*beta) -
                                                             rhw_sigma[2]*sin(4*beta) + pow(fPsi,2)*rhw_sigma[2]*sin(4*beta) + 2*rhw_sigma[1]*sin(5*beta) -
                                                             4*fPsi*rhw_sigma[1]*sin(5*beta) + 2*pow(fPsi,2)*rhw_sigma[1]*sin(5*beta) + rhw_sigma[1]*sin(7*beta) -
                                                             2*fPsi*rhw_sigma[1]*sin(7*beta) + pow(fPsi,2)*rhw_sigma[1]*sin(7*beta)))/
                             (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),3)));
    
    jacobianf1(1,1) =      ((fA - exp(fB*i1)*fC)*fPsi*(8*(fA - exp(fB*i1)*fC)*fPsi*pow(cos(beta),2)*
                                                       pow(1 + fPsi + 8*(-1 + fPsi)*pow(sin(beta),3),2) +
                                                       2*sqrt(2)*cos(beta)*(15 - 34*fPsi + 15*pow(fPsi,2) - 6*pow(-1 + fPsi,2)*cos(2*beta) +
                                                                            6*pow(-1 + fPsi,2)*cos(4*beta) + 2*cos(6*beta) - 4*fPsi*cos(6*beta) +
                                                                            2*pow(fPsi,2)*cos(6*beta) + 12*sin(beta) - 12*pow(fPsi,2)*sin(beta) - 13*sin(3*beta) +
                                                                            13*pow(fPsi,2)*sin(3*beta))*(2*sqrt(2)*(fA - exp(fB*i1)*fC)*fPsi*cos(beta) -
                                                                                                         rhw_sigma[1]*(1 + fPsi + (-1 + fPsi)*sin(3*beta))) +
                                                       8*(fA - exp(fB*i1)*fC)*fPsi*pow(3*(-1 + fPsi)*cos(beta)*cos(3*beta) +
                                                                                       sin(beta)*(1 + fPsi + (-1 + fPsi)*sin(3*beta)),2) +
                                                       sqrt(2)*(2*sqrt(2)*(fA - exp(fB*i1)*fC)*fPsi*sin(beta) -
                                                                rhw_sigma[2]*(1 + fPsi + (-1 + fPsi)*sin(3*beta)))*
                                                       ((-1 + pow(fPsi,2))*cos(2*beta) - 13*(-1 + pow(fPsi,2))*cos(4*beta) +
                                                        8*(3 - 7*fPsi + 3*pow(fPsi,2))*sin(beta) +
                                                        2*pow(-1 + fPsi,2)*(-4*sin(5*beta) + sin(7*beta)))))/
    (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),4));
    
    jacobianf1(1,2) =    0;
    
    
    jacobianf1(2,0) =     -1;
    
    
    jacobianf1(2,1) =     0;
    
    jacobianf1(2,2) =     -3*CK*(CPer + fD/exp(fD*(CX0 + fA*fR - exp(fB*k)*fC*fR - k)))*(1 + exp(fB*k)*fB*fC*fR)*fW;
    
    
    

    
    
}


void TPZSandlerExtended::Jacobianf2(const TPZVec<STATE> &trial_stress, STATE theta, STATE beta, STATE k, TPZFMatrix<STATE> &jacobianf2)const {

   
    jacobianf2.Resize(3, 3);
    STATE CX0   = X_0();
    STATE CPer  = CPerturbation();
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE CG    = fE/(2.0*(1.0 + fnu));
    
    TPZManVector<REAL,3> rhw_sigma(3);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_sigma);
    

    jacobianf2(0,0) = (2*(fA - exp(fB*k)*fC)*(108*CK*(fA - exp(fB*k)*fC)*pow(fPsi,2)*pow(cos(beta),2)*
                                              pow(cos(theta),2) + 108*CK*(fA - exp(fB*k)*fC)*pow(fPsi,2)*pow(cos(theta),2)*
                                              pow(sin(beta),2) - sqrt(3)*CG*fR*cos(theta)*
                                              (3*rhw_sigma[0] + sqrt(3)*(-k + (fA - exp(fB*k)*fC)*fR*cos(theta)))*
                                              pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2) +
                                              3*CG*(fA - exp(fB*k)*fC)*pow(fR,2)*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2)*
                                              pow(sin(theta),2) + 27*sqrt(2)*CK*fPsi*cos(beta)*sin(theta)*
                                              (rhw_sigma[1]*(1 + fPsi + (-1 + fPsi)*sin(3*beta)) -
                                               2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*cos(beta)*sin(theta)) +
                                              27*sqrt(2)*CK*fPsi*sin(beta)*sin(theta)*
                                              (rhw_sigma[2]*(1 + fPsi + (-1 + fPsi)*sin(3*beta)) -
                                               2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*sin(beta)*sin(theta))))/
    (27.*CG*CK*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    
    
    
    jacobianf2(0,1) =   (sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*cos(theta)*
                         (-((3 + 2*fPsi + 3*pow(fPsi,2))*rhw_sigma[2]*cos(beta)) + 5*(-1 + pow(fPsi,2))*rhw_sigma[1]*cos(2*beta) -
                          rhw_sigma[1]*cos(4*beta) + pow(fPsi,2)*rhw_sigma[1]*cos(4*beta) + 2*rhw_sigma[2]*cos(5*beta) -
                          4*fPsi*rhw_sigma[2]*cos(5*beta) + 2*pow(fPsi,2)*rhw_sigma[2]*cos(5*beta) - rhw_sigma[2]*cos(7*beta) +
                          2*fPsi*rhw_sigma[2]*cos(7*beta) - pow(fPsi,2)*rhw_sigma[2]*cos(7*beta) + 3*rhw_sigma[1]*sin(beta) +
                          2*fPsi*rhw_sigma[1]*sin(beta) + 3*pow(fPsi,2)*rhw_sigma[1]*sin(beta) + 5*rhw_sigma[2]*sin(2*beta) -
                          5*pow(fPsi,2)*rhw_sigma[2]*sin(2*beta) - rhw_sigma[2]*sin(4*beta) + pow(fPsi,2)*rhw_sigma[2]*sin(4*beta) +
                          2*rhw_sigma[1]*sin(5*beta) - 4*fPsi*rhw_sigma[1]*sin(5*beta) + 2*pow(fPsi,2)*rhw_sigma[1]*sin(5*beta) +
                          rhw_sigma[1]*sin(7*beta) - 2*fPsi*rhw_sigma[1]*sin(7*beta) + pow(fPsi,2)*rhw_sigma[1]*sin(7*beta) -
                          12*sqrt(2)*fA*fPsi*sin(3*beta - theta) + 12*sqrt(2)*exp(fB*k)*fC*fPsi*sin(3*beta - theta) +
                          12*sqrt(2)*fA*pow(fPsi,2)*sin(3*beta - theta) -
                          12*sqrt(2)*exp(fB*k)*fC*pow(fPsi,2)*sin(3*beta - theta) +
                          12*sqrt(2)*fA*fPsi*sin(3*beta + theta) - 12*sqrt(2)*exp(fB*k)*fC*fPsi*sin(3*beta + theta) -
                          12*sqrt(2)*fA*pow(fPsi,2)*sin(3*beta + theta) +
                          12*sqrt(2)*exp(fB*k)*fC*pow(fPsi,2)*sin(3*beta + theta)))/
    (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),3));
    
    
    jacobianf2(0,2) =   (2*((9*sqrt(2)*exp(fB*k)*fB*fC*fPsi*rhw_sigma[1]*cos(beta)*cos(theta))/
                            (CG*(1 + fPsi + (-1 + fPsi)*sin(3*beta))) +
                            (9*sqrt(2)*exp(fB*k)*fB*fC*fPsi*rhw_sigma[2]*cos(theta)*sin(beta))/
                            (CG*(1 + fPsi + (-1 + fPsi)*sin(3*beta))) +
                            (sqrt(3)*exp(fB*k)*fB*fC*fR*rhw_sigma[0]*sin(theta))/CK +
                            ((fA - exp(fB*k)*fC)*fR*(1 + exp(fB*k)*fB*fC*fR*cos(theta))*sin(theta))/CK -
                            (exp(fB*k)*fB*fC*fR*(k - (fA - exp(fB*k)*fC)*fR*cos(theta))*sin(theta))/CK +
                            (36*exp(fB*k)*fB*fC*(-fA + exp(fB*k)*fC)*pow(fPsi,2)*pow(cos(beta),2)*sin(2*theta))/
                            (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2)) +
                            (36*exp(fB*k)*fB*fC*(-fA + exp(fB*k)*fC)*pow(fPsi,2)*pow(sin(beta),2)*sin(2*theta))/
                            (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2))))/9.0;
    
    jacobianf2(1,0) = (sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*cos(theta)*
                       (-((3 + 2*fPsi + 3*pow(fPsi,2))*rhw_sigma[2]*cos(beta)) + 5*(-1 + pow(fPsi,2))*rhw_sigma[1]*cos(2*beta) -
                        rhw_sigma[1]*cos(4*beta) + pow(fPsi,2)*rhw_sigma[1]*cos(4*beta) + 2*rhw_sigma[2]*cos(5*beta) -
                        4*fPsi*rhw_sigma[2]*cos(5*beta) + 2*pow(fPsi,2)*rhw_sigma[2]*cos(5*beta) - rhw_sigma[2]*cos(7*beta) +
                        2*fPsi*rhw_sigma[2]*cos(7*beta) - pow(fPsi,2)*rhw_sigma[2]*cos(7*beta) + 3*rhw_sigma[1]*sin(beta) +
                        2*fPsi*rhw_sigma[1]*sin(beta) + 3*pow(fPsi,2)*rhw_sigma[1]*sin(beta) + 5*rhw_sigma[2]*sin(2*beta) -
                        5*pow(fPsi,2)*rhw_sigma[2]*sin(2*beta) - rhw_sigma[2]*sin(4*beta) + pow(fPsi,2)*rhw_sigma[2]*sin(4*beta) +
                        2*rhw_sigma[1]*sin(5*beta) - 4*fPsi*rhw_sigma[1]*sin(5*beta) + 2*pow(fPsi,2)*rhw_sigma[1]*sin(5*beta) +
                        rhw_sigma[1]*sin(7*beta) - 2*fPsi*rhw_sigma[1]*sin(7*beta) + pow(fPsi,2)*rhw_sigma[1]*sin(7*beta) -
                        12*sqrt(2)*fA*fPsi*sin(3*beta - theta) + 12*sqrt(2)*exp(fB*k)*fC*fPsi*sin(3*beta - theta) +
                        12*sqrt(2)*fA*pow(fPsi,2)*sin(3*beta - theta) -
                        12*sqrt(2)*exp(fB*k)*fC*pow(fPsi,2)*sin(3*beta - theta) +
                        12*sqrt(2)*fA*fPsi*sin(3*beta + theta) - 12*sqrt(2)*exp(fB*k)*fC*fPsi*sin(3*beta + theta) -
                        12*sqrt(2)*fA*pow(fPsi,2)*sin(3*beta + theta) +
                        12*sqrt(2)*exp(fB*k)*fC*pow(fPsi,2)*sin(3*beta + theta)))/
    (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),3));
    
    
    jacobianf2(1,1) =   ((fA - exp(fB*k)*fC)*fPsi*sin(theta)*(8*(fA - exp(fB*k)*fC)*fPsi*
                                                              pow(2*(-1 + fPsi)*cos(2*beta) + (-1 + fPsi)*cos(4*beta) + (1 + fPsi)*sin(beta),2)*sin(theta)\
                                                              + 8*(fA - exp(fB*k)*fC)*fPsi*pow(cos(beta),2)*
                                                              pow(1 + fPsi + 8*(-1 + fPsi)*pow(sin(beta),3),2)*sin(theta) -
                                                              2*sqrt(2)*(18*pow(-1 + fPsi,2)*cos(beta)*pow(cos(3*beta),2) +
                                                                         6*(-1 + fPsi)*cos(3*beta)*sin(beta)*(1 + fPsi + (-1 + fPsi)*sin(3*beta)) +
                                                                         9*(-1 + fPsi)*cos(beta)*sin(3*beta)*(1 + fPsi + (-1 + fPsi)*sin(3*beta)) -
                                                                         cos(beta)*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2))*
                                                              (rhw_sigma[1]*(1 + fPsi + (-1 + fPsi)*sin(3*beta)) -
                                                               2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*cos(beta)*sin(theta)) +
                                                              sqrt(2)*((-1 + pow(fPsi,2))*cos(2*beta) - 13*(-1 + pow(fPsi,2))*cos(4*beta) +
                                                                       8*(3 - 7*fPsi + 3*pow(fPsi,2))*sin(beta) +
                                                                       2*pow(-1 + fPsi,2)*(-4*sin(5*beta) + sin(7*beta)))*
                                                              (-((1 + fPsi)*rhw_sigma[2]) - 3*(-1 + fPsi)*rhw_sigma[2]*pow(cos(beta),2)*sin(beta) +
                                                               (-1 + fPsi)*rhw_sigma[2]*pow(sin(beta),3) +
                                                               2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*sin(beta)*sin(theta))))/
    (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),4));
    
    
    jacobianf2(1,2) =   -((sqrt(2)*exp(fB*k)*fB*fC*fPsi*sin(theta)*
                           (-((3 + 2*fPsi + 3*pow(fPsi,2))*rhw_sigma[2]*cos(beta)) +
                            5*(-1 + pow(fPsi,2))*rhw_sigma[1]*cos(2*beta) - rhw_sigma[1]*cos(4*beta) +
                            pow(fPsi,2)*rhw_sigma[1]*cos(4*beta) + 2*rhw_sigma[2]*cos(5*beta) - 4*fPsi*rhw_sigma[2]*cos(5*beta) +
                            2*pow(fPsi,2)*rhw_sigma[2]*cos(5*beta) - rhw_sigma[2]*cos(7*beta) + 2*fPsi*rhw_sigma[2]*cos(7*beta) -
                            pow(fPsi,2)*rhw_sigma[2]*cos(7*beta) + 3*rhw_sigma[1]*sin(beta) + 2*fPsi*rhw_sigma[1]*sin(beta) +
                            3*pow(fPsi,2)*rhw_sigma[1]*sin(beta) + 5*rhw_sigma[2]*sin(2*beta) - 5*pow(fPsi,2)*rhw_sigma[2]*sin(2*beta) -
                            rhw_sigma[2]*sin(4*beta) + pow(fPsi,2)*rhw_sigma[2]*sin(4*beta) + 2*rhw_sigma[1]*sin(5*beta) -
                            4*fPsi*rhw_sigma[1]*sin(5*beta) + 2*pow(fPsi,2)*rhw_sigma[1]*sin(5*beta) + rhw_sigma[1]*sin(7*beta) -
                            2*fPsi*rhw_sigma[1]*sin(7*beta) + pow(fPsi,2)*rhw_sigma[1]*sin(7*beta) -
                            12*sqrt(2)*fA*fPsi*sin(3*beta - theta) +
                            12*sqrt(2)*exp(fB*k)*fC*fPsi*sin(3*beta - theta) +
                            12*sqrt(2)*fA*pow(fPsi,2)*sin(3*beta - theta) -
                            12*sqrt(2)*exp(fB*k)*fC*pow(fPsi,2)*sin(3*beta - theta) +
                            12*sqrt(2)*fA*fPsi*sin(3*beta + theta) -
                            12*sqrt(2)*exp(fB*k)*fC*fPsi*sin(3*beta + theta) -
                            12*sqrt(2)*fA*pow(fPsi,2)*sin(3*beta + theta) +
                            12*sqrt(2)*exp(fB*k)*fC*pow(fPsi,2)*sin(3*beta + theta)))/
                          (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),3)));
    
    jacobianf2(2,0) =     -((fA - exp(fB*k)*fC)*fR*sin(theta));
    
    jacobianf2(2,1) =     0;
    
    
    jacobianf2(2,2) =     -1 - 3*CK*(CPer + fD/exp(fD*(CX0 + fA*fR - exp(fB*k)*fC*fR - k)))*(1 + exp(fB*k)*fB*fC*fR)*
    fW - exp(fB*k)*fB*fC*fR*cos(theta);
    
    
}

/// Compute the jacobian function of the f2 (cap) distance as a function of beta and k
void TPZSandlerExtended::Jacobianf2CoVertex(const TPZVec<STATE> &trial_stress, STATE beta, STATE k, TPZFMatrix<STATE> &jacobianf2_covertex)const{
    
    jacobianf2_covertex.Resize(2, 2);
    STATE CX0   = X_0();
    STATE CPer  = CPerturbation();
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE CG    = fE/(2.0*(1.0 + fnu));
    STATE theta = M_PI_2;
    
    TPZManVector<REAL,3> rhw_sigma(3);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_sigma);
    
    jacobianf2_covertex(0,0) = pow((2*sqrt(2)*(fA - exp(fB*k)*fC)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)*sin(beta))/
        pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) -
        (2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)),2)/CG +
    pow((2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi))/
        pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) +
        (2*sqrt(2)*(fA - exp(fB*k)*fC)*sin(beta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)),2)/CG +
    ((rhw_sigma[1] - (2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)))*
     ((-4*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*pow(3*cos(3*beta) - (3*cos(3*beta))/fPsi,2))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),3) -
      (4*sqrt(2)*(fA - exp(fB*k)*fC)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)*sin(beta))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) +
      (2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)) +
      (2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*(-9*sin(3*beta) + (9*sin(3*beta))/fPsi))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2)))/CG +
    ((rhw_sigma[2] - (2*sqrt(2)*(fA - exp(fB*k)*fC)*sin(beta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)))*
     ((-4*sqrt(2)*(fA - exp(fB*k)*fC)*pow(3*cos(3*beta) - (3*cos(3*beta))/fPsi,2)*sin(beta))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),3) +
      (4*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) +
      (2*sqrt(2)*(fA - exp(fB*k)*fC)*sin(beta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)) +
      (2*sqrt(2)*(fA - exp(fB*k)*fC)*sin(beta)*(-9*sin(3*beta) + (9*sin(3*beta))/fPsi))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2)))/CG;
    
    
    jacobianf2_covertex(0,1) = (2*sqrt(2)*exp(fB*k)*fB*fC*sin(beta)*((2*sqrt(2)*(fA - exp(fB*k)*fC)*
                                           (3*cos(3*beta) - (3*cos(3*beta))/fPsi)*sin(beta))/
                                          pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) -
                                          (2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta))))/
    (CG*(1 + (1 - sin(3*beta))/fPsi + sin(3*beta))) +
    ((rhw_sigma[1] - (2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)))*
     ((-2*sqrt(2)*exp(fB*k)*fB*fC*cos(beta)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) -
      (2*sqrt(2)*exp(fB*k)*fB*fC*sin(beta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta))))/CG +
    (((-2*sqrt(2)*exp(fB*k)*fB*fC*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)*sin(beta))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) +
      (2*sqrt(2)*exp(fB*k)*fB*fC*cos(beta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)))*
     (rhw_sigma[2] - (2*sqrt(2)*(fA - exp(fB*k)*fC)*sin(beta))/
      (1 + (1 - sin(3*beta))/fPsi + sin(3*beta))))/CG +
    (2*sqrt(2)*exp(fB*k)*fB*fC*cos(beta)*
     ((2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) +
      (2*sqrt(2)*(fA - exp(fB*k)*fC)*sin(beta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta))))/
    (CG*(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)));
    
    jacobianf2_covertex(1,0) = 0;
    
    jacobianf2_covertex(1,1) = -1 - 3*CK*(CPer*(1 + exp(fB*k)*fB*fC*fR) +
               exp(fD*(-CX0 - (fA - exp(fB*k)*fC)*fR + k))*fD*(1 + exp(fB*k)*fB*fC*fR))*fW;
    
}

/// Compute the jacobian function of the vertex on f2 (cap) distance as a function of k
void TPZSandlerExtended::Jacobianf2Vertex(const TPZVec<STATE> &trial_stress, STATE k, STATE &jacobianf2_vertex)const{

    STATE CX0   = X_0();
    STATE CPer  = CPerturbation();
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE theta = 0.0;
    
//    TPZManVector<REAL,3> rhw_sigma(3);
//    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_sigma);
    
    jacobianf2_vertex = -1 - 3*CK*(CPer + fD/exp(fD*(CX0 + fA*fR - exp(fB*k)*fC*fR - k)))*(1 + exp(fB*k)*fB*fC*fR)*fW -
    exp(fB*k)*fB*fC*fR*cos(theta);
}

void TPZSandlerExtended::JacobianVertex(const TPZVec<STATE> &trial_stress, STATE k, STATE &jacobian_vertex) const {
    
    STATE theta = 0.0;
    STATE CX0   = X_0();
    STATE CPer  = CPerturbation();
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE CG    = fE/(2.0*(1.0 + fnu));
    
    jacobian_vertex = -1 - 3*CK*(CPer*(1 + exp(fB*k)*fB*fC*fR) +
               exp(fD*(-CX0 - (fA - exp(fB*k)*fC)*fR + k))*fD*(1 + exp(fB*k)*fB*fC*fR))*fW - exp(fB*k)*fB*fC*fR*cos(theta);
    
}

void TPZSandlerExtended::JacobianCoVertex(const TPZVec<STATE> &trial_stress, STATE beta, STATE k, TPZFMatrix<STATE> &jacobian_covertex) const {
    
    // In this implementation the definition for theta is given by the angle formed from -I1 axis to sqrt(J2) axis with origin on the damage variable kappa.
    // Thus, covertex is located at theta  = M_PI_2;
    STATE theta = M_PI_2;
    
    jacobian_covertex.Resize(2, 2);
    STATE CX0   = X_0();
    STATE CPer  = CPerturbation();
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE CG    = fE/(2.0*(1.0 + fnu));
    
    TPZManVector<REAL,3> rhw_sigma(3);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_sigma);

    jacobian_covertex(0,0) = pow((2*sqrt(2)*(fA - exp(fB*k)*fC)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)*sin(beta)*sin(theta))/
        pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) -
        (2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*sin(theta))/
        (1 + (1 - sin(3*beta))/fPsi + sin(3*beta)),2)/CG +
    pow((2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)*
         sin(theta))/pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) +
        (2*sqrt(2)*(fA - exp(fB*k)*fC)*sin(beta)*sin(theta))/
        (1 + (1 - sin(3*beta))/fPsi + sin(3*beta)),2)/CG +
    ((rhw_sigma[1] - (2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*sin(theta))/
      (1 + (1 - sin(3*beta))/fPsi + sin(3*beta)))*
     ((-4*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*pow(3*cos(3*beta) - (3*cos(3*beta))/fPsi,2)*
       sin(theta))/pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),3) -
      (4*sqrt(2)*(fA - exp(fB*k)*fC)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)*sin(beta)*sin(theta))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) +
      (2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*sin(theta))/
      (1 + (1 - sin(3*beta))/fPsi + sin(3*beta)) +
      (2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*(-9*sin(3*beta) + (9*sin(3*beta))/fPsi)*
       sin(theta))/pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2)))/CG +
    ((rhw_sigma[2] - (2*sqrt(2)*(fA - exp(fB*k)*fC)*sin(beta)*sin(theta))/
      (1 + (1 - sin(3*beta))/fPsi + sin(3*beta)))*
     ((-4*sqrt(2)*(fA - exp(fB*k)*fC)*pow(3*cos(3*beta) - (3*cos(3*beta))/fPsi,2)*sin(beta)*
       sin(theta))/pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),3) +
      (4*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)*sin(theta))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) +
      (2*sqrt(2)*(fA - exp(fB*k)*fC)*sin(beta)*sin(theta))/
      (1 + (1 - sin(3*beta))/fPsi + sin(3*beta)) +
      (2*sqrt(2)*(fA - exp(fB*k)*fC)*sin(beta)*(-9*sin(3*beta) + (9*sin(3*beta))/fPsi)*
       sin(theta))/pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2)))/CG;
    
    
    jacobian_covertex(0,1) = (2*sqrt(2)*exp(fB*k)*fB*fC*sin(beta)*sin(theta)*
     ((2*sqrt(2)*(fA - exp(fB*k)*fC)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)*sin(beta)*sin(theta))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) -
      (2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*sin(theta))/
      (1 + (1 - sin(3*beta))/fPsi + sin(3*beta))))/(CG*(1 + (1 - sin(3*beta))/fPsi + sin(3*beta))) +
    ((rhw_sigma[1] - (2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*sin(theta))/
      (1 + (1 - sin(3*beta))/fPsi + sin(3*beta)))*
     ((-2*sqrt(2)*exp(fB*k)*fB*fC*cos(beta)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)*sin(theta))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) -
      (2*sqrt(2)*exp(fB*k)*fB*fC*sin(beta)*sin(theta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)))
     )/CG + (((-2*sqrt(2)*exp(fB*k)*fB*fC*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)*sin(beta)*
               sin(theta))/pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) +
              (2*sqrt(2)*exp(fB*k)*fB*fC*cos(beta)*sin(theta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)))
             *(rhw_sigma[2] - (2*sqrt(2)*(fA - exp(fB*k)*fC)*sin(beta)*sin(theta))/
               (1 + (1 - sin(3*beta))/fPsi + sin(3*beta))))/CG +
    (2*sqrt(2)*exp(fB*k)*fB*fC*cos(beta)*sin(theta)*
     ((2*sqrt(2)*(fA - exp(fB*k)*fC)*cos(beta)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)*sin(theta))/
      pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) +
      (2*sqrt(2)*(fA - exp(fB*k)*fC)*sin(beta)*sin(theta))/
      (1 + (1 - sin(3*beta))/fPsi + sin(3*beta))))/(CG*(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)));
    
    
    jacobian_covertex(1,0) = 0;
    
    jacobian_covertex(1,1) = -1 - 3*CK*(CPer*(1 + exp(fB*k)*fB*fC*fR) +
               exp(fD*(-CX0 - (fA - exp(fB*k)*fC)*fR + k))*fD*(1 + exp(fB*k)*fB*fC*fR))*fW -
    exp(fB*k)*fB*fC*fR*cos(theta);
    
}

void TPZSandlerExtended::YieldFunction(const TPZVec<STATE> &sigma, STATE kprev, TPZVec<STATE> &yield) const {

    yield.resize(3);
    STATE II1, JJ2, ggamma, temp1, temp3, f1, f2, phi, sqrtj2, X, xi, rho, beta;
    TPZManVector<STATE, 3> cylstress(3);
    TPZHWTools::FromPrincipalToHWCyl(sigma, cylstress);
    
    // Zylinderkoordinaten
    xi = cylstress[0];
    rho = cylstress[1];
    beta = cylstress[2];

    II1 = sqrt(3.0)*xi;
    JJ2 = 0.5*rho*rho;

    if (IsZero(JJ2)) {
        JJ2 = 0.0;
    }
    
    sqrtj2 = sqrt(JJ2);
    ggamma = 0.5 * (1. + (1. - sin(3. * beta)) / fPsi + sin(3. * beta));

    temp1 = (kprev-II1) / (fR * F(kprev));
    temp3 = (ggamma * sqrtj2) / (F(kprev));

    f1 = sqrtj2 - F(II1)/ggamma;
    f2 = (temp1 * temp1 + temp3 * temp3 - 1);
    
    X = kprev - fR * F(kprev);
    
    // hardcoded
    if (II1 > kprev) {
        phi = f1;
    }else if (II1 > X || IsZero(II1-X) ) {
        phi = f2;
    }else{
        phi = 0.0;
    }
    
    yield[0] = f1;
    yield[1] = f2;
    yield[2] = phi;

}

void TPZSandlerExtended::Print(std::ostream& out) const {
    out << "TPZSandlerExtended\n";
    out << "A: " << fA << std::endl;
    out << "B: " << fB << std::endl;
    out << "C: " << fC << std::endl;
    out << "D: " << fD << std::endl;
    out << "R: " << fR << std::endl;
    out << "W: " << fW << std::endl;
    out << "k_0: " << fkappa_0 << std::endl;
}


void TPZSandlerExtended::Phi(TPZVec<REAL> sigma, STATE alpha, TPZVec<STATE> &phi)const {

    //    TPZTensor<REAL>::TPZDecomposed DecompSig;
    //    TPZTensor<STATE> sig;
    //    TPZElasticResponse ER;
    //    ER.Compute(eps,sig);
    //    sig.EigenSystem(DecompSig);
    YieldFunction(sigma, alpha, phi);
}

std::map<int, int64_t> gF1Stat;
std::map<int, int64_t> gF2Stat;
std::vector<int64_t> gYield;

STATE TPZSandlerExtended::NormalFunctionToF1(STATE & I1, STATE & k) const{
    
    STATE normal_to_failure = I1/(exp(fB*k)*fB*fC) + (sqrt(1 + exp(2*fB*k)*pow(fB,2)*pow(fC,2))*
                            (-((fA - exp(fB*k)*fC + 1/sqrt(1 + exp(2*fB*k)*pow(fB,2)*pow(fC,2)))*k) +
                             (fA - exp(fB*k)*fC)*((exp(fB*k)*fB*fC)/
                                                  sqrt(1 + exp(2*fB*k)*pow(fB,2)*pow(fC,2)) + k)))/(exp(fB*k)*fB*fC);
    return normal_to_failure;
    
}

void TPZSandlerExtended::ProjectApex(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigproj, STATE &kproj) const {
    
    REAL K = fElasticResponse.K();
    REAL xi_apex = Apex()/3.0;
    
    // Trial
    REAL ptr_np1 = 0.;
    for (int i = 0; i < 3; i++) {
        ptr_np1 += sigmatrial[i];
    }
    ptr_np1 /= 3.;
    
    REAL p_np1;
    REAL delta_eps_np1 = 0;
    REAL res = xi_apex - ptr_np1;
    REAL jac;
    
    bool stop_criterion;
    int n_iterations = 50; // @TODO : Define a numeric controls manager object and use it to obtain this information
    int i;
    for (i = 0; i < n_iterations; i++) {
        jac = K;
        delta_eps_np1 -= res / jac;
        p_np1 = ptr_np1 - K * delta_eps_np1;
        res = xi_apex - p_np1;
        stop_criterion = IsZero(res);
        if (stop_criterion) {
            break;
        }
    }
    
#ifdef PZDEBUG
    if (i == n_iterations) {
#ifdef WIN32
	REAL tol = 1.e-10;
#else
	REAL tol = 1.e-12;
#endif
        throw TPZConvergenceException(tol, n_iterations, res, i, "TPZSandlerExtended::ProjectApex:: Newton process did not converge.");
    }
#endif
    
    for (int i = 0; i < 3; i++) {
        sigproj[i] = p_np1;
    }
    kproj = kprev;
}

void TPZSandlerExtended::ProjectF1(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> & projected_stress, STATE &kproj) const {
#ifdef LOG4CXX
    if (loggerConvTest->isDebugEnabled()) {
        std::stringstream outfile;
        outfile << "\n projection over F1 " << endl;
        LOGPZ_DEBUG(loggerConvTest, outfile.str());
    }
#endif

    TPZManVector<STATE, 3> hw_space_xi_rho_beta(3);
    TPZManVector<STATE, 3> rhw_space_s1_s2_s3(3);;
    TPZHWTools::FromPrincipalToHWCyl(trial_stress, hw_space_xi_rho_beta);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_space_s1_s2_s3);
    
    STATE i1_guess = sqrt(3.0)*rhw_space_s1_s2_s3[0];
    STATE beta_guess = atan2(rhw_space_s1_s2_s3[2],rhw_space_s1_s2_s3[1]);
    STATE k_guess = sqrt(3.0)*rhw_space_s1_s2_s3[0];
    
    TPZFNMatrix<3, STATE> delta_par(3, 1, 0.), par(3, 1, 0.), residue(3, 1, 0.);
    par(0, 0) = i1_guess;
    par(1, 0) = beta_guess;
    par(2, 0) = k_guess;
    
    TPZFNMatrix<9, STATE> jac(3, 3);
    TPZFNMatrix<9, STATE> jac_inv(3,3);
    
    TPZManVector<STATE,3> residue_vec(3);
    STATE residue_norm;
    bool stop_criterion_res;
    int max_iterations = 50;
    int it;
    for (it = 0; it < max_iterations; it++) {
        // Computing the Residue vector for a Newton step
        Res1(trial_stress, par(0), par(1), par(2), kprev, residue_vec); // Residue
        for (int k = 0; k < 3; k++) residue(k, 0) = - 1.0 * residue_vec[k]; // Transfering to a Matrix object
        residue_norm = Norm(residue);
        stop_criterion_res = std::fabs(residue_norm) <= ftol;
        if (stop_criterion_res) {
            break;
        }
        
        Jacobianf1(trial_stress, par(0), par(1), par(2), jac); // Jacobian
        TPZHWTools::A3x3Inverse(jac, jac_inv);
        jac_inv.Multiply(residue, delta_par);
        
        par += delta_par;
    }
    
#ifdef PZDEBUG
    if (it == max_iterations) {
        throw TPZConvergenceException(ftol, max_iterations, residue_norm, it, "TPZSandlerExtended::ProjectF1:: Newton process did not converge.");
    }
#endif
    
    STATE xi, beta;
    xi      = par(0)/sqrt(3.0);
    beta    = par(1);
    kproj   = par(2);
    
    TPZManVector<STATE, 3> f1cyl(3);
    F1Cyl(xi, beta, f1cyl);
    TPZHWTools::FromHWCylToPrincipal(f1cyl, projected_stress);
    
}

void TPZSandlerExtended::ProjectF2(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const {

    TPZManVector<STATE, 3> hw_space_xi_rho_beta(3);
    TPZManVector<STATE, 3> rhw_space_s1_s2_s3(3);;
    TPZHWTools::FromPrincipalToHWCyl(trial_stress, hw_space_xi_rho_beta);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_space_s1_s2_s3);
    
    STATE k_guess = kprev;
    STATE beta_guess = atan2(rhw_space_s1_s2_s3[2],rhw_space_s1_s2_s3[1]);
    STATE theta_guess = atan2(hw_space_xi_rho_beta[1],(k_guess/sqrt(3))-sqrt(3)*hw_space_xi_rho_beta[0]);

    bool cap_vertex_validity_Q = theta_guess < ftol;
    if (cap_vertex_validity_Q) {
//        std::cout << "Projecting on Cap Vertex " << std::endl;
        ProjectCapVertex(trial_stress, kprev, projected_stress, kproj);
        return;
    }
    
    if (theta_guess > M_PI_2) { // Restriction for theta
//        std::cout << "Reached restriction for theta guess = " << theta_guess*(180.0/M_PI) <<  std::endl;
//        std::cout << "Reached restriction for beta  guess = " << beta_guess*(180.0/M_PI) <<  std::endl;
        theta_guess = M_PI_2;
    }
    
    TPZFNMatrix<3, STATE> delta_par(3, 1, 0.), par(3, 1, 0.), residue(3, 1, 0.);
    par(0, 0) = theta_guess;
    par(1, 0) = beta_guess;
    par(2, 0) = k_guess;

    TPZFNMatrix<9, STATE> jac(3, 3);
    TPZFNMatrix<9, STATE> jac_inv(3,3);
    
    TPZManVector<STATE,3> residue_vec(3);
    STATE residue_norm;
    bool stop_criterion_res;
    int max_iterations = 50;
    int it;
    for (it = 0; it < max_iterations; it++) {
        // Computing the Residue vector for a Newton step
        Res2(trial_stress, par(0), par(1), par(2), kprev, residue_vec); // Residue
        for (int k = 0; k < 3; k++) residue(k, 0) = - 1.0 * residue_vec[k]; // Transfering to a Matrix object
        
        residue_norm = Norm(residue);
        stop_criterion_res = std::fabs(residue_norm) <= ftol;
        if (stop_criterion_res) {
            break;
        }
        
        Jacobianf2(trial_stress, par(0), par(1), par(2), jac); // Jacobian
        TPZHWTools::A3x3Inverse(jac, jac_inv);
        jac_inv.Multiply(residue, delta_par);
        par += delta_par;
        
        if (par(0) > M_PI_2) { // Restriction for theta
//            std::cout << "Reached restriction for theta = " << par(0)*(180.0/M_PI) <<  std::endl;
//            std::cout << "Reached restriction for beta = " << par(1)*(180.0/M_PI) <<  std::endl;
//            std::cout << "Reached restriction for k = " << par(2) <<  std::endl;
            par(0) = M_PI_2;
        }
    }
#ifdef PZDEBUG
    if (it == max_iterations) {
        throw TPZConvergenceException(ftol, max_iterations, residue_norm, it, "TPZSandlerExtended::ProjectF2:: Newton process did not converge.");
    }
#endif

    STATE theta, beta;
    theta   = par(0);
    beta    = par(1);
    kproj   = par(2);

    TPZManVector<STATE, 3> f2cyl(3);
    F2Cyl(theta, beta, kproj, f2cyl);
    TPZHWTools::FromHWCylToPrincipal(f2cyl, projected_stress);
}

void TPZSandlerExtended::ProjectCoVertex(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const{

    TPZManVector<STATE, 3> hw_space_xi_rho_beta(3);
    TPZManVector<STATE, 3> rhw_space_s1_s2_s3(3);;
    TPZHWTools::FromPrincipalToHWCyl(trial_stress, hw_space_xi_rho_beta);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_space_s1_s2_s3);
    
    STATE k_guess = kprev;
    STATE beta_guess = atan2(rhw_space_s1_s2_s3[2],rhw_space_s1_s2_s3[1]);
    STATE theta_guess = M_PI_2;
    
    TPZFNMatrix<2, STATE> delta_par(2, 1, 0.), par(2, 1, 0.), residue(2, 1, 0.);
    par(0, 0) = beta_guess;
    par(1, 0) = k_guess;
    
    TPZFNMatrix<9, STATE> jac(2, 2);
    TPZFNMatrix<9, STATE> jac_inv(2,2);
    
    TPZManVector<STATE,3> residue_vec(2);
    STATE residue_norm;
    bool stop_criterion_res;
    int max_iterations = 50;
    int it;
    for (it = 0; it < max_iterations; it++) {
        
        // Computing the Residue vector for a Newton step
        Res2CoVertex(trial_stress, par(0), par(1), kprev, residue_vec); // Residue
        for (int k = 0; k < 2; k++) residue(k, 0) = - 1.0 * residue_vec[k]; // Transfering to a Matrix object
        
        residue_norm = Norm(residue);
        stop_criterion_res = std::fabs(residue_norm) <= ftol;
        if (stop_criterion_res) {
            break;
        }
        
        JacobianCoVertex(trial_stress, par(0), par(1), jac); // Jacobian
        TPZHWTools::A2x2Inverse(jac, jac_inv);
        jac_inv.Multiply(residue, delta_par);
        par += delta_par;
    }
    
#ifdef PZDEBUG
    if (it == max_iterations) {
        throw TPZConvergenceException(ftol, max_iterations, residue_norm, it, "TPZSandlerExtended::ProjectCoVertex:: Newton process did not converge.");
    }
#endif
    
    STATE theta, beta;
    theta   = M_PI_2;
    beta    = par(0);
    kproj   = par(1);
    
    TPZManVector<STATE, 3> f2cyl(3);
    F2Cyl(theta, beta, kproj, f2cyl);
    TPZHWTools::FromHWCylToPrincipal(f2cyl, projected_stress);
    
}

void TPZSandlerExtended::ProjectVertex(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const{
    
    TPZManVector<STATE, 3> hw_space_xi_rho_beta(3);
    TPZManVector<STATE, 3> rhw_space_s1_s2_s3(3);;
    TPZHWTools::FromPrincipalToHWCyl(trial_stress, hw_space_xi_rho_beta);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_space_s1_s2_s3);
    
    STATE k = kprev;
    STATE dk;

    STATE jac, jac_inv;
    STATE residue;
    bool stop_criterion_res;
    int max_iterations = 50;
    int it;
    for (it = 0; it < max_iterations; it++) {
        // Computing the Residue vector for a Newton step
        Res2Vertex(trial_stress, k, kprev, residue); // Residue
        stop_criterion_res = std::fabs(residue) <= ftol;
        if (stop_criterion_res) {
            break;
        }
        
        JacobianVertex(trial_stress, k, jac); // Jacobian
        jac_inv = 1.0/jac;
        dk = jac_inv * residue;
        k += dk;
    }
#ifdef PZDEBUG
    if (it == max_iterations) {
        throw TPZConvergenceException(ftol, max_iterations, stop_criterion_res, it, "TPZSandlerExtended::ProjectVertex:: Newton process did not converge.");
    }
#endif
    
    TPZManVector<STATE, 3> f2cyl(3);
    f2cyl[0] = (k - fR * F(k))/sqrt(3.0); // xi
    f2cyl[1] = 0.0; // rho
    f2cyl[2] = 0.0; // beta
    TPZHWTools::FromHWCylToPrincipal(f2cyl, projected_stress);
    
}

void TPZSandlerExtended::ProjectCapVertex(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const{
    
    TPZManVector<STATE, 3> hw_space_xi_rho_beta(3);
    TPZManVector<STATE, 3> rhw_space_s1_s2_s3(3);;
    TPZHWTools::FromPrincipalToHWCyl(trial_stress, hw_space_xi_rho_beta);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_space_s1_s2_s3);
    
    STATE k = kprev, dk;
    
    TPZManVector<STATE,3> residue_vec(3);
    STATE residue, jac, jac_inv;
    bool stop_criterion_res;
    int max_iterations = 50;
    int it;
    for (it = 0; it < max_iterations; it++) {
        // Computing the Residue vector for a Newton step
        Res2Vertex(trial_stress, k, kprev, residue); // Residue
        
        stop_criterion_res = std::fabs(residue) <= ftol;
        if (stop_criterion_res) {
            break;
        }
        
        Jacobianf2Vertex(trial_stress, k, jac); // Jacobian
        jac_inv = 1.0 / jac;
        dk = - residue * jac_inv;
        k += dk;
    }
#ifdef PZDEBUG
    if (it == max_iterations) {        
        throw TPZConvergenceException(ftol, max_iterations, residue, it, "TPZSandlerExtended::ProjectCapVertex:: Newton process did not converge.");
    }
#endif
    
    
    STATE theta, beta;
    theta   = 0.0;
    beta    = atan2(0.0, 0.0);
    kproj   = k;

    TPZManVector<STATE, 3> f2cyl(3);
    F2Cyl(theta, beta, kproj, f2cyl); // The definition for theta was corrected.
    TPZHWTools::FromHWCylToPrincipal(f2cyl, projected_stress);
    
}

void TPZSandlerExtended::ProjectCapCoVertex(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const{
    
    TPZManVector<STATE, 3> hw_space_xi_rho_beta(3);
    TPZManVector<STATE, 3> rhw_space_s1_s2_s3(3);;
    TPZHWTools::FromPrincipalToHWCyl(trial_stress, hw_space_xi_rho_beta);
    TPZHWTools::FromPrincipalToHWCart(trial_stress, rhw_space_s1_s2_s3);
    
    STATE k_guess = kprev;
    STATE beta_guess = atan2(rhw_space_s1_s2_s3[2],rhw_space_s1_s2_s3[1]);
    
    TPZFNMatrix<2, STATE> delta_par(2, 1, 0.), par(2, 1, 0.), residue(2, 1, 0.);
    par(0, 0) = beta_guess;
    par(1, 0) = k_guess;
    
    TPZFNMatrix<9, STATE> jac(2, 2);
    TPZFNMatrix<9, STATE> jac_inv(2,2);
    
    TPZManVector<STATE,2> residue_vec(2);
    STATE residue_norm;
    bool stop_criterion_res;
    int max_terations = 50;
    int it;
    for (it = 0; it < max_terations; it++) {
        // Computing the Residue vector for a Newton step
        Res2CoVertex(trial_stress, par(0), par(1), kprev, residue_vec); // Residue
        for (int k = 0; k < 2; k++) residue(k, 0) = - 1.0 * residue_vec[k]; // Transfering to a Matrix object
        
        residue_norm = Norm(residue);
        stop_criterion_res = std::fabs(residue_norm) <= ftol;
        if (stop_criterion_res) {
            break;
        }
        
        Jacobianf2CoVertex(trial_stress, par(0), par(1), jac); // Jacobian
        TPZHWTools::A2x2Inverse(jac, jac_inv);
        jac_inv.Multiply(residue, delta_par);
        par += delta_par;
    }
#ifdef PZDEBUG
    if (it == max_terations) {
        throw TPZConvergenceException(ftol, max_terations, residue_norm, it, "TPZSandlerExtended::ProjectCapCoVertex:: Newton process did not converge.");
    }
#endif
    
    STATE theta, beta;
    theta   = M_PI_2;
    beta    = par(0);
    kproj   = par(1);
    
    TPZManVector<STATE, 3> f2cyl(3);
    F2Cyl(theta, beta, kproj, f2cyl); // The definition for theta was corrected.
    TPZHWTools::FromHWCylToPrincipal(f2cyl, projected_stress);
    
}

void TPZSandlerExtended::ProjectRing(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigproj, STATE &kproj) const {
#ifdef LOG4CXX
	if (loggerConvTest->isDebugEnabled()) {
		std::stringstream outfile;
		outfile << "\n projection over Ring " << endl;
		LOGPZ_DEBUG(loggerConvTest, outfile.str());
	}
#endif
	STATE beta = 0., distnew;
    STATE resnorm, disttheta;
    disttheta = 1.e8;

    for (STATE betaguess = 0; betaguess <= 2 * M_PI; betaguess += M_PI / 20.) {
        distnew = DistF2(sigmatrial, M_PI / 2, betaguess, kprev);
        if (fabs(distnew) < fabs(disttheta)) {
            //            STATE theta = M_PI/2;
            beta = betaguess;
            disttheta = distnew;
        }
    }
    resnorm = 1;
    int64_t counter = 1;
    TPZFMatrix<STATE> xn1(3, 1, 0.), xn(3, 1, 0.), fxn(3, 1, 0.), diff(3, 1, 0.);
    TPZFNMatrix<3, STATE> sol(3, 1, 0.);
    xn(0, 0) = M_PI / 2;
    xn(1, 0) = beta;
    xn(2, 0) = kprev;
    while (resnorm > ftol && counter < 30) {
        TPZFNMatrix<9, STATE> jac(3, 3);
        Jacobianf2(sigmatrial, xn[0], xn[1], xn[2], jac);
        TPZManVector<STATE> fxnvec(3);
        Res2(sigmatrial, xn(0), xn(1), xn(2), kprev, fxnvec);
        for (int k = 0; k < 3; k++) fxn(k, 0) = fxnvec[k];

        for (int i = 0; i < 3; i++) {
            jac(i, 0) = 0.;
            jac(0, i) = 0.;
        }
        jac(0, 0) = 1.;
        fxn(0, 0) = 0.;
        sol = fxn;
        resnorm = Norm(sol);
        jac.Solve_LU(&sol);

        xn1(0, 0) = xn(0, 0);
        xn1(1, 0) = xn(1, 0) - sol(1, 0);
        xn1(2, 0) = xn(2, 0) - sol(2, 0);

#ifdef LOG4CXX
		if (loggerConvTest->isDebugEnabled()) {
			std::stringstream outfile; //("convergencF1.txt");
			outfile << counter << " " << log(resnorm) << endl;
			//jac.Print(outfile);
			//outfile<< "\n xn " << " "<<fxnvec <<endl;
			//outfile<< "\n res " << " "<<fxnvec <<endl;
			LOGPZ_DEBUG(loggerConvTest, outfile.str());
		}
#endif
        //diff=xn1-xn;
        //resnorm=Norm(diff);
        xn = xn1;
        counter++;

    }
#ifdef PZDEBUG
	if (counter == 30) {
            throw TPZConvergenceException(ftol, 30, resnorm, counter, "TPZSandlerExtended::ProjectRing:: Newton process did not converge.");
	}
#endif
    //    cout<< "\n resnorm = "<<resnorm <<endl;
    //    cout<< "\n counter = "<<xn1 <<endl;
    //    cout<< "\n k = "<<xn1[2] <<endl;
    STATE thetasol, betasol, ksol;

    thetasol = xn1[0];
    betasol = xn1[1];
    ksol = xn1[2];

    TPZManVector<STATE, 3> f2cyl(3);
    F2Cyl(thetasol, betasol, ksol, f2cyl);
    TPZHWTools::FromHWCylToPrincipal(f2cyl, sigproj);

    kproj = ksol;

}

void TPZSandlerExtended::ProjectBetaConstF2(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigproj, STATE &kproj) const {
    //#ifdef LOG4CXX
    //    {
    //        std::stringstream outfile;
    //        outfile << "\n projection over F2 " <<endl;
    //        LOGPZ_DEBUG(logger,outfile.str());
    //
    //    }
    //#endif

    STATE theta = 0., beta = 0., distnew;
    STATE resnorm, disttheta;
    disttheta = 1.e8;
    STATE betaconst = 0;
    for (STATE thetaguess = 0; thetaguess <= M_PI; thetaguess += M_PI / 20.) {
        distnew = DistF2(sigmatrial, thetaguess, betaconst, kprev);
        if (fabs(distnew) < fabs(disttheta)) {
            theta = thetaguess;
            beta = betaconst;
            disttheta = distnew;
        }
    }

    resnorm = 1;
    int counter = 1;
    TPZFNMatrix<3, STATE> xn1(3, 1, 0.), xn(3, 1, 0.), sol(3, 1, 0.), fxn(3, 1, 0.), diff(3, 1, 0.);
    xn(0, 0) = theta;
    xn(1, 0) = beta;
    xn(2, 0) = kprev;
    while (resnorm > ftol && counter < 30) {
        TPZFNMatrix<9, STATE> jac(3, 3);
        Jacobianf2(sigmatrial, xn(0), xn(1), xn(2), jac);
        TPZManVector<STATE> fxnvec(3);
        Res2(sigmatrial, xn(0), xn(1), xn(2), kprev, fxnvec);
        for (int k = 0; k < 3; k++) fxn(k, 0) = fxnvec[k];
        for (int i = 0; i < 3; i++) {
            jac(i, 1) = 0.;
            jac(1, i) = 0.;
        }
        jac(1, 1) = 1.;
        fxn(1, 0) = 0.;
        sol = fxn;
        resnorm = Norm(sol);
        jac.Solve_LU(&sol);

        xn1(0, 0) = xn(0, 0) - sol(0, 0);
        xn1(1, 0) = xn(1, 0);
        xn1(2, 0) = xn(2, 0) - sol(2, 0);

        //        diff=xn1-xn;
        //        resnorm=Norm(diff);
        xn = xn1;
        counter++;

    }
    //if(counter == 30) cout << "resnorm = " << resnorm << std::endl;
    STATE thetasol, betasol, ksol;




    thetasol = xn1(0);
    betasol = xn1(1);
    ksol = xn1(2);
    kproj = ksol;

    TPZManVector<STATE, 3> f2cyl(3);
    F2Cyl(thetasol, betasol, ksol, f2cyl);
    TPZHWTools::FromHWCylToPrincipal(f2cyl, sigproj);
}

void TPZSandlerExtended::ComputeI1(TPZVec<STATE> stress, STATE &I1)const {
    STATE sig1, sig2, sig3;
    sig1 = stress[0];
    sig2 = stress[1];
    sig3 = stress[2];
    I1 = sig1 + sig2 + sig3;

}

void TPZSandlerExtended::ComputeJ2(TPZVec<STATE> stress, STATE &J2)const {
    STATE sig1, sig2, sig3;
    sig1 = stress[0];
    sig2 = stress[1];
    sig3 = stress[2];
    J2 = (2. * sig1 * sig2 + pow(sig1 + (-sig1 - sig2 - sig3) / 3., 2.) +
            pow(sig2 + (-sig1 - sig2 - sig3) / 3., 2.) + 2 * sig1 * sig3 + 2. * sig2 * sig3 +
            pow((-sig1 - sig2 - sig3) / 3. + sig3, 2.)) / 2.;
}

void TPZSandlerExtended::ApplyStrainComputeElasticStress(TPZVec<STATE> &strain, TPZVec<STATE> &stress)const {
    STATE sig1, sig2, sig3, s1, s2, s3;
    sig1 = strain[0];
    sig2 = strain[1];
    sig3 = strain[2];

    s1 = sig1 - (1. / 3.)*(sig1 + sig2 + sig3);
    s2 = sig2 - (1. / 3.)*(sig1 + sig2 + sig3);
    s3 = sig3 - (1. / 3.)*(sig1 + sig2 + sig3);

    stress[0] = s1 * (2 * fG) + fK * (sig1 + sig2 + sig3);
    stress[1] = s2 * (2 * fG) + fK * (sig1 + sig2 + sig3);
    stress[2] = s3 * (2 * fG) + fK * (sig1 + sig2 + sig3);
}

void TPZSandlerExtended::ApplyStressComputeElasticStrain(TPZVec<STATE> &stress, TPZVec<STATE> &strain)const {
    STATE sig1, sig2, sig3, s1, s2, s3;
    sig1 = stress[0];
    sig2 = stress[1];
    sig3 = stress[2];

    s1 = sig1 - (1. / 3.)*(sig1 + sig2 + sig3);
    s2 = sig2 - (1. / 3.)*(sig1 + sig2 + sig3);
    s3 = sig3 - (1. / 3.)*(sig1 + sig2 + sig3);

    strain[0] = s1 / (2. * fG)+(sig1 + sig2 + sig3) / (9. * fK);
    strain[1] = s2 / (2. * fG)+(sig1 + sig2 + sig3) / (9. * fK);
    strain[2] = s3 / (2. * fG)+(sig1 + sig2 + sig3) / (9. * fK);
}

/**
 * Imposes the specified strain tensor and returns the correspondent stress state.
 *
 * @param[in] epsTotal Imposed total strain tensor
 * @param[out] sigma Resultant stress
 */
void TPZSandlerExtended::ApplyStrainComputeSigma(TPZVec<STATE> &epst, TPZVec<STATE> &epsp, STATE & kprev, TPZVec<STATE> &epspnext, TPZVec<STATE> &stressnext, STATE & knext, TPZFMatrix<REAL> * tangent) const {

    bool require_tangent_Q = true;
    if (!tangent) {
        require_tangent_Q = false;
    }
    
#ifdef PZDEBUG
    // Check for required dimensions of tangent
    if (!(tangent->Rows() == 6 && tangent->Cols() == 6)) {
        std::cerr << "Unable to compute the tangent operator. Required tangent array dimensions are 6x6." << std::endl;
        DebugStop();
    }
#endif
    
    if (require_tangent_Q) {
        DebugStop(); // implemented this functionality.
    }
    
    STATE trial_I1, I1proj;

    TPZManVector<STATE, 3> trial_stress(3), yield(2), deltastress(3), delepsp(3), epsT(epst);
    epsT[0] -= epsp[0];
    epsT[1] -= epsp[1];
    epsT[2] -= epsp[2];
    ApplyStrainComputeElasticStress(trial_stress, epsT);
    YieldFunction(trial_stress, kprev, yield);
    ComputeI1(trial_stress, trial_I1);

    if ((yield[1] <= 0 && trial_I1 <= kprev) || (yield[0] <= 0 && trial_I1 > kprev)) {
        epspnext = epsp;
        knext = kprev;
        stressnext = trial_stress;
        //cout<<"\n elastic "<<endl;
    } else {
        //cout<<"\n plastic "<<endl;
        if (yield[1] > 0 && trial_I1 < kprev) {
            //cout<<"\n F2 "<<endl;
            ProjectF2(trial_stress, kprev, stressnext, knext);
        } else {
            //cout<<"\n F1 "<<endl;
            ProjectF1(trial_stress, kprev, stressnext, knext);
            ComputeI1(stressnext, I1proj);
            if (I1proj < knext) {
                //cout<<"\n Ring "<<endl;
                ProjectRing(trial_stress, kprev, stressnext, knext);
            }
        }

        for (int i = 0; i < 3; i++) {
            deltastress[i] = trial_stress[i] - stressnext[i];
        }

        ApplyStressComputeElasticStrain(deltastress, delepsp);

        for (int i = 0; i < 3; i++) {
            epspnext[i] = epsp[i] + delepsp[i];
        }
    }
}

void TPZSandlerExtended::ProjectSigma(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &sigproj, STATE &kproj, int &m_type, TPZFMatrix<REAL> * gradient) const {
    
    bool require_gradient_Q = (gradient != NULL);
    
#ifdef PZDEBUG
    if (require_gradient_Q) {
        // Check for required dimensions of tangent
        if (!(gradient->Rows() == 3 && gradient->Cols() == 3)) {
            std::cerr << "Unable to compute the gradient operator. Required gradient array dimensions are 3x3." << std::endl;
            DebugStop();
        }
    }
#endif
    
    //Firstk(epspv,k0);
    TPZManVector<STATE, 2> yield(2);
    STATE I1 = sigtrial[0] + sigtrial[1] + sigtrial[2];
    REAL J2 = (1.0/3.0) * (sigtrial[0]*sigtrial[0] + sigtrial[1]*sigtrial[1] + sigtrial[2]*sigtrial[2] - sigtrial[1]*sigtrial[2] - sigtrial[0]*sigtrial[2] - sigtrial[0]*sigtrial[1]);
    
    YieldFunction(sigtrial, kprev, yield);

    if (I1 < kprev) {
        if (yield[1] > 0.0 || IsZero(yield[1]) ) {
            
            m_type = 2; // cap behavior
//            std::cout << "Projecting on Cap " << std::endl;
            ProjectF2(sigtrial, kprev, sigproj, kproj);
            
            STATE proj_i1 = sigproj[0] + sigproj[1] + sigproj[2];
#ifdef PZDEBUG
            if (proj_i1 > kproj) {
                std::ostringstream stringbuilder("TPZSandlerExtended::ProjectSigma: proj_i1 > kproj with proj_i1 = ");
                stringbuilder << proj_i1 << " and kproj = " << kproj << ".";
                throw TPZInconsistentStateException(stringbuilder.str());
            }
#endif
            if (require_gradient_Q) {
                ComputeCapTangent(sigtrial, kprev, sigproj, kproj, gradient);
            }
            return;
    
        } else {
            m_type = 0; // Elastic behaviour
            sigproj = sigtrial;
            kproj = kprev;
            if (require_gradient_Q) {
                gradient->Identity();
            }
            return;
        }
    } else {

        if (yield[0] > 0.0 || IsZero(yield[0]) ) {

            bool failure_validity_Q = I1 > kprev || IsZero(I1 - kprev);
            if (failure_validity_Q) {
                STATE normal_to_f1_at_last_k = NormalToF1(I1, kprev);
                bool covertex_validity_Q = normal_to_f1_at_last_k < sqrt(J2) || IsZero(normal_to_f1_at_last_k - sqrt(J2));
                if (covertex_validity_Q) {
                    m_type = 2; // cap behavior
                    ProjectCapCoVertex(sigtrial, kprev, sigproj, kproj);
                    if (require_gradient_Q) {
                        ComputeCapCoVertexTangent(sigtrial, kprev, sigproj, kproj, gradient);
                    }
                    return;
                }else{
                    
                    REAL apex_i1 = Apex();
                    bool inside_apex_region_Q = I1 > apex_i1 || IsZero(I1 - apex_i1);
                    if (inside_apex_region_Q) {
                        STATE normal_to_f1_at_appex = NormalToF1(I1, apex_i1);
                        bool apex_validity_Q = normal_to_f1_at_appex > sqrt(J2) || IsZero(normal_to_f1_at_appex - sqrt(J2));
                        if (apex_validity_Q) { // Tensile behavior
                            m_type = -1;
//                            std::cout << "Projecting on Apex " << std::endl;
                            ProjectApex(sigtrial, kprev, sigproj, kproj);
                            if (require_gradient_Q) {
                                gradient->Identity(); //  Because tangent here is zero matrix, it is preferred use the elastic one.
                            }
                            return;
                        }
                        
                        m_type = 1; // failure behavior
//                        std::cout << "Projecting on Failure " << std::endl;
                        ProjectF1(sigtrial, kprev, sigproj, kproj);
                        if (require_gradient_Q) {
                            ComputeFailureTangent(sigtrial, kprev, sigproj, kproj, gradient);
                        }
                        return;
                        
                    }
                    else{
                        m_type = 1; // failure behavior
                        ProjectF1(sigtrial, kprev, sigproj, kproj);
                        if (require_gradient_Q) {
                            ComputeFailureTangent(sigtrial, kprev, sigproj, kproj, gradient);
                        }
                        return;
                    }
                }
            }
            
            DebugStop();

        } else {
            m_type = 0; // elastic behaviour
            sigproj = sigtrial;
            kproj = kprev;
            if (require_gradient_Q) {
                gradient->Identity();
            }
            return;
        }
    }
}


void TPZSandlerExtended::ComputeCapTangent(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj, TPZFMatrix<REAL> * gradient) const {
    
    TPZFMatrix<STATE> Rot(3,3,0.0);
    TPZHWTools::GetRotMatrix(Rot);
    
    
    STATE fv    = F(kproj);
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE CG    = fE/(2.0*(1.0 + fnu));
    
    TPZManVector<STATE, 3> hw_space_xi_rho_beta(3);
    TPZManVector<STATE, 3> rhw_space_s1_s2_s3(3);;
    TPZHWTools::FromPrincipalToHWCyl(projected_stress, hw_space_xi_rho_beta);
    TPZHWTools::FromPrincipalToHWCart(projected_stress, rhw_space_s1_s2_s3);

    STATE k = kproj;
    STATE i1v = sqrt(3)*hw_space_xi_rho_beta[0];
    STATE ratio = (k-i1v)/(fR*fv);
    
    bool cap_vertex_validity_Q = fabs(ratio - 1.0) <= ftol;
    if (cap_vertex_validity_Q) {
        ComputeCapVertexTangent(trial_stress, kprev, projected_stress, kproj, gradient);
        return;
    }
    
    STATE theta = acos(ratio);
    STATE beta = atan2(rhw_space_s1_s2_s3[2],rhw_space_s1_s2_s3[1]);
    
    // Computing the jacobian and the inverse
    TPZFMatrix<STATE> jac(3,3,0.0),jac_inv(3,3,0.0);
    Jacobianf2(trial_stress, theta, beta, k, jac); // Jacobian
    TPZHWTools::A3x3Inverse(jac, jac_inv); // Jacobian inverse
    
    TPZFMatrix<STATE> d_res_d_sig_trial(3,3,0.0);
    TPZFMatrix<STATE> d_sig_proj_d_theta_beta_kappa(3,3,0.0);
    TPZFMatrix<STATE> d_theta_beta_kappa_d_sig_trial(3,3,0.0);
    TPZFMatrix<STATE> d_sig_proj_d_sig_trial(3,3,0.0);

    // Derivative for the cap residue respect to trial stresses
    d_res_d_sig_trial(0,0) = (-2*(fA - exp(fB*k)*fC)*fR*sin(theta))/(3.*sqrt(3)*CK);
    d_res_d_sig_trial(0,1) = (-2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*cos(beta)*cos(theta))/(CG*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));
    d_res_d_sig_trial(0,2) = (-2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*cos(theta)*sin(beta))/(CG*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));
    d_res_d_sig_trial(1,0) = 0;
    d_res_d_sig_trial(1,1) = (2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*(2*(-1 + fPsi)*cos(2*beta) + (-1 + fPsi)*cos(4*beta) +
                                         (1 + fPsi)*sin(beta))*sin(theta))/(CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    d_res_d_sig_trial(1,2) = (-2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*cos(beta)*
     (1 + fPsi + 6*(-1 + fPsi)*sin(beta) - 2*(-1 + fPsi)*sin(3*beta))*sin(theta))/
    (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    d_res_d_sig_trial(2,0) = sqrt(3);
    d_res_d_sig_trial(2,1) = 0;
    d_res_d_sig_trial(2,2) = 0;
    
    
    
    // Derivative for the projected stresses respect to internal variables
    
    d_sig_proj_d_theta_beta_kappa(0,0) = ((fA - exp(fB*k)*fC)*(4*sqrt(3)*fPsi*cos(beta)*cos(theta) +
                                                               fR*(1 + fPsi + (-1 + fPsi)*sin(3*beta))*sin(theta)))/(3.*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));
    d_sig_proj_d_theta_beta_kappa(0,1) = (-4*(fA - exp(fB*k)*fC)*fPsi*(2*(-1 + fPsi)*cos(2*beta) + (-1 + fPsi)*cos(4*beta) +
                                                                       (1 + fPsi)*sin(beta))*sin(theta))/(sqrt(3)*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    d_sig_proj_d_theta_beta_kappa(0,2) = (1 + exp(fB*k)*fB*fC*fR*cos(theta) - (4*sqrt(3)*exp(fB*k)*fB*fC*fPsi*cos(beta)*sin(theta))/
                                          (1 + fPsi + (-1 + fPsi)*sin(3*beta)))/3.;
    
    
    d_sig_proj_d_theta_beta_kappa(1,0) = ((fA - exp(fB*k)*fC)*(-2*sqrt(3)*fPsi*cos(beta)*cos(theta) + 6*fPsi*cos(theta)*sin(beta) +
                                                               fR*(1 + fPsi + (-1 + fPsi)*sin(3*beta))*sin(theta)))/(3.*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));
    d_sig_proj_d_theta_beta_kappa(1,1) = (2*(fA - exp(fB*k)*fC)*fPsi*(3*(1 + fPsi)*cos(beta) + 2*sqrt(3)*(-1 + fPsi)*cos(2*beta) -
                                                                      sqrt(3)*cos(4*beta) + sqrt(3)*fPsi*cos(4*beta) + sqrt(3)*sin(beta) + sqrt(3)*fPsi*sin(beta) -
                                                                      6*sin(2*beta) + 6*fPsi*sin(2*beta) + 3*sin(4*beta) - 3*fPsi*sin(4*beta))*sin(theta))/
    (3.*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    d_sig_proj_d_theta_beta_kappa(1,2) = (1 + fPsi + (-1 + fPsi)*sin(3*beta) + exp(fB*k)*fB*fC*fR*cos(theta)*
                                          (1 + fPsi + (-1 + fPsi)*sin(3*beta)) + 2*sqrt(3)*exp(fB*k)*fB*fC*fPsi*cos(beta)*sin(theta) -
                                          6*exp(fB*k)*fB*fC*fPsi*sin(beta)*sin(theta))/(3.*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));

    
    d_sig_proj_d_theta_beta_kappa(2,0) = ((fA - exp(fB*k)*fC)*(-2*sqrt(3)*fPsi*cos(beta)*cos(theta) - 6*fPsi*cos(theta)*sin(beta) +
                                                               fR*(1 + fPsi + (-1 + fPsi)*sin(3*beta))*sin(theta)))/(3.*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));
    d_sig_proj_d_theta_beta_kappa(2,1) = (2*(fA - exp(fB*k)*fC)*fPsi*(-3*(1 + fPsi)*cos(beta) + 2*sqrt(3)*(-1 + fPsi)*cos(2*beta) -
                                 sqrt(3)*cos(4*beta) + sqrt(3)*fPsi*cos(4*beta) + sqrt(3)*sin(beta) + sqrt(3)*fPsi*sin(beta) +
                                 6*sin(2*beta) - 6*fPsi*sin(2*beta) - 3*sin(4*beta) + 3*fPsi*sin(4*beta))*sin(theta))/
    (3.*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    d_sig_proj_d_theta_beta_kappa(2,2) = (1 + fPsi + (-1 + fPsi)*sin(3*beta) + exp(fB*k)*fB*fC*fR*cos(theta)*
                                          (1 + fPsi + (-1 + fPsi)*sin(3*beta)) + 2*sqrt(3)*exp(fB*k)*fB*fC*fPsi*cos(beta)*sin(theta) +
                                          6*exp(fB*k)*fB*fC*fPsi*sin(beta)*sin(theta))/(3.*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));

    // Composing the gradient
    jac_inv *= -1.0;
    jac_inv.Multiply(d_res_d_sig_trial, d_theta_beta_kappa_d_sig_trial);//(ok)
    d_sig_proj_d_theta_beta_kappa.Multiply(d_theta_beta_kappa_d_sig_trial, d_sig_proj_d_sig_trial);
    d_sig_proj_d_sig_trial.Multiply(Rot, *gradient);

    
}

/// Compute the derivative of the projected stresses respect to trial stresses (tangent) over the cap
void TPZSandlerExtended::ComputeCapVertexTangent(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj, TPZFMatrix<REAL> * gradient) const{

    TPZFMatrix<STATE> Rot(3,3,0.0);
    TPZHWTools::GetRotMatrix(Rot);
    
    STATE fv    = F(kproj);
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE CG    = fE/(2.0*(1.0 + fnu));
    
    TPZManVector<STATE, 3> hw_space_xi_rho_beta(3);
    TPZManVector<STATE, 3> rhw_space_s1_s2_s3(3);;
    TPZHWTools::FromPrincipalToHWCyl(projected_stress, hw_space_xi_rho_beta);
    TPZHWTools::FromPrincipalToHWCart(projected_stress, rhw_space_s1_s2_s3);
    
    STATE k = kproj;
    STATE theta = 0.0;
    
    // Computing the jacobian and the inverse
    STATE jac,jac_inv;
    Jacobianf2Vertex(trial_stress, k, jac); // Jacobian
    jac_inv = 1.0 / jac; // Jacobian inverse
    
    TPZFMatrix<STATE> d_res_d_sig_trial(1,3,0.0);
    TPZFMatrix<STATE> d_sig_proj_d_kappa(3,1,0.0);
    TPZFMatrix<STATE> d_kappa_d_sig_trial(1,3,0.0);
    TPZFMatrix<STATE> d_sig_proj_d_sig_trial(3,3,0.0);
    
    // Derivative for the Covertex residue respect to trial stresses (1x3)
    d_res_d_sig_trial(0,0) = sqrt(3);
    d_res_d_sig_trial(0,1) =0;
    d_res_d_sig_trial(0,2) =0;
    
    // Derivative for the projected stresses respect to internal variables beta and k (3x1)
    d_sig_proj_d_kappa(0,0) = (1 + exp(fB*k)*fB*fC*fR*cos(theta))/3.;
    d_sig_proj_d_kappa(1,0) = (1 + exp(fB*k)*fB*fC*fR*cos(theta))/3.;
    d_sig_proj_d_kappa(2,0) = (1 + exp(fB*k)*fB*fC*fR*cos(theta))/3.;
    
    // Composing the gradient
    jac_inv *= -1.0;
    d_kappa_d_sig_trial = jac_inv * d_res_d_sig_trial;
    d_sig_proj_d_kappa.Multiply(d_kappa_d_sig_trial, d_sig_proj_d_sig_trial);
    d_sig_proj_d_sig_trial.Multiply(Rot, *gradient);
    
}

/// Compute the derivative of the projected stresses respect to trial stresses (tangent) over the cap
void TPZSandlerExtended::ComputeCapCoVertexTangent(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj, TPZFMatrix<REAL> * gradient) const{
    
    TPZFMatrix<STATE> Rot(3,3,0.0);
    TPZHWTools::GetRotMatrix(Rot);
    
    STATE fv    = F(kproj);
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE CG    = fE/(2.0*(1.0 + fnu));
    
    TPZManVector<STATE, 3> hw_space_xi_rho_beta(3);
    TPZManVector<STATE, 3> rhw_space_s1_s2_s3(3);;
    TPZHWTools::FromPrincipalToHWCyl(projected_stress, hw_space_xi_rho_beta);
    TPZHWTools::FromPrincipalToHWCart(projected_stress, rhw_space_s1_s2_s3);
    
    STATE k = kproj;
    STATE i1v = sqrt(3)*hw_space_xi_rho_beta[0];
    STATE theta = M_PI_2;
    STATE beta = atan2(rhw_space_s1_s2_s3[2],rhw_space_s1_s2_s3[1]);
    
    // Computing the jacobian and the inverse
    TPZFNMatrix<4, STATE> jac(2,2,0.0),jac_inv(2,2,0.0);
    Jacobianf2CoVertex(trial_stress, beta, k, jac); // Jacobian
    TPZHWTools::A2x2Inverse(jac, jac_inv); // Jacobian inverse
    
    TPZFMatrix<STATE> d_res_d_sig_trial(2,3,0.0);
    TPZFMatrix<STATE> d_sig_proj_d_beta_kappa(3,2,0.0);
    TPZFMatrix<STATE> d_beta_kappa_d_sig_trial(2,3,0.0);
    TPZFMatrix<STATE> d_sig_proj_d_sig_trial(3,3,0.0);
    
    // Derivative for the Covertex residue respect to trial stresses (2x3)
    d_res_d_sig_trial(0,0) = 0;
    d_res_d_sig_trial(0,1) = (2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*(2*(-1 + fPsi)*cos(2*beta) + (-1 + fPsi)*cos(4*beta) +
                                         (1 + fPsi)*sin(beta)))/(CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    d_res_d_sig_trial(0,2) = (-2*sqrt(2)*(fA - exp(fB*k)*fC)*fPsi*cos(beta)*
     (1 + fPsi + 6*(-1 + fPsi)*sin(beta) - 2*(-1 + fPsi)*sin(3*beta)))/
    (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    
    d_res_d_sig_trial(1,0) = sqrt(3);
    d_res_d_sig_trial(1,1) = 0;
    d_res_d_sig_trial(1,2) = 0;
    
    // Derivative for the projected stresses respect to internal variables beta and k (3x2)
    
    d_sig_proj_d_beta_kappa(0,0) = (-4*(fA - exp(fB*k)*fC)*fPsi*(2*(-1 + fPsi)*cos(2*beta) + (-1 + fPsi)*cos(4*beta) +
                                  (1 + fPsi)*sin(beta)))/(sqrt(3)*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    d_sig_proj_d_beta_kappa(0,1) = (1.0/3.0) - (4*exp(fB*k)*fB*fC*fPsi*cos(beta))/
    (sqrt(3)*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));
    
    d_sig_proj_d_beta_kappa(1,0) = (2*(fA - exp(fB*k)*fC)*fPsi*(3*(1 + fPsi)*cos(beta) + 2*sqrt(3)*(-1 + fPsi)*cos(2*beta) -
                                 sqrt(3)*cos(4*beta) + sqrt(3)*fPsi*cos(4*beta) + sqrt(3)*sin(beta) + sqrt(3)*fPsi*sin(beta) -
                                 6*sin(2*beta) + 6*fPsi*sin(2*beta) + 3*sin(4*beta) - 3*fPsi*sin(4*beta)))/
    (3.*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    d_sig_proj_d_beta_kappa(1,1) = (1 + fPsi + 2*sqrt(3)*exp(fB*k)*fB*fC*fPsi*cos(beta) - 6*exp(fB*k)*fB*fC*fPsi*sin(beta) -
     sin(3*beta) + fPsi*sin(3*beta))/(3.*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));
    
    
    d_sig_proj_d_beta_kappa(2,0) = (2*(fA - exp(fB*k)*fC)*fPsi*(-3*(1 + fPsi)*cos(beta) + 2*sqrt(3)*(-1 + fPsi)*cos(2*beta) -
                                 sqrt(3)*cos(4*beta) + sqrt(3)*fPsi*cos(4*beta) + sqrt(3)*sin(beta) + sqrt(3)*fPsi*sin(beta) +
                                 6*sin(2*beta) - 6*fPsi*sin(2*beta) - 3*sin(4*beta) + 3*fPsi*sin(4*beta)))/
    (3.*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    d_sig_proj_d_beta_kappa(2,1) = (1 + fPsi + 2*sqrt(3)*exp(fB*k)*fB*fC*fPsi*cos(beta) + 6*exp(fB*k)*fB*fC*fPsi*sin(beta) -
     sin(3*beta) + fPsi*sin(3*beta))/(3.*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));
    
    
    // Composing the gradient
    jac_inv *= -1.0;
    jac_inv.Multiply(d_res_d_sig_trial, d_beta_kappa_d_sig_trial);//(ok)
    d_sig_proj_d_beta_kappa.Multiply(d_beta_kappa_d_sig_trial, d_sig_proj_d_sig_trial);
    d_sig_proj_d_sig_trial.Multiply(Rot, *gradient);
    
}

/// Compute the derivative of the projected stresses respect to trial stresses (tangent) over the failure
void TPZSandlerExtended::ComputeFailureTangent(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj, TPZFMatrix<REAL> * gradient) const{
    
    TPZFMatrix<STATE> Rot(3,3,0.0);
    TPZHWTools::GetRotMatrix(Rot);
    
    
    STATE fv    = F(kproj);
    STATE CK    = fE/(3.0*(1.0 - 2.0 *fnu));
    STATE CG    = fE/(2.0*(1.0 + fnu));
    
    TPZManVector<STATE, 3> hw_space_xi_rho_beta(3);
    TPZManVector<STATE, 3> rhw_space_s1_s2_s3(3);
    TPZHWTools::FromPrincipalToHWCyl(projected_stress, hw_space_xi_rho_beta);
    TPZHWTools::FromPrincipalToHWCart(projected_stress, rhw_space_s1_s2_s3);
    
    STATE i1 = sqrt(3.0)*rhw_space_s1_s2_s3[0];
    STATE beta = atan2(rhw_space_s1_s2_s3[2],rhw_space_s1_s2_s3[1]);
    STATE k = kproj;
    
    
    // Computing the jacobian and the inverse
    TPZFMatrix<STATE> jac(3,3,0.0),jac_inv(3,3,0.0);
    Jacobianf1(trial_stress, i1, beta, k, jac); // Jacobian
    TPZHWTools::A3x3Inverse(jac, jac_inv); // Jacobian inverse
    
    TPZFMatrix<STATE> d_res_d_sig_trial(3,3,0.0);
    TPZFMatrix<STATE> d_sig_proj_d_i1_beta_kappa(3,3,0.0);
    TPZFMatrix<STATE> d_i1_beta_kappa_d_sig_trial(3,3,0.0);
    TPZFMatrix<STATE> d_sig_proj_d_sig_trial(3,3,0.0);
    
    // Derivative for the cap residue respect to trial stresses
    d_res_d_sig_trial(0,0) =     -2/(3.*sqrt(3)*CK);
    d_res_d_sig_trial(0,1) = (2*sqrt(2)*exp(fB*i1)*fB*fC*fPsi*cos(beta))/(CG*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));
    d_res_d_sig_trial(0,2) = (2*sqrt(2)*exp(fB*i1)*fB*fC*fPsi*sin(beta))/(CG*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));
    
    d_res_d_sig_trial(1,0) = 0;
    d_res_d_sig_trial(1,1) = (2*sqrt(2)*(fA - exp(fB*i1)*fC)*fPsi*(2*(-1 + fPsi)*cos(2*beta) + (-1 + fPsi)*cos(4*beta) +
                                          (1 + fPsi)*sin(beta)))/(CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    d_res_d_sig_trial(1,2) = (-2*sqrt(2)*(fA - exp(fB*i1)*fC)*fPsi*cos(beta)*
     (1 + fPsi + 6*(-1 + fPsi)*sin(beta) - 2*(-1 + fPsi)*sin(3*beta)))/
    (CG*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    
    d_res_d_sig_trial(2,0) = sqrt(3);
    d_res_d_sig_trial(2,1) = 0;
    d_res_d_sig_trial(2,2) = 0;
    
    // Derivative for the projected stresses respect to internal variables
    d_sig_proj_d_i1_beta_kappa(0,0) = 1.0/3.0 - (4*exp(fB*i1)*fB*fC*fPsi*cos(beta))/(sqrt(3)*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));
    d_sig_proj_d_i1_beta_kappa(0,1) = (-4*(fA - exp(fB*i1)*fC)*fPsi*(2*(-1 + fPsi)*cos(2*beta) + (-1 + fPsi)*cos(4*beta) +
                                                                      (1 + fPsi)*sin(beta)))/(sqrt(3)*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    d_sig_proj_d_i1_beta_kappa(0,2) = 0;
    
    d_sig_proj_d_i1_beta_kappa(1,0) = (1 + fPsi + 2*sqrt(3)*exp(fB*i1)*fB*fC*fPsi*cos(beta) - 6*exp(fB*i1)*fB*fC*fPsi*sin(beta) -
                                        sin(3*beta) + fPsi*sin(3*beta))/(3.*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));
    d_sig_proj_d_i1_beta_kappa(1,1) = (2*(fA - exp(fB*i1)*fC)*fPsi*(3*(1 + fPsi)*cos(beta) + 2*sqrt(3)*(-1 + fPsi)*cos(2*beta) -
                                                                     sqrt(3)*cos(4*beta) + sqrt(3)*fPsi*cos(4*beta) + sqrt(3)*sin(beta) + sqrt(3)*fPsi*sin(beta) -
                                                                     6*sin(2*beta) + 6*fPsi*sin(2*beta) + 3*sin(4*beta) - 3*fPsi*sin(4*beta)))/
    (3.*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    d_sig_proj_d_i1_beta_kappa(1,2) = 0;
    
    d_sig_proj_d_i1_beta_kappa(2,0) = (1 + fPsi + 2*sqrt(3)*exp(fB*i1)*fB*fC*fPsi*cos(beta) + 6*exp(fB*i1)*fB*fC*fPsi*sin(beta) -
                                        sin(3*beta) + fPsi*sin(3*beta))/(3.*(1 + fPsi + (-1 + fPsi)*sin(3*beta)));
    d_sig_proj_d_i1_beta_kappa(2,1) = (2*(fA - exp(fB*i1)*fC)*fPsi*(-3*(1 + fPsi)*cos(beta) + 2*sqrt(3)*(-1 + fPsi)*cos(2*beta) -
                                                                     sqrt(3)*cos(4*beta) + sqrt(3)*fPsi*cos(4*beta) + sqrt(3)*sin(beta) + sqrt(3)*fPsi*sin(beta) +
                                                                     6*sin(2*beta) - 6*fPsi*sin(2*beta) - 3*sin(4*beta) + 3*fPsi*sin(4*beta)))/
    (3.*pow(1 + fPsi + (-1 + fPsi)*sin(3*beta),2));
    d_sig_proj_d_i1_beta_kappa(2,2) = 0;
    
    // Composing the gradient
    jac_inv *= -1.0;
    jac_inv.Multiply(d_res_d_sig_trial, d_i1_beta_kappa_d_sig_trial);//(ok)
    d_sig_proj_d_i1_beta_kappa.Multiply(d_i1_beta_kappa_d_sig_trial, d_sig_proj_d_sig_trial);
    d_sig_proj_d_sig_trial.Multiply(Rot, *gradient);

//    jac.Print(std::cout);
//    jac_inv.Print(std::cout);
//    d_res_d_sig_trial.Print(std::cout);
//    d_i1_beta_kappa_d_sig_trial.Print(std::cout);
//    d_sig_proj_d_i1_beta_kappa.Print(std::cout);
//    d_sig_proj_d_sig_trial.Print(std::cout);
//    gradient->Print(std::cout);
    
}

void TPZSandlerExtended::ProjectSigmaDep(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &sigproj, STATE &kproj, TPZFMatrix<STATE> &GradSigma) const {
    
    DebugStop(); //  Deprecated method.
    
    STATE I1;
    //Firstk(epspv,k0);
    TPZManVector<STATE, 2> yield(2);
    I1 = sigtrial[0] + sigtrial[1] + sigtrial[2];

    YieldFunction(sigtrial, kprev, yield);
    bool threeEigEqual = (fabs(sigtrial[0] - sigtrial[1]) < ftol && fabs(sigtrial[1] - sigtrial[2]) < ftol);
    
    // kprev corresponde a L do artigo
    if (I1 < kprev) {
        if (yield[1] > 0. && !threeEigEqual) {
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "Projecting on F2, distinct eigenvalues " << sigtrial;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            ProjectF2(sigtrial, kprev, sigproj, kproj);
            // we can compute the tangent matrix
            TPZFNMatrix<9, STATE> dbetadsigtrial(3, 3), jacF2(3, 3), DF2cart(3, 3);
            STATE theta, beta;
            SurfaceParamF2(sigproj, kproj, theta, beta);
            GradF2SigmaTrial(sigtrial, theta, beta, kproj, kprev, dbetadsigtrial);
            Jacobianf2(sigtrial, theta, beta, kproj, jacF2);
            jacF2.Solve_LU(&dbetadsigtrial);
            DF2Cart(theta, beta, kproj, DF2cart);
            DF2cart.Multiply(dbetadsigtrial, GradSigma);
            GradSigma *= -1.;
        } else if (yield[1] > 0. && threeEigEqual) {
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "Projecting on F2, equal eigenvalues " << sigtrial;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            ProjectBetaConstF2(sigtrial, kprev, sigproj, kproj);
            // we can compute the tangent matrix
            STATE theta, beta;
            // compute theta beta as a function of sigproj, kproj
            SurfaceParamF2(sigproj, kproj, theta, beta);
            // theta should be Pi
            // for hydrostatic stress beta doesn't mean anything
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "Surface parameters for sigproj = " << sigproj << " kproj " << kproj << " theta " << theta;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            TPZManVector<STATE, 2> sigtrialIJ(2, 0.), sigprojIJ(2), ddistf2(2);
            sigtrialIJ[0] = sigtrial[0] + sigtrial[1] + sigtrial[2];
            DDistF2IJ(sigtrialIJ, theta, kproj, kprev, ddistf2);
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "Derivative of the distance function (should be zero) = " << ddistf2;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            TPZFNMatrix<4, STATE> JacThetaK(2, 2), JacSigtrIJ(2, 2), JacSigprojThetaK(2, 2);
            {
                TPZManVector<TFad<2, STATE>, 2> sigtrialIJFAD(2), ddistf2fad(2);
                TFad<2, STATE> thetafad(theta, 0), kprojfad(kproj, 1);
                sigtrialIJFAD[0].val() = sigtrialIJ[0];
                sigtrialIJFAD[1].val() = sigtrialIJ[1];
                DDistF2IJ(sigtrialIJFAD, thetafad, kprojfad, kprev, ddistf2fad);
                JacThetaK(0, 0) = ddistf2fad[0].d(0);
                JacThetaK(0, 1) = ddistf2fad[0].d(1);
                JacThetaK(1, 0) = ddistf2fad[1].d(0);
                JacThetaK(1, 1) = ddistf2fad[1].d(1);
            }
            {
                TPZManVector<TFad<2, STATE>, 2> sigtrialIJFAD(2), ddistf2fad(2);
                TFad<2, STATE> thetafad(theta), kprojfad(kproj);
                sigtrialIJFAD[0].val() = sigtrialIJ[0];
                sigtrialIJFAD[0].fastAccessDx(0) = 1.;
                sigtrialIJFAD[1].val() = sigtrialIJ[1];
                sigtrialIJFAD[1].fastAccessDx(1) = 1.;
                DDistF2IJ(sigtrialIJFAD, thetafad, kprojfad, kprev, ddistf2fad);
                JacSigtrIJ(0, 0) = ddistf2fad[0].d(0);
                JacSigtrIJ(0, 1) = ddistf2fad[0].d(1);
                JacSigtrIJ(1, 0) = ddistf2fad[1].d(0);
                JacSigtrIJ(1, 1) = ddistf2fad[1].d(1);
            }
            FromThetaKToSigIJ(theta, kproj, sigprojIJ);
            {
                TPZManVector<TFad<2, STATE>, 2> sigprojIJFAD(2);
                TFad<2, STATE> thetafad(theta, 0), kprojfad(kproj, 1);
                FromThetaKToSigIJ(thetafad, kprojfad, sigprojIJFAD);
                JacSigprojThetaK(0, 0) = sigprojIJFAD[0].d(0);
                JacSigprojThetaK(0, 1) = sigprojIJFAD[0].d(1);
                JacSigprojThetaK(1, 0) = sigprojIJFAD[1].d(0);
                JacSigprojThetaK(1, 1) = sigprojIJFAD[1].d(1);
            }
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "Derivative of the distanceIJ " << ddistf2 << std::endl;
                JacThetaK.Print("Derivative of distanceIJ with respect to theta, K", sout);
                JacSigtrIJ.Print("Derivative of distanceIJ with respect to sigtialIJ", sout);
                sout << "SigmaProjected IJ " << sigprojIJ << std::endl;
                JacSigprojThetaK.Print("Derivative of sigproj with respect to theta K", sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            std::list<int64_t> singular;
            JacThetaK.Solve_LU(&JacSigtrIJ, singular);
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "Negative of derivative of Theta,K with respect to sigtrIJ" << std::endl;
                JacSigtrIJ.Print("Derivative = ", sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            TPZFNMatrix<4, STATE> dsigprojdsigtr(2, 2);
            JacSigprojThetaK.Multiply(JacSigtrIJ, dsigprojdsigtr);
            dsigprojdsigtr *= -1.;
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                dsigprojdsigtr.Print("Derivative of SigprojIJ with respect to SigtrialIJ", sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    GradSigma(i, j) = dsigprojdsigtr(0, 0) / 3.;
                    if (i == j) {
                        GradSigma(i, j) += dsigprojdsigtr(1, 1)*2. / 3.;
                    } else {
                        GradSigma(i, j) -= dsigprojdsigtr(1, 1) / 3.;
                    }
                }
            }
        } else {
            sigproj = sigtrial;
            kproj = kprev;
            GradSigma.Identity();

        }
    } else {
        if (yield[0] > 0.) {
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "Projecting on F1";
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            ProjectF1(sigtrial, kprev, sigproj, kproj);
#ifdef PZDEBUG
            {
                TPZManVector<STATE, 2> yieldcopy(2);
                YieldFunction(sigproj, kproj, yieldcopy);
                if (fabs(yieldcopy[0]) > 1.e-3) {
                    DebugStop();
                }
            }
#endif

            I1 = 0.;
            for (int i = 0; i < 3; i++) {
                I1 += sigproj[i];
            }
            if (I1 < kproj) {
                ProjectRing(sigtrial, kprev, sigproj, kproj);


                // we can compute the tangent matrix
                TPZFNMatrix<9, STATE> dbetadsigtrial(3, 3), jacF2(3, 3), DF2cart(3, 3);
                STATE theta, beta;
                SurfaceParamF2(sigproj, kproj, theta, beta);
#ifdef PZDEBUG
                if (fabs(theta) - M_PI_2 > 1.e-8) {
                    DebugStop();
                }
#endif
                GradF2SigmaTrial(sigtrial, theta, beta, kproj, kprev, dbetadsigtrial);
                for (int i = 0; i < 3; i++) dbetadsigtrial(0, i) = 0.;
                Jacobianf2(sigtrial, theta, beta, kproj, jacF2);
                for (int i = 0; i < 3; i++) {
                    jacF2(i, 0) = 0.;
                    jacF2(0, 1) = 0.;
                }
                jacF2(0, 0) = 1.;
                jacF2.Solve_LU(&dbetadsigtrial);
                DF2Cart(theta, beta, kproj, DF2cart);
                for (int i = 0; i < 3; i++) DF2cart(i, 0) = 0.;
                DF2cart.Multiply(dbetadsigtrial, GradSigma);
                GradSigma *= -1.;

            } else {
                // we can compute the tangent matrix
                TPZFNMatrix<9, STATE> dbetadsigtrial(2, 3), jacF1(2, 2), DF1cart(3, 2);
                STATE xi, beta;
                SurfaceParamF1(sigproj, xi, beta);
                GradF1SigmaTrial(sigtrial, xi, beta, dbetadsigtrial);
                D2DistFunc1(sigtrial, xi, beta, jacF1);
                jacF1.Solve_LU(&dbetadsigtrial);
                DF1Cart(xi, beta, DF1cart);
                DF1cart.Multiply(dbetadsigtrial, GradSigma);
                GradSigma *= -1.;

            }
        } else {
#ifdef LOG4CXX
            {
                if (logger->isDebugEnabled()) {
                    std::stringstream sout;
                    sout << "Elastic Behaviour";
                    LOGPZ_DEBUG(logger, sout.str())
                }
            }
#endif
            // elastic behaviour
            sigproj = sigtrial;
            kproj = kprev;
            GradSigma.Identity();
        }
    }
}


#define Cos cos
#define Sin sin
#define Sqrt sqrt


/// Compute the derivative of the residual with respect to sigtrial

void TPZSandlerExtended::GradF1SigmaTrial(const TPZVec<STATE> &sigtrial, STATE xi, STATE beta, TPZFMatrix<STATE> &deriv) const {
    STATE s3beta = sin(3. * beta);
    STATE c3beta = cos(3. * beta);
    STATE Gamma = (1 + s3beta)+(1. - s3beta) / fPsi;
    STATE DGamma = 3. * c3beta * (1. - 1. / fPsi);
    STATE SQR3 = sqrt(3.);
    STATE FFI = fA - fPhi * SQR3 * xi - fC * exp(fB * SQR3 * xi);
    STATE NN = fN;
    deriv.Redim(2, 3);
    deriv(0, 0) = (-2 * (fG * Gamma - 6 * fK * (SQR3 * fA * fB - SQR3 * fB * FFI + SQR3 * fPhi - 3 * fB * fPhi * xi) * Cos(beta))) / (3. * SQR3 * fG * fK * Gamma);
    deriv(1, 0) = (4 * (FFI - NN)*(DGamma * Cos(beta) + Gamma * Sin(beta))) / (SQR3 * fG * Gamma * Gamma);

    deriv(0, 1) = (-2 * (SQR3 * fG * Gamma + 9 * fK * (fA * fB + fPhi - fB * (FFI + SQR3 * fPhi * xi)) * Cos(beta) -
            9 * fK * (SQR3 * fA * fB + SQR3 * fPhi - fB * (SQR3 * FFI + 3 * fPhi * xi)) * Sin(beta))) / (9. * fG * fK * Gamma);
    deriv(1, 1) = (-2 * (FFI - NN)*((SQR3 * DGamma + 3 * Gamma) * Cos(beta) + (-3 * DGamma + SQR3 * Gamma) * Sin(beta))) / (3. * fG * Gamma * Gamma);

    deriv(0, 2) = (-2 * (SQR3 * fG * Gamma + 9 * fK * (fA * fB + fPhi - fB * (FFI + SQR3 * fPhi * xi)) * Cos(beta) +
            9 * fK * (SQR3 * fA * fB + SQR3 * fPhi - fB * (SQR3 * FFI + 3 * fPhi * xi)) * Sin(beta))) / (9. * fG * fK * Gamma);
    deriv(1, 2) = (-2 * (FFI - NN)*((SQR3 * DGamma - 3 * Gamma) * Cos(beta) + (3 * DGamma + SQR3 * Gamma) * Sin(beta))) / (3. * fG * Gamma * Gamma);
}

/// Compute the derivative of the F2 residual with respect to sigtrial

void TPZSandlerExtended::GradF2SigmaTrial(const TPZVec<STATE> &sigtrial, STATE theta, STATE beta, STATE k, STATE kprev, TPZFMatrix<STATE> &deriv) const {
    STATE s3beta = sin(3. * beta);
    STATE c3beta = cos(3. * beta);
    STATE Gamma = (1 + s3beta)+(1. - s3beta) / fPsi;
    STATE DGamma = 3. * c3beta * (1. - 1. / fPsi);
    STATE SQR3 = sqrt(3.);
    STATE FFK = fA - fPhi * k - fC * exp(fB * k);
    deriv.Redim(3, 3);
    deriv(0, 0) = (-4 * (FFK - fN) * Cos(beta) * Cos(theta)) / (SQR3 * fG * Gamma) + (2 * FFK * fR * Sin(theta)) / (9. * fK);
    deriv(1, 0) = (4 * (FFK - fN)*(DGamma * Cos(beta) + Gamma * Sin(beta)) * Sin(theta)) / (SQR3 * fG * Gamma * Gamma);
    deriv(2, 0) = -1.;

    deriv(0, 1) = (2 * (3 * SQR3 * fK * (FFK - fN) * Cos(beta) * Cos(theta) - 9 * fK * (FFK - fN) * Cos(theta) * Sin(beta) + FFK * fG * fR * Gamma * Sin(theta))) / (9. * fG * fK * Gamma);
    deriv(1, 1) = (-2 * (FFK - fN)*((SQR3 * DGamma + 3 * Gamma) * Cos(beta) + (-3 * DGamma + SQR3 * Gamma) * Sin(beta)) * Sin(theta)) / (3. * fG * Gamma * Gamma);
    deriv(2, 1) = -1.;

    deriv(0, 2) = (2 * (3 * SQR3 * fK * (FFK - fN) * Cos(beta) * Cos(theta) + 9 * fK * (FFK - fN) * Cos(theta) * Sin(beta) + FFK * fG * fR * Gamma * Sin(theta))) / (9. * fG * fK * Gamma);
    deriv(1, 2) = (-2 * (FFK - fN)*((SQR3 * DGamma - 3 * Gamma) * Cos(beta) + (3 * DGamma + SQR3 * Gamma) * Sin(beta)) * Sin(theta)) / (3. * fG * Gamma * Gamma);
    deriv(2, 2) = -1.;

}

void TPZSandlerExtended::TaylorCheckDistF1(const TPZVec<STATE> &sigmatrial, STATE xi, STATE beta, TPZVec<STATE> &xnorm,
        TPZVec<STATE> &errnorm) const {
    STATE deltaxi = 0.4;
    STATE deltabeta = 0.05;
    STATE dist0 = DistF1(sigmatrial, xi, beta);
    TPZFNMatrix<4, STATE> jac(2, 1);
    DDistFunc1(sigmatrial, xi, beta, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i = 1; i <= 10; i++) {
        STATE diffxi = deltaxi * i / 10.;
        STATE xinext = xi + diffxi;
        STATE diffbeta = deltabeta * i / 10.;
        STATE betanext = beta + diffbeta;
        STATE distnext = DistF1(sigmatrial, xinext, betanext);
        STATE distguess = dist0 + jac(0, 0) * diffxi + jac(1, 0) * diffbeta;
        xnorm[i - 1] = sqrt(diffxi * diffxi + diffbeta * diffbeta);
        errnorm[i - 1] = fabs(distnext - distguess);
    }

}

void TPZSandlerExtended::TaylorCheckDDistF1(const TPZVec<STATE> &sigmatrial, STATE xi, STATE beta, TPZVec<STATE> &xnorm,
        TPZVec<STATE> &errnorm) const {
    STATE deltaxi = 0.4;
    STATE deltabeta = 0.05;
    TPZFNMatrix<2, STATE> res0(2, 1), resid(2, 1), residguess(2, 1), diff(2, 1);
    TPZFNMatrix<4, STATE> jac(2, 2);
    DDistFunc1(sigmatrial, xi, beta, res0);
    D2DistFunc1(sigmatrial, xi, beta, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i = 1; i <= 10; i++) {
        STATE diffxi = deltaxi * i / 10.;
        STATE xinext = xi + diffxi;
        STATE diffbeta = deltabeta * i / 10.;
        diff(0) = diffxi;
        diff(1) = diffbeta;
        STATE betanext = beta + diffbeta;
        DDistFunc1(sigmatrial, xinext, betanext, resid);
        jac.Multiply(diff, residguess);
        residguess += res0;
        xnorm[i - 1] = Norm(diff);
        errnorm[i - 1] = Norm(resid - residguess);
    }

}

void TPZSandlerExtended::TaylorCheckDDistF1DSigtrial(const TPZVec<STATE> &sigmatrial, STATE xi, STATE beta, TPZVec<STATE> &xnorm,
        TPZVec<STATE> &errnorm) const {
    TPZManVector<STATE, 3> deltasigma(3, 0.3);
    deltasigma[1] *= -1.;
    deltasigma[2] *= 0.7;

    TPZFNMatrix<2, STATE> res0(2, 1), resid(2, 1), residguess(2, 1), diff(3, 1);
    TPZFNMatrix<6, STATE> jac(2, 3);
    DDistFunc1(sigmatrial, xi, beta, res0);
    GradF1SigmaTrial(sigmatrial, xi, beta, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i = 1; i <= 10; i++) {
        for (int j = 0; j < 3; j++) diff(j) = deltasigma[j] * i / 10.;
        TPZManVector<STATE, 3> sigmanext(3);
        for (int j = 0; j < 3; j++) sigmanext[j] = sigmatrial[j] + diff(j);
        DDistFunc1(sigmanext, xi, beta, resid);
        jac.Multiply(diff, residguess);
        residguess += res0;
        xnorm[i - 1] = Norm(diff);
        errnorm[i - 1] = Norm(resid - residguess);
    }

}

void TPZSandlerExtended::ConvergenceRate(TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm, TPZVec<STATE> &convergence) {
    convergence.resize(xnorm.size() - 1);
    for (int i = 1; i < xnorm.size(); i++) {
        convergence[i - 1] = log(errnorm[i] / errnorm[i - 1]) / log(xnorm[i] / xnorm[i - 1]);
    }
}

void TPZSandlerExtended::TaylorCheckDistF2(const TPZVec<STATE> &sigmatrial, STATE theta, STATE beta, STATE k, STATE kprev, TPZVec<STATE> &xnorm,
        TPZVec<STATE> &errnorm) const {
    STATE deltatheta = 0.4;
    STATE deltabeta = 0.05;
    STATE deltak = 0.;
    STATE dist0 = DistF2(sigmatrial, theta, beta, k);
    TPZFNMatrix<4, STATE> jac(3, 1);
    TPZManVector<STATE> fxnvec(3);
    Res2(sigmatrial, theta, beta, k, kprev, fxnvec);
    for (int kk = 0; kk < 3; kk++) jac(kk, 0) = fxnvec[kk];
    //    DDistFunc2(sigmatrial, theta, beta, k, kprev, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i = 1; i <= 10; i++) {
        STATE difftheta = deltatheta * i / 10.;
        STATE thetanext = theta + difftheta;
        STATE diffbeta = deltabeta * i / 10.;
        STATE betanext = beta + diffbeta;
        STATE diffk = deltak * i / 10.;
        STATE knext = k + deltak;
        STATE distnext = DistF2(sigmatrial, thetanext, betanext, knext);
        STATE distguess = dist0 + jac(0, 0) * difftheta + jac(1, 0) * diffbeta + jac(2, 0) * diffk;
        xnorm[i - 1] = sqrt(difftheta * difftheta + diffbeta * diffbeta + diffk * diffk);
        errnorm[i - 1] = fabs(distnext - distguess);
    }

}

void TPZSandlerExtended::TaylorCheckDDistF2(const TPZVec<STATE> &sigmatrial, STATE theta, STATE beta, STATE k, STATE kprev, TPZVec<STATE> &xnorm,
        TPZVec<STATE> &errnorm) const {
    STATE deltatheta = 0.4;
    STATE deltabeta = 0.05;
    STATE deltak = 0.5;
    TPZFNMatrix<3, STATE> res0(3, 1), resid(3, 1), residguess(3, 1), diff(3, 1);
    TPZFNMatrix<9, STATE> jac(3, 3);
    TPZManVector<STATE> fxnvec(3);
    Res2(sigmatrial, theta, beta, k, kprev, fxnvec);
    for (int kk = 0; kk < 3; kk++) res0(kk, 0) = fxnvec[kk];
    //    DDistFunc2(sigmatrial, theta, beta, k, kprev, res0);
    Jacobianf2(sigmatrial, theta, beta, k, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i = 1; i <= 10; i++) {
        STATE difftheta = deltatheta * i / 10.;
        STATE thetanext = theta + difftheta;
        STATE diffbeta = deltabeta * i / 10.;
        STATE betanext = beta + diffbeta;
        STATE diffk = deltak * i / 10.;
        STATE knext = k + diffk;
        diff(0) = difftheta;
        diff(1) = diffbeta;
        diff(2) = diffk;
        TPZManVector<STATE> fxnvec(3);
        Res2(sigmatrial, thetanext, betanext, knext, kprev, fxnvec);
        for (int k = 0; k < 3; k++) resid(k, 0) = fxnvec[k];
        //        DDistFunc2(sigmatrial, thetanext, betanext, knext, kprev, resid);
        jac.Multiply(diff, residguess);
        residguess += res0;
        xnorm[i - 1] = Norm(diff);
        errnorm[i - 1] = Norm(resid - residguess);
    }

}

/// teste da derivada D(ResF2)/D(sigtrial)

void TPZSandlerExtended::TaylorCheckDDistF2DSigtrial(const TPZVec<STATE> &sigmatrial, STATE theta, STATE beta, STATE k, STATE kprev, TPZVec<STATE> &xnorm,
        TPZVec<STATE> &errnorm) const {
    TPZManVector<STATE, 3> deltasigma(3, 0.3);
    deltasigma[1] *= -1.;
    deltasigma[2] *= 0.7;

    TPZFNMatrix<3, STATE> res0(3, 1), resid(3, 1), residguess(3, 1), diff(3, 1);
    TPZFNMatrix<9, STATE> jac(3, 3);
    TPZManVector<STATE> fxnvec(3);
    Res2(sigmatrial, theta, beta, k, kprev, fxnvec);
    for (int kk = 0; kk < 3; kk++) res0(kk, 0) = fxnvec[kk];
    //    DDistFunc2(sigmatrial, theta, beta, k, kprev, res0);
    GradF2SigmaTrial(sigmatrial, theta, beta, k, kprev, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i = 1; i <= 10; i++) {
        for (int j = 0; j < 3; j++) diff(j) = deltasigma[j] * i / 10.;
        TPZManVector<STATE, 3> sigmanext(3);
        for (int j = 0; j < 3; j++) sigmanext[j] = sigmatrial[j] + diff(j);
        TPZManVector<STATE> fxnvec(3);
        Res2(sigmanext, theta, beta, k, kprev, fxnvec);
        for (int k = 0; k < 3; k++) resid(k, 0) = fxnvec[k];
        //        DDistFunc2(sigmanext, theta, beta,k,kprev,resid);
        jac.Multiply(diff, residguess);
        residguess += res0;
        xnorm[i - 1] = Norm(diff);
        errnorm[i - 1] = Norm(resid - residguess);
    }

}

#include "pzvec_extras.h"

void TPZSandlerExtended::CheckCoordinateTransformation(TPZVec<STATE> &cart) {
    TPZManVector<STATE, 3> HWCart(3), HWCyl(3), HWCart2(3), Cart2(3);
    TPZHWTools::FromPrincipalToHWCart(cart, HWCart);
    TPZHWTools::FromPrincipalToHWCyl(cart, HWCyl);
    TPZHWTools::FromHWCylToHWCart(HWCyl, HWCart2);
    TPZHWTools::FromHWCylToPrincipal(HWCyl, Cart2);
    REAL dist1 = dist(cart, Cart2);
    REAL dist2 = dist(HWCart, HWCart2);
    cout << __FUNCTION__ << " dist1 = " << dist1 << " dist2 = " << dist2 << endl;

}

/// Compute the derivative of the stress (principal s;tresses) as a function of xi and beta

void TPZSandlerExtended::DF1Cart(STATE xi, STATE beta, TPZFMatrix<STATE> &DF1) const {
    STATE s3beta = sin(3. * beta);
    STATE c3beta = cos(3. * beta);
    STATE Gamma = (1 + s3beta)+(1. - s3beta) / fPsi;
    STATE DGamma = 3. * c3beta * (1. - 1. / fPsi);
    STATE SQR3 = sqrt(3.);
    STATE FFI = fA - fPhi * SQR3 * xi - fC * exp(fB * SQR3 * xi);
    DF1.Resize(3, 2);
    DF1(0, 0) = 1 / SQR3 + (4 * (-(SQR3 * fPhi) - SQR3 * fB * (fA - FFI - SQR3 * fPhi * xi)) * Cos(beta)) /
            (SQR3 * Gamma);
    DF1(1, 0) = 1 / SQR3 + (4 * (-(SQR3 * fPhi) -
            SQR3 * fB * (fA - FFI - SQR3 * fPhi * xi)) * Sin(beta - M_PI / 6.)) / (SQR3 * Gamma);
    DF1(2, 0) =
            1 / SQR3 - (4 * (-(SQR3 * fPhi) - SQR3 * fB * (fA - FFI - SQR3 * fPhi * xi)) *
            Sin(beta + M_PI / 6.)) / (SQR3 * Gamma);

    DF1(0, 1) = (-4 * DGamma * (FFI - fN) * Cos(beta)) / (SQR3 * Gamma * Gamma) -
            (4 * (FFI - fN) * Sin(beta)) / (SQR3 * Gamma);
    DF1(1, 1) = (4 * (FFI - fN) * Cos(beta - M_PI / 6.)) / (SQR3 * Gamma) -
            (4 * DGamma * (FFI - fN) * Sin(beta - M_PI / 6.)) / (SQR3 * Gamma * Gamma);

    DF1(2, 1) = (-4 * (FFI - fN) * Cos(beta + M_PI / 6.)) / (SQR3 * Gamma) +
            (4 * DGamma * (FFI - fN) * Sin(beta + M_PI / 6.)) / (SQR3 * Gamma * Gamma);
}

/// Compute the derivative of the stress (principal stresses) as a function of xi and beta

void TPZSandlerExtended::DF2Cart(STATE theta, STATE beta, STATE k, TPZFMatrix<STATE> &DF2) const {
    STATE s3beta = sin(3. * beta);
    STATE c3beta = cos(3. * beta);
    STATE Gamma = (1 + s3beta)+(1. - s3beta) / fPsi;
    STATE DGamma = 3. * c3beta * (1. - 1. / fPsi);
    STATE SQR3 = sqrt(3.);
    STATE FFK = fA - fPhi * k - fC * exp(fB * k);
    DF2.Resize(3, 3);
    DF2(0, 0) = (4 * (FFK - fN) * Cos(beta) * Cos(theta)) / (SQR3 * Gamma) - (FFK * fR * Sin(theta)) / 3.;
    DF2(1, 0) = (4 * (FFK - fN) * Cos(theta) * Sin(beta - M_PI / 6.)) / (SQR3 * Gamma) - (FFK * fR * Sin(theta)) / 3.;
    DF2(2, 0) = (4 * (-FFK + fN) * Cos(theta) * Sin(beta + M_PI / 6.)) / (SQR3 * Gamma) - (FFK * fR * Sin(theta)) / 3.;

    DF2(0, 1) = (-4 * (FFK - fN)*(DGamma * Cos(beta) + Gamma * Sin(beta)) * Sin(theta)) / (SQR3 * Gamma * Gamma);
    DF2(1, 1) = (4 * (FFK - fN)*(Gamma * Cos(beta - M_PI / 6.) - DGamma * Sin(beta - M_PI / 6.)) * Sin(theta)) / (SQR3 * Gamma * Gamma);
    DF2(2, 1) = (-4 * (FFK - fN)*(Gamma * Cos(beta + M_PI / 6.) - DGamma * Sin(beta + M_PI / 6.)) * Sin(theta)) / (SQR3 * Gamma * Gamma);

    DF2(0, 2) = (Gamma - (fA * fB + fPhi - fB * (FFK + fPhi * k))*
            (fR * Gamma * Cos(theta) + 4 * SQR3 * Cos(beta) * Sin(theta))) / (3. * Gamma);
    DF2(1, 2) = (Gamma - (fA * fB + fPhi - fB * (FFK + fPhi * k))*
            (fR * Gamma * Cos(theta) + 4 * SQR3 * Sin(beta - M_PI / 6.) * Sin(theta))) / (3. * Gamma);
    DF2(2, 2) = (Gamma - (fA * fB + fPhi - fB * (FFK + fPhi * k))*
            (fR * Gamma * Cos(theta) - 4 * SQR3 * Sin(beta + M_PI / 6.) * Sin(theta))) / (3. * Gamma);
}

void TPZSandlerExtended::TaylorCheckDF1Cart(STATE xi, STATE beta, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const {
    STATE deltaxi = 0.4;
    STATE deltabeta = 0.05;
    TPZFNMatrix<3, STATE> res0(3, 1), resid(3, 1), residguess(3, 1), diff(2, 1);
    TPZFNMatrix<9, STATE> jac(3, 2);
    TPZManVector<STATE> sigHWCyl(3), sigCart(3);
    F1Cyl(xi, beta, sigHWCyl);
    TPZHWTools::FromHWCylToPrincipal(sigHWCyl, sigCart);
    for (int i = 0; i < 3; i++) {
        res0(i) = sigCart[i];
    }
    DF1Cart(xi, beta, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i = 1; i <= 10; i++) {
        STATE diffxi = deltaxi * i / 10.;
        STATE xinext = xi + diffxi;
        STATE diffbeta = deltabeta * i / 10.;
        diff(0) = diffxi;
        diff(1) = diffbeta;
        STATE betanext = beta + diffbeta;
        F1Cyl(xinext, betanext, sigHWCyl);
        TPZHWTools::FromHWCylToPrincipal(sigHWCyl, sigCart);
        for (int ii = 0; ii < 3; ii++) {
            resid(ii) = sigCart[ii];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
        xnorm[i - 1] = Norm(diff);
        errnorm[i - 1] = Norm(resid - residguess);
    }
}

void TPZSandlerExtended::TaylorCheckDF2Cart(STATE theta, STATE beta, STATE k, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const {
    STATE deltatheta = 0.4;
    STATE deltabeta = 0.05;
    STATE deltak = 0.5;
    TPZFNMatrix<3, STATE> res0(3, 1), resid(3, 1), residguess(3, 1), diff(3, 1);
    TPZFNMatrix<9, STATE> jac(3, 3);
    DF2Cart(theta, beta, k, jac);
    TPZManVector<STATE> sigHWCyl(3), sigCart(3);
    F2Cyl(theta, beta, k, sigHWCyl);
    TPZHWTools::FromHWCylToPrincipal(sigHWCyl, sigCart);
    for (int i = 0; i < 3; i++) {
        res0(i) = sigCart[i];
    }
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i = 1; i <= 10; i++) {
        STATE difftheta = deltatheta * i / 10.;
        STATE thetanext = theta + difftheta;
        STATE diffbeta = deltabeta * i / 10.;
        STATE betanext = beta + diffbeta;
        STATE diffk = deltak * i / 10.;
        STATE knext = k + diffk;
        diff(0) = difftheta;
        diff(1) = diffbeta;
        diff(2) = diffk;
        F2Cyl(thetanext, betanext, knext, sigHWCyl);
        TPZHWTools::FromHWCylToPrincipal(sigHWCyl, sigCart);
        for (int ii = 0; ii < 3; ii++) {
            resid(ii) = sigCart[ii];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
        xnorm[i - 1] = Norm(diff);
        errnorm[i - 1] = Norm(resid - residguess);
    }

}

void TPZSandlerExtended::TaylorCheckProjectSigma(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const {
    TPZManVector<STATE, 3> deltasigma(3, -0.000012), sigproj(3);
    deltasigma[1] = -4.e-6;
    deltasigma[2] = -4.e-6;
    //    TPZManVector<STATE,3> deltasigma(3,-0.01), sigproj(3);
    //    deltasigma[1] *= 3;
    //    deltasigma[2] *= -2;

    TPZFNMatrix<3, STATE> res0(3, 1), diff(3, 1), resid(3, 1), residguess(3, 1);
    TPZFNMatrix<9, STATE> jac(3, 3);
    STATE kproj;
    int m_type;
    ProjectSigmaDep(sigtrial, kprev, sigproj, kproj, jac);
    for (int j = 0; j < 3; j++) res0(j) = sigproj[j];
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i = 1; i <= 10; i++) {
        TPZManVector<STATE, 3> diffsigma(3), nextsigma(3);
        for (int j = 0; j < 3; j++) {
            diffsigma[j] = deltasigma[j] * i / 10.;
            nextsigma[j] = sigtrial[j] + diffsigma[j];
            diff(j) = diffsigma[j];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
        ProjectSigma(nextsigma, kprev, sigproj, kproj, m_type);
        for (int j = 0; j < 3; j++) resid(j) = sigproj[j];
        xnorm[i - 1] = Norm(diff);
        errnorm[i - 1] = Norm(resid - residguess);
    }
}

/// verifies the validity of dxi/dsigtrial and dbeta/dsigtrial

void TPZSandlerExtended::TaylorCheckParamF1Sigtrial(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const {
    TPZManVector<STATE, 3> deltasigma(3, -0.01), sigproj(3);
    deltasigma[1] *= 3.;
    deltasigma[2] *= -2.;
    TPZFNMatrix<3, STATE> res0(2, 1), diff(3, 1), resid(2, 1), residguess(2, 1);
    TPZFNMatrix<9, STATE> jac(2, 3), jacF1(2, 2), gradF1(2, 3);
    STATE kproj;
    ProjectF1(sigtrial, kprev, sigproj, kproj);
    STATE xi, beta;
    SurfaceParamF1(sigproj, xi, beta);
    res0(0) = xi;
    res0(1) = beta;
    D2DistFunc1(sigtrial, xi, beta, jacF1);
    GradF1SigmaTrial(sigtrial, xi, beta, gradF1);
    //gradF1.Transpose();
    jacF1.Solve_LU(&gradF1);
    jac = -1. * gradF1;
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i = 1; i <= 10; i++) {
        TPZManVector<STATE, 3> diffsigma(3), nextsigma(3);
        for (int j = 0; j < 3; j++) {
            diffsigma[j] = deltasigma[j] * i / 10.;
            nextsigma[j] = sigtrial[j] + diffsigma[j];
            diff(j) = diffsigma[j];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
        ProjectF1(nextsigma, kprev, sigproj, kproj);
        SurfaceParamF1(sigproj, xi, beta);
        resid(0) = xi;
        resid(1) = beta;
        xnorm[i - 1] = Norm(diff);
        errnorm[i - 1] = Norm(resid - residguess);
    }

}

void TPZSandlerExtended::TaylorCheckProjectF1(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const {
    TPZManVector<STATE, 3> deltasigma(3, -0.01), sigproj(3);
    deltasigma[1] *= 3.;
    deltasigma[2] *= -2.;
    TPZFNMatrix<3, STATE> res0(3, 1), diff(3, 1), resid(3, 1), residguess(3, 1);
    TPZFNMatrix<9, STATE> jac(3, 3), jacF1(2, 2), gradF1(2, 3), DF1cart(2, 3), GradSigma(3, 3);
    STATE kproj;
    ProjectF1(sigtrial, kprev, sigproj, kproj);
    STATE xi, beta;
    SurfaceParamF1(sigproj, xi, beta);
    res0(0) = sigproj[0];
    ;
    res0(1) = sigproj[1];
    res0(2) = sigproj[2];
    D2DistFunc1(sigtrial, xi, beta, jacF1);
    GradF1SigmaTrial(sigtrial, xi, beta, gradF1);
    //gradF1.Transpose();
    jacF1.Solve_LU(&gradF1);
    DF1Cart(xi, beta, DF1cart);
    DF1cart.Multiply(gradF1, GradSigma);
    GradSigma *= -1.;

    jac = GradSigma;
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i = 1; i <= 10; i++) {
        TPZManVector<STATE, 3> diffsigma(3), nextsigma(3);
        for (int j = 0; j < 3; j++) {
            diffsigma[j] = deltasigma[j] * i / 10.;
            nextsigma[j] = sigtrial[j] + diffsigma[j];
            diff(j) = diffsigma[j];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
        ProjectF1(nextsigma, kprev, sigproj, kproj);
        resid(0) = sigproj[0];
        resid(1) = sigproj[1];
        resid(2) = sigproj[2];
        xnorm[i - 1] = Norm(diff);
        errnorm[i - 1] = Norm(resid - residguess);
    }

}

void TPZSandlerExtended::TaylorCheckProjectF2(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const {
    //    TPZManVector<STATE,3> deltasigma(3,-0.01), sigproj(3);
    //    deltasigma[1] *= 3.;
    //    deltasigma[2] *= -2.;
    TPZManVector<STATE, 3> deltasigma(3, -0.000012), sigproj(3);
    deltasigma[1] = -4.e-6;
    deltasigma[2] = -4.e-6;
    TPZFNMatrix<3, STATE> res0(3, 1), diff(3, 1), resid(3, 1), residguess(3, 1);
    TPZFNMatrix<9, STATE> jac(3, 3), jacF2(3, 3), gradF2(3, 3), DF2cart(3, 3), GradSigma(3, 3);
    STATE kproj;
    ProjectF2(sigtrial, kprev, sigproj, kproj);
    STATE theta, beta;
    SurfaceParamF2(sigproj, kproj, theta, beta);
    res0(0) = sigproj[0];
    res0(1) = sigproj[1];
    res0(2) = sigproj[2];
    Jacobianf2(sigtrial, theta, beta, kproj, jacF2);

    TFad<3, STATE> thetafad(theta, 0), betafad(beta, 1), kprojfad(kproj, 2);
    TPZManVector<TFad<3, STATE>, 3> sigtrialfad(3), ddistf2(3);
    for (int m = 0; m < 3; m++) {
        sigtrialfad[m] = sigtrial[m];
    }
    //DDistFunc2(sigtrialfad, thetafad, betafad, kprojfad, kprev, ddistf2);

    TPZFNMatrix<9, STATE> diffjac(3, 3);
    for (int m = 0; m < 3; m++) {
        for (int n = 0; n < 3; n++) {
            diffjac(m, n) = jacF2(m, n) - ddistf2[m].fastAccessDx(n);
            jacF2(m, n) = ddistf2[m].fastAccessDx(n);
        }
    }


    GradF2SigmaTrial(sigtrial, theta, beta, kproj, kprev, gradF2);
    //gradF1.Transpose();
    jacF2.Solve_LU(&gradF2);
    DF2Cart(theta, beta, kproj, DF2cart);
    DF2cart.Multiply(gradF2, GradSigma);
    GradSigma *= -1.;

    jac = GradSigma;
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i = 1; i <= 10; i++) {
        TPZManVector<STATE, 3> diffsigma(3), nextsigma(3);
        for (int j = 0; j < 3; j++) {
            diffsigma[j] = deltasigma[j] * i / 10.;
            nextsigma[j] = sigtrial[j] + diffsigma[j];
            diff(j) = diffsigma[j];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
        ProjectF2(nextsigma, kprev, sigproj, kproj);
        resid(0) = sigproj[0];
        resid(1) = sigproj[1];
        resid(2) = sigproj[2];
        xnorm[i - 1] = Norm(diff);
        errnorm[i - 1] = Norm(resid - residguess);
    }

}

void TPZSandlerExtended::TaylorCheckDtbkDsigtrial(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const {
    //    TPZManVector<STATE,3> deltasigma(3,-0.01), sigproj(3);
    //    deltasigma[1] *= 3.;
    //    deltasigma[2] *= -2.;
    TPZManVector<STATE, 3> deltasigma(3, -0.000012), sigproj(3);
    deltasigma[1] = -4.e-6;
    deltasigma[2] = -4.e-6;
    TPZFNMatrix<3, STATE> res0(3, 1), diff(3, 1), resid(3, 1), residguess(3, 1);
    TPZFNMatrix<9, STATE> jac(3, 3), jacF2(3, 3), gradF2(3, 3), GradTBK(3, 3);
    STATE kproj;
    ProjectF2(sigtrial, kprev, sigproj, kproj);
    STATE theta, beta;
    SurfaceParamF2(sigproj, kproj, theta, beta);
    res0(0) = theta;
    res0(1) = beta;
    res0(2) = kproj;
    Jacobianf2(sigtrial, theta, beta, kproj, jacF2);
    TFad<3, STATE> thetafad(theta, 0), betafad(beta, 1), kprojfad(kproj, 2);
    TPZManVector<TFad<3, STATE>, 3> sigtrialfad(3), ddistf2(3);
    for (int m = 0; m < 3; m++) {
        sigtrialfad[m] = sigtrial[m];
    }
    //DDistFunc2(sigtrialfad, thetafad, betafad, kprojfad, kprev, ddistf2);

    TPZFNMatrix<9, STATE> diffjac(3, 3);
    for (int m = 0; m < 3; m++) {
        for (int n = 0; n < 3; n++) {
            diffjac(m, n) = jacF2(m, n) - ddistf2[m].fastAccessDx(n);
            jacF2(m, n) = ddistf2[m].fastAccessDx(n);
        }
    }
    diffjac.Print("DiffMatrix");

    GradF2SigmaTrial(sigtrial, theta, beta, kproj, kprev, gradF2);
    //gradF1.Transpose();
    TPZFMatrix<STATE> jacF2temp(jacF2), gradF2temp(gradF2);
    jacF2temp.Solve_LU(&gradF2temp);
    GradTBK = gradF2temp;
    GradTBK *= -1.;

    jac = GradTBK;

    TPZFNMatrix<3, STATE> tbk(3, 1, 0.), rhs;
    tbk(1, 0) = 0.1;
    rhs = tbk;
    GradTBK.Solve_LU(&rhs);
    for (int i = 0; i < 3; i++) {
        deltasigma[i] = rhs(i, 0);
    }


    xnorm.resize(10);
    errnorm.resize(10);
    TPZFNMatrix<30, STATE> erros(3, 10);
    for (int i = 1; i <= 10; i++) {
        TPZManVector<STATE, 3> diffsigma(3), nextsigma(3);
        TPZFNMatrix<3, STATE> difftbk(3, 1);
        for (int j = 0; j < 3; j++) {
            diffsigma[j] = deltasigma[j] * i / 10.;
            difftbk(j) = tbk[j] * i / 10.;
            nextsigma[j] = sigtrial[j] + diffsigma[j];
            diff(j) = diffsigma[j];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
        //        TPZFNMatrix<3,STATE> multsig(3,0),multbk(3,0);
        //        jacF2.Multiply(difftbk, multbk);
        //        gradF2.Multiply(diff, multsig);
        //        residguess = multbk;
        //        residguess = multsig;
        //        TPZFNMatrix<3,STATE> fxn(3,1);
        //        DDistFunc2(nextsigma, theta+difftbk[0],beta+difftbk[1],kproj+difftbk[2],kprev,fxn);
        //        TPZManVector<STATE> fxnvec(3);
        //        DDistFunc2<STATE>(sigtrial,theta+difftbk[0],beta+difftbk[1],kproj+difftbk[2],kprev,fxnvec);
        //        for(int k=0; k<3; k++) fxn(k,0) = fxnvec[k];
        //        DDistFunc2(sigtrial, theta+difftbk[0],beta+difftbk[1],kproj+difftbk[2],kprev,fxn);
        //        DDistFunc2(nextsigma, theta,beta,kproj,kprev,fxn);
        ProjectF2(nextsigma, kprev, sigproj, kproj);
        STATE theta, beta;
        SurfaceParamF2(sigproj, kproj, theta, beta);
        resid(0) = theta;
        resid(1) = beta;
        resid(2) = kproj;
        xnorm[i - 1] = Norm(diff);
        errnorm[i - 1] = Norm(resid - residguess);
        for (int a = 0; a < 3; a++) erros(a, i - 1) = fabs(resid[a] - residguess[a]);
    }
    erros.Print(cout);

}

void TPZSandlerExtended::MCormicRanchSand(TPZSandlerExtended &mat)//em ksi
{
    STATE E = 100, nu = 0.25, A = 0.25, B = 0.67, C = 0.18, D = 0.67, R = 2.5, W = 0.066, N = 0., phi = 0, psi = 1.0;
    STATE G = E / (2. * (1. + nu));
    STATE K = E / (3. * (1. - 2 * nu));
    mat.fA = A;
    mat.fB = B;
    mat.fC = C;
    mat.fD = D;
    mat.fK = K;
    mat.fG = G;
    mat.fW = W;
    mat.fR = R;
    mat.fPhi = phi;
    mat.fN = N;
    mat.fPsi = psi;
    mat.fE = E;
    mat.fnu = nu;
    TPZElasticResponse ER;
    ER.SetUp(E, nu);
    mat.fElasticResponse = ER;

}

void TPZSandlerExtended::ReservoirSandstone(TPZSandlerExtended &mat)//em ksi
{
    STATE E = 1305, nu = 0.25, A = 2.61, B = 0.169, C = 2.57, D = 0.05069, R = 1.5, W = 0.0908, N = 0., phi = 0, psi = 1.0;
    STATE G = E / (2. * (1. + nu));
    STATE K = E / (3. * (1. - 2 * nu));
    mat.fA = A;
    mat.fB = B;
    mat.fC = C;
    mat.fD = D;
    mat.fK = K;
    mat.fG = G;
    mat.fW = W;
    mat.fR = R;
    mat.fPhi = phi;
    mat.fN = N;
    mat.fPsi = psi;
    mat.fE = E;
    mat.fnu = nu;
    TPZElasticResponse ER;
    ER.SetUp(E, nu);
    mat.fElasticResponse = ER;


}

void TPZSandlerExtended::SalemLimestone(TPZSandlerExtended &mat)// em MPa
{
    STATE E = 22547., nu = 0.2524, A = 689.2,
            B = 3.94e-4, C = 675.2, D = 1.47e-3, R = 28, W = 0.08, N = 6., phi = 0, psi = 1.0;
    STATE G = E / (2. * (1. + nu));
    STATE K = E / (3. * (1. - 2 * nu));
    mat.fA = A;
    mat.fB = B;
    mat.fC = C;
    mat.fD = D;
    mat.fK = K;
    mat.fG = G;
    mat.fW = W;
    mat.fR = R;
    mat.fPhi = phi;
    mat.fN = N;
    mat.fPsi = psi;
    mat.fE = E;
    mat.fnu = nu;
    TPZElasticResponse ER;
    ER.SetUp(E, nu);
    mat.fElasticResponse = ER;
}

void TPZSandlerExtended::PreSMat(TPZSandlerExtended &mat)// em MPa
{

    STATE E = 29269, nu = 0.203, A = 116.67,
            B = 0.0036895, C = 111.48, D = 0.018768, R = 0.91969, W = 0.006605, N = 0., phi = 0, psi = 1.0;
    STATE G = E / (2. * (1. + nu));
    STATE K = E / (3. * (1. - 2 * nu));
    mat.fA = A;
    mat.fB = B;
    mat.fC = C;
    mat.fD = D;
    mat.fK = K;
    mat.fG = G;
    mat.fW = W;
    mat.fR = R;
    mat.fPhi = phi;
    mat.fN = N;
    mat.fPsi = psi;
    mat.fE = E;
    mat.fnu = nu;
    TPZElasticResponse ER;
    ER.SetUp(E, nu);
    mat.fElasticResponse = ER;
}

int TPZSandlerExtended::ClassId() const {
    return Hash("TPZSandlerExtended");
}
