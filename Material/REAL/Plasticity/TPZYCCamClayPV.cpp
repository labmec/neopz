/* 
 * File:   TPZYCCamClayPV.cpp
 * Author: quinelato
 * 
 * Created on August 14, 2017, 5:18 PM
 */

#include "pzerror.h"
#include "TPZYCCamClayPV.h"
#include "TPZHWTools.h"

using namespace std;

TPZYCCamClayPV::TPZYCCamClayPV() : fA(0.), fGamma(0.), fM(0.), fPt(0.) {
}

TPZYCCamClayPV::TPZYCCamClayPV(const TPZYCCamClayPV& other) : fER(other.fER), fA(other.fA), fGamma(other.fGamma), fM(other.fM), fPt(other.fPt) {
}

void TPZYCCamClayPV::SetUp(const TPZElasticResponse &ER, REAL a, REAL gamma, REAL m, REAL pt) {
    SetElasticResponse(ER);
    fA = a;
    fGamma = gamma;
    fM = m;
    fPt = pt;
}

void TPZYCCamClayPV::SetElasticResponse(const TPZElasticResponse &ER) {
    fER = ER;
}

void TPZYCCamClayPV::Read(TPZStream &buf) {
    buf.Read(&fA);
    buf.Read(&fGamma);
    buf.Read(&fM);
    buf.Read(&fPt);
    fER.Read(buf);
}

void TPZYCCamClayPV::Write(TPZStream &buf) const {
    buf.Write(&fA);
    buf.Write(&fGamma);
    buf.Write(&fM);
    buf.Write(&fPt);
    fER.Write(buf);
}

REAL TPZYCCamClayPV::bFromP(REAL p) const {
    return (p >= fPt - fA) ? 1. : fGamma;
}

REAL TPZYCCamClayPV::bFromTheta(REAL theta) const {
    return (theta >= M_PI_2) ? 1. : fGamma;
}

void TPZYCCamClayPV::Phi(TPZVec<REAL> sigmaPV, REAL alpha, TPZVec<REAL> &phi)const {
    phi.resize(NYield);

    TPZTensor<REAL> sigma;
    sigma.XX() = sigmaPV[0];
    sigma.YY() = sigmaPV[1];
    sigma.ZZ() = sigmaPV[2];

    REAL p = sigma.I1() / 3.;
    REAL q = sqrt(3. * sigma.J2());

    phi[0] = 1 / (pow(bFromP(p), 2.)) * pow(p - fPt + fA, 2.) + pow(q / fM, 2.) - pow(fA, 2);
}

void TPZYCCamClayPV::SurfaceInCyl(const REAL theta, const REAL beta, TPZVec<REAL> &returnValue) const {
    const REAL SQRT1_3 = sqrt(1 / 3.);
    const REAL I1 = 3. * (fPt - fA * (1 + bFromTheta(theta) * cos(theta)));
    REAL sqrtj2 = fM * SQRT1_3 * fA * sin(theta);
    REAL xi = SQRT1_3*I1;
    REAL rho = M_SQRT2*sqrtj2;
    returnValue[0] = xi; // hydrostatic component
    returnValue[1] = rho; // deviatoric component
    returnValue[2] = beta; // Lode angle
}

void TPZYCCamClayPV::DDistanceToSurface(const TPZVec<STATE> &pt, const STATE theta, const STATE beta, TPZFMatrix<STATE> &fxn) const{

}

void TPZYCCamClayPV::D2DistanceToSurface(const TPZVec<STATE> &pt, const STATE theta, const STATE beta, TPZFMatrix<STATE> &jac) const{

}

REAL TPZYCCamClayPV::DistanceToSurface(const TPZVec<REAL> &sigma_trial_pv, const REAL theta, const REAL beta) const {
    TPZManVector<REAL, 3> projection_cyl(3);
    SurfaceInCyl(theta, beta, projection_cyl);
    TPZManVector<REAL, 3> projection_cart(3);
    TPZHWTools::FromHWCylToHWCart(projection_cyl, projection_cart);
    TPZManVector<REAL, 3> trial_cart(3);
    TPZHWTools::FromPrincipalToHWCart(sigma_trial_pv, trial_cart);
    REAL k = fER.K();
    REAL g = fER.G();
    return ((1. / (3. * k))*(trial_cart[0] - projection_cart[0])*(trial_cart[0] - projection_cart[0]))
            +(1. / (2. * g))*((trial_cart[1] - projection_cart[1])*(trial_cart[1] - projection_cart[1])+(trial_cart[2] - projection_cart[2])*(trial_cart[2] - projection_cart[2]));

}

void TPZYCCamClayPV::ProjectToSurface(const TPZVec<REAL> &sigma_trial_pv, const REAL kprev, TPZVec<REAL> &sigma_pv, REAL &kproj, const REAL tol) const {
    REAL theta = 0.;
    REAL theta_distance = std::numeric_limits<REAL>::max();
    TPZManVector<REAL> sigma_cyl;
    TPZHWTools::FromPrincipalToHWCyl(sigma_trial_pv, sigma_cyl);
    REAL beta = sigma_cyl[2];
    {
        const REAL initial_theta_guess = 0; // initial_xi_guess = fPt*sqrt3; // multiplying by sqrt converts from p to xi coordinates
        const REAL final_theta_guess = M_PI; // final_xi_guess (fPt - (1 + fGamma * fA)) * sqrt3;
        const unsigned int n_steps_theta = 40;
        const REAL theta_interval = (final_theta_guess - initial_theta_guess) / n_steps_theta;
        for (unsigned int i = 0; i < n_steps_theta; ++i) {
            REAL theta_guess = initial_theta_guess + i*theta_interval;
            REAL distance = DistanceToSurface(sigma_trial_pv, theta_guess, beta);
            if (fabs(distance) < fabs(theta_distance)) {
                theta = theta_guess;
                theta_distance = distance;
            }
        }
    }

    REAL residual_norm = std::numeric_limits<REAL>::max();
    TPZFNMatrix<4, STATE> xn1(2, 1, 0.), xn(2, 1, 0.), sol(2, 1, 0.), fxn(2, 1, 0.);
    xn(0, 0) = theta;
    xn(1, 0) = beta;
    for (unsigned int i = 0; i < 30; ++i) {
        TPZFNMatrix<4, STATE> jac(2, 2);
        D2DistanceToSurface(sigma_trial_pv, xn(0), xn(1), jac);
        DDistanceToSurface(sigma_trial_pv, xn(0), xn(1), fxn);
        sol = fxn;
        residual_norm = Norm(sol);

#ifdef LOG4CXX
        if (loggerConvTest->isDebugEnabled()) {
            std::stringstream outfile; //("convergencF1.txt");
            outfile << i << " " << log(residual_norm) << endl;
            //jac.Print(outfile);
            //outfile<< "\n xn " << " "<<fxnvec <<endl;
            //outfile<< "\n res " << " "<<fxnvec <<endl;
            LOGPZ_DEBUG(loggerConvTest, outfile.str());
        }
#endif

        jac.Solve_LU(&sol);
        xn1 = xn - sol;
        xn = xn1;
        if (residual_norm < tol) break;
    }

    TPZManVector<REAL, 3> sigprojcyl(3);
    SurfaceInCyl(xn[0], xn[1], sigprojcyl);

    TPZHWTools::FromHWCylToPrincipal(sigprojcyl, sigma_pv);

    STATE kguess = kprev;
    STATE resl = ResLF1(sigma_trial_pv, sigma_pv, kguess, kprev);
    while (resl < 0.) {
        kguess += 1.;
        resl = ResLF1(sigma_trial_pv, sigma_pv, kguess, kprev);
    }

    for (unsigned int count = 0; count < 30; ++count) {
        REAL dresl = DResLF1(sigma_trial_pv, sigma_pv, kguess, kprev);

#ifdef PZDEBUG
        if (IsZero(dresl)) {
            DebugStop();
        }
#endif
        kguess -= resl / dresl;
        resl = ResLF1(sigma_trial_pv, sigma_pv, kguess, kprev);
        if (fabs(resl) < tol) break;
    }

#ifdef PZDEBUG
    {
        if (count >= 30) {
            DebugStop();
        }
    }
#endif

    kproj = kguess;
}

void TPZYCCamClayPV::ProjectSigma(const TPZVec<REAL> &sigma_trial, const REAL kprev, TPZVec<REAL> &sigma, REAL &kproj) const {
    TPZVec<REAL> yield(NYield);

    this->Phi(sigma_trial, kprev, yield);

    if (yield[0] <= 0.) {
        sigma = sigma_trial;
        kproj = kprev;
    } else {
        ProjectToSurface(sigma_trial, kprev, sigma, kproj, 1.e-5);
    }
    DebugStop();
}

void TPZYCCamClayPV::ProjectSigmaDep(const TPZVec<REAL> &sigma_trial, const REAL kprev, TPZVec<REAL> &sigma, REAL &kproj, TPZFMatrix<REAL> &GradSigma) const {
    DebugStop();
}

TPZYCCamClayPV::~TPZYCCamClayPV() {
}

