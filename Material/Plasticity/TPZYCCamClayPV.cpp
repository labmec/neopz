/* 
 * File:   TPZYCCamClayPV.cpp
 * Author: quinelato
 * 
 * Created on August 14, 2017, 5:18 PM
 */

#include "pzlog.h"
#include "pzerror.h"
#include "TPZYCCamClayPV.h"
#include "TPZHWTools.h"


#ifdef PZ_LOG
static TPZLogger loggerConvTest("ConvTest");
#endif

using namespace std;

TPZYCCamClayPV::TPZYCCamClayPV() : fGamma(0.), fM(0.), fPt(0.), fLogHardening(0.), fLogBulkModulus(0.), fA0(0.), fE0(0.) {
}

TPZYCCamClayPV::TPZYCCamClayPV(const TPZYCCamClayPV& other) : fER(other.fER), fGamma(other.fGamma), fM(other.fM), fPt(other.fPt), fLogHardening(other.fLogHardening), fLogBulkModulus(other.fLogBulkModulus), fA0(other.fA0), fE0(other.fE0) {
}

void TPZYCCamClayPV::SetUp(const TPZElasticResponse &ER, REAL gamma, REAL m, REAL pt, REAL logHardening, REAL logBulkModulus, REAL a0, REAL e0) {
    SetElasticResponse(ER);
    fGamma = gamma;
    fM = m;
    fPt = pt;
    fLogHardening = logHardening;
    fLogBulkModulus = logBulkModulus;
    fA0 = a0;
    fE0 = e0;
}

void TPZYCCamClayPV::SetElasticResponse(const TPZElasticResponse &ER) {
    fER = ER;
}

int TPZYCCamClayPV::ClassId() const{
    return Hash("TPZYCCamClayPV");
}

void TPZYCCamClayPV::Read(TPZStream& buf, void* context) {
    buf.Read(&fGamma);
    buf.Read(&fM);
    buf.Read(&fPt);
    buf.Read(&fLogHardening);
    buf.Read(&fLogBulkModulus);
    buf.Read(&fA0);
    buf.Read(&fE0);
    fER.Read(buf, context);
}

void TPZYCCamClayPV::Write(TPZStream& buf, int withclassid) const {
    buf.Write(&fGamma);
    buf.Write(&fM);
    buf.Write(&fPt);
    buf.Write(&fLogHardening);
    buf.Write(&fLogBulkModulus);
    buf.Write(&fA0);
    buf.Write(&fE0);
    fER.Write(buf, withclassid);
}

REAL TPZYCCamClayPV::InitialDamage(const TPZVec<REAL> &stress_p) const{
    std::cout << "Method not implemented." << std::endl;
    DebugStop();
    return -1.0;
}

REAL TPZYCCamClayPV::bFromP(const REAL p, const REAL a) const {
    return (p >= fPt - a) ? 1. : fGamma;
}

REAL TPZYCCamClayPV::bFromTheta(REAL theta) const {
    return (theta >= M_PI_2) ? 1. : fGamma;
}

STATE TPZYCCamClayPV::PlasticVolumetricStrain(REAL a) const {
    STATE log_a_a0 = log(a / fA0);
    return log((1 + fE0 - fLogHardening * log_a_a0) / (1 + fE0 - fLogBulkModulus * log_a_a0));
}

void TPZYCCamClayPV::Phi(TPZVec<REAL> sigmaPV, REAL a, TPZVec<REAL> &phi)const {
    phi.resize(NYield);

    TPZTensor<REAL> sigma;
    sigma.XX() = sigmaPV[0];
    sigma.YY() = sigmaPV[1];
    sigma.ZZ() = sigmaPV[2];

    const REAL p = sigma.I1() / 3.;
    const REAL q = sqrt(3. * sigma.J2());

    phi[0] = 1 / (pow(bFromP(p, a), 2.)) * pow(p - fPt + a, 2.) + pow(q / fM, 2.) - pow(a, 2);
}

void TPZYCCamClayPV::SurfaceInCyl(const REAL theta, const REAL beta, const REAL a, TPZVec<REAL> &returnValue) const {
    const REAL SQRT1_3 = sqrt(1 / 3.);
    const REAL I1 = 3. * (fPt - a * (1 + bFromTheta(theta) * cos(theta)));
    REAL sqrtj2 = fM * SQRT1_3 * a * sin(theta);
    REAL xi = SQRT1_3*I1;
    REAL rho = M_SQRT2*sqrtj2;
    returnValue[0] = xi; // hydrostatic component
    returnValue[1] = rho; // deviatoric component
    returnValue[2] = beta; // Lode angle
}

REAL TPZYCCamClayPV::ResLFunc(const TPZVec<STATE> &sigma_trial_pv, STATE theta, STATE beta, REAL a, REAL aPrev) const {
    const STATE ctheta = cos(theta);
    const STATE b = bFromTheta(theta);
    const STATE I1Trial = (sigma_trial_pv[0])+(sigma_trial_pv[1])+(sigma_trial_pv[2]);
    const STATE I1Proj = 3. * (fPt - a * (1 + b * ctheta));
    return (I1Trial - I1Proj) - 3 * fER.K()*(PlasticVolumetricStrain(a) - PlasticVolumetricStrain(aPrev));
}

//REAL TPZYCCamClayPV::DResLFunc(const TPZVec<STATE> &sigma_trial_pv, STATE theta, STATE beta, REAL a, REAL aPrev) const {
//    const STATE ctheta = cos(theta);
//    const STATE b = bFromTheta(theta);
//    const STATE I1Trial = (sigma_trial_pv[0])+(sigma_trial_pv[1])+(sigma_trial_pv[2]);
//    const STATE I1Proj = 3. * (fPt - a * (1 + b * ctheta));
//    return (I1Trial - I1Proj) - 3 * fER.K()*(PlasticVolumetricStrain(a) - PlasticVolumetricStrain(aPrev));
//}

REAL TPZYCCamClayPV::DistanceToSurface(const TPZVec<REAL> &sigma_trial_pv, const REAL theta, const REAL beta, const REAL a) const {
    TPZManVector<REAL, 3> projection_cyl(3);
    SurfaceInCyl(theta, beta, a, projection_cyl);
    TPZManVector<REAL, 3> projection_cart(3);
    TPZHWTools::FromHWCylToHWCart(projection_cyl, projection_cart);
    TPZManVector<REAL, 3> trial_cart(3);
    TPZHWTools::FromPrincipalToHWCart(sigma_trial_pv, trial_cart);
    REAL k = fER.K();
    REAL g = fER.G();
    return ((1. / (3. * k))*(trial_cart[0] - projection_cart[0])*(trial_cart[0] - projection_cart[0]))
            +(1. / (2. * g))*((trial_cart[1] - projection_cart[1])*(trial_cart[1] - projection_cart[1])+(trial_cart[2] - projection_cart[2])*(trial_cart[2] - projection_cart[2]));
}

void TPZYCCamClayPV::DDistanceToSurface(const TPZVec<STATE> &sigma_trial_pv, const STATE theta, const STATE beta, const REAL a, const REAL aPrev, TPZVec<STATE> &fxn) const {
    const STATE stheta = sin(theta);
    const STATE ctheta = cos(theta);
    const STATE sbeta = sin(beta);
    const STATE cbeta = cos(beta);
    const REAL sqrt2_sqrt3 = sqrt(2. / 3.);
    const REAL sqrt3_3 = sqrt(3.) / 3;
    const STATE b = bFromTheta(theta);
    TPZVec<STATE> ptcart(3);
    TPZHWTools::FromPrincipalToHWCart(sigma_trial_pv, ptcart);
    const STATE sig1 = ptcart[0];
    const STATE sig2 = ptcart[1];
    const STATE sig3 = ptcart[2];
    fxn.Resize(3);
    fxn[0] = 2. * a * b * stheta / fER.K()*(-sqrt3_3 * sig1 + fPt - a * (1 + b * ctheta)) + a * fM * sqrt2_sqrt3 * ctheta * (-sig2 * cbeta - sig3 * sbeta + a * fM * sqrt2_sqrt3 * stheta) / fER.G();
    fxn[1] = a * fM * sqrt2_sqrt3 * stheta * (sbeta * sig2 - cbeta * sig3) / fER.G();
    fxn[2] = ResLFunc(sigma_trial_pv, theta, beta, a, aPrev);
}

void TPZYCCamClayPV::D2DistanceToSurface(const TPZVec<STATE> &sigma_trial_pv, const STATE theta, const STATE beta, const REAL a, TPZFNMatrix<9, STATE> &jac) const {
    const STATE b = bFromTheta(theta);
    const STATE a2 = pow(a, 2.);
    const STATE b2 = pow(b, 2.);
    const STATE M2 = pow(fM, 2.);
    const STATE stheta = sin(theta);
    const STATE stheta2 = pow(stheta, 2.);
    const STATE ctheta = cos(theta);
    const STATE ctheta2 = pow(ctheta, 2);
    const STATE sbeta = sin(beta);
    const STATE cbeta = cos(beta);
    const REAL sqrt2_sqrt3 = sqrt(2. / 3.);
    const REAL sqrt3_3 = sqrt(3.) / 3;
    TPZVec<STATE> ptcart(3);
    TPZHWTools::FromPrincipalToHWCart(sigma_trial_pv, ptcart);
    const STATE sig1 = ptcart[0];
    const STATE sig2 = ptcart[1];
    const STATE sig3 = ptcart[2];
    const STATE loga_a0 = log(a / fA0);

    jac.Resize(3, 3);
    jac(0, 0) = 2. * a * b * ctheta * (-sqrt3_3 * sig1 + fPt - a * (1 + b * ctheta)) / fER.K() + 2 * a2 * b2 * stheta2 / fER.K() - a * fM * sqrt2_sqrt3 * stheta * (-sig2 * cbeta - sig3 * sbeta + a * fM * sqrt2_sqrt3 * stheta) / fER.G() + 2. * a2 * M2 * ctheta2 / (3. * fER.G());
    jac(0, 1) = a * fM * sqrt2_sqrt3 * ctheta * (sbeta * sig2 - cbeta * sig3) / fER.G();
    jac(0, 2) = 2. * b * stheta * (-sqrt3_3 * sig1 + fPt - 2. * a * (1. + b * ctheta)) / fER.K() + fM * sqrt2_sqrt3 * ctheta * (-sig2 * cbeta - sig3 * sbeta) / fER.G() + 4. * a * M2 * stheta * ctheta / (3. * fER.G());
    jac(1, 0) = jac(0, 1);
    jac(1, 1) = a * fM * sqrt2_sqrt3 * stheta * (sbeta * sig2 + cbeta * sig3) / fER.G();
    jac(1, 2) = fM * sqrt2_sqrt3 * stheta * (sbeta * sig2 - cbeta * sig3) / fER.G();
    jac(2, 0) = -3. * a * b*stheta;
    jac(2, 1) = 0.;
    jac(2, 2) = 3. * (1 + b * ctheta + fER.K() * fLogHardening / (a * (1. + fE0 - fLogHardening * loga_a0)) - fER.K() * fLogBulkModulus / (a * (1. + fE0 - fLogBulkModulus * loga_a0)));
}

void TPZYCCamClayPV::ProjectToSurfaceConstantBeta(const TPZVec<REAL> &sigma_trial_pv, const REAL aPrev, TPZVec<REAL> &sigma_pv, REAL &aProj, const REAL tol) const {
    REAL theta = 0.;
    REAL theta_distance = std::numeric_limits<REAL>::max();
    REAL beta = 0;
    {
        const REAL initial_theta_guess = 0; // initial_xi_guess = fPt*sqrt3; // multiplying by sqrt converts from p to xi coordinates
        const REAL final_theta_guess = M_PI; // final_xi_guess (fPt - (1 + fGamma * fA)) * sqrt3;
        const unsigned int n_steps_theta = 40;
        const REAL theta_interval = (final_theta_guess - initial_theta_guess) / n_steps_theta;
        for (unsigned int i = 0; i < n_steps_theta; ++i) {
            REAL theta_guess = initial_theta_guess + i*theta_interval;
            REAL distance = DistanceToSurface(sigma_trial_pv, theta_guess, beta, aPrev);
            if (fabs(distance) < fabs(theta_distance)) {
                theta = theta_guess;
                theta_distance = distance;
            }
        }
    }

    REAL residual_norm = std::numeric_limits<REAL>::max();
    TPZFNMatrix<3, STATE> xn(3, 1, 0.), sol(3, 1, 0.);
    TPZManVector<STATE> fxn(3);
    xn(0, 0) = theta;
    xn(1, 0) = beta;
    xn(2, 0) = aPrev;
    for (unsigned int i = 0; i < 30; ++i) {
        TPZFNMatrix<9, STATE> jac(3, 3);
        D2DistanceToSurface(sigma_trial_pv, xn(0), xn(1), xn(2), jac);
        DDistanceToSurface(sigma_trial_pv, xn(0), xn(1), xn(2), aPrev, fxn);
        for (unsigned int k = 0; k < 3; ++k) {
            sol(k, 0) = fxn[k];
        }
        residual_norm = Norm(sol);
        
        for (unsigned int k = 0; k < 3; k++) {
            jac(k, 1) = 0.;
            jac(1, k) = 0.;
        }
        jac(1, 1) = 1.;
        
#ifdef PZ_LOG
        if (loggerConvTest.isDebugEnabled()) {
            std::stringstream outfile; //("convergencF1.txt");
            outfile << i << " " << log(residual_norm) << endl;
            //jac.Print(outfile);
            //outfile<< "\n xn " << " "<<fxnvec <<endl;
            //outfile<< "\n res " << " "<<fxnvec <<endl;
            LOGPZ_DEBUG(loggerConvTest, outfile.str());
        }
#endif

        jac.Solve_LU(&sol);
        xn(0) = xn(0) - sol(0);
        xn(1) = beta;
        xn(2) = xn(2) - sol(2);
        if (residual_norm < tol) break;
    }

    STATE thetasol, betasol, asol;

    thetasol = xn(0);
    betasol = xn(1);
    asol = xn(2);
    aProj = asol;

    TPZManVector<STATE, 3> surfaceCyl(3);
    SurfaceInCyl(thetasol, betasol, asol, surfaceCyl);
    TPZHWTools::FromHWCylToPrincipal(surfaceCyl, sigma_pv);
}

void TPZYCCamClayPV::ProjectToSurface(const TPZVec<REAL> &sigma_trial_pv, const REAL aPrev, TPZVec<REAL> &sigma_pv, REAL &aProj, const REAL tol) const {
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
            REAL distance = DistanceToSurface(sigma_trial_pv, theta_guess, beta, aPrev);
            if (fabs(distance) < fabs(theta_distance)) {
                theta = theta_guess;
                theta_distance = distance;
            }
        }
    }

    REAL residual_norm = std::numeric_limits<REAL>::max();
    TPZFNMatrix<3, STATE> xn(3, 1, 0.), sol(3, 1, 0.);
    TPZManVector<STATE> fxn(3);
    xn(0, 0) = theta;
    xn(1, 0) = beta;
    xn(2, 0) = aPrev;
    for (unsigned int i = 0; i < 30; ++i) {
        TPZFNMatrix<9, STATE> jac(3, 3);
        D2DistanceToSurface(sigma_trial_pv, xn(0), xn(1), xn(2), jac);
        DDistanceToSurface(sigma_trial_pv, xn(0), xn(1), xn(2), aPrev, fxn);
        for (unsigned int k = 0; k < 3; ++k) {
            sol(k, 0) = fxn[k];
        }
        residual_norm = Norm(sol);

#ifdef PZ_LOG
        if (loggerConvTest.isDebugEnabled()) {
            std::stringstream outfile; //("convergencF1.txt");
            outfile << i << " " << log(residual_norm) << endl;
            //jac.Print(outfile);
            //outfile<< "\n xn " << " "<<fxnvec <<endl;
            //outfile<< "\n res " << " "<<fxnvec <<endl;
            LOGPZ_DEBUG(loggerConvTest, outfile.str());
        }
#endif

        jac.Solve_LU(&sol);
        xn = xn - sol;
        if (residual_norm < tol) break;
    }

    STATE thetasol, betasol, asol;

    thetasol = xn(0);
    betasol = xn(1);
    asol = xn(2);
    aProj = asol;

    TPZManVector<STATE, 3> surfaceCyl(3);
    SurfaceInCyl(thetasol, betasol, asol, surfaceCyl);
    TPZHWTools::FromHWCylToPrincipal(surfaceCyl, sigma_pv);
}

void TPZYCCamClayPV::ProjectSigma(const TPZVec<REAL> &sigma_trial_pv, const REAL aPrev, TPZVec<REAL> &sigma_pv, REAL &aProj, int &m_type,  TPZFMatrix<REAL> * gradient) const {
    
    bool require_gradient_Q = true;
    if (!gradient) {
        require_gradient_Q = false;
    }
    
#ifdef PZDEBUG
    // Check for required dimensions of tangent
    if (!(gradient->Rows() == 3 && gradient->Cols() == 3)) {
        std::cerr << "Unable to compute the gradient operator. Required gradient array dimensions are 3x3." << std::endl;
        DebugStop();
    }
    
    if (require_gradient_Q) {
        DebugStop(); // implemented this functionality.
    }
    
#endif
    
    TPZVec<REAL> yield(NYield);
    this->Phi(sigma_trial_pv, aPrev, yield);

    if (yield[0] <= 0.) {
        sigma_pv = sigma_trial_pv;
        aProj = aPrev;
    } else {
        ProjectToSurface(sigma_trial_pv, aPrev, sigma_pv, aProj, 1.e-5);
    }
}

void TPZYCCamClayPV::SurfaceParam(const TPZVec<STATE> &sigma_pv, const STATE a, STATE &theta, STATE &beta) const {
    TPZManVector<STATE> sigmaHWCyl(3);
    TPZHWTools::FromPrincipalToHWCyl(sigma_pv, sigmaHWCyl);
    STATE I1 = sigmaHWCyl[0] * sqrt(3.);
    const REAL p = I1 / 3.;
    const STATE b = bFromP(p, a);
    STATE costheta = (-I1 / 3. + fPt - a) / (a * b);
    STATE sqrtj2 = M_SQRT1_2 * sigmaHWCyl[1];
    STATE sintheta = sqrt(3.) * sqrtj2 / (a * fM);
    theta = atan2(sintheta, costheta);
    beta = sigmaHWCyl[2];
}

void TPZYCCamClayPV::GradSigmaTrial(const TPZVec<REAL> &sigma_trial_pv, const REAL theta, const REAL beta, const REAL a, TPZFNMatrix<9, STATE> &ddist_dsigmatrial) const {
    const STATE stheta = sin(theta);
    const STATE ctheta = cos(theta);
    const STATE sbeta = sin(beta);
    const STATE cbeta = cos(beta);
    const REAL b = bFromTheta(theta);
    const REAL sqrt3 = sqrt(3.);
    const REAL sqrt27 = pow(3, 1.5);

    ddist_dsigmatrial(0, 0) = -(2 * a * b * stheta * fER.G() + 2 * a * cbeta * ctheta * fER.K() * fM) / (3 * fER.G() * fER.K());
    ddist_dsigmatrial(0, 1) = -(2 * sqrt3 * a * b * stheta * fER.G()+(3 * a * sbeta - sqrt3 * a * cbeta) * ctheta * fER.K() * fM) / (sqrt27 * fER.G() * fER.K());
    ddist_dsigmatrial(0, 2) = ((3 * a * sbeta + sqrt3 * a * cbeta) * ctheta * fER.K() * fM - 2 * sqrt3 * a * b * stheta * fER.G()) / (sqrt27 * fER.G() * fER.K());

    ddist_dsigmatrial(1, 0) = (2 * a * sbeta * stheta * fM) / (3 * fER.G());
    ddist_dsigmatrial(1, 1) = -((3 * a * cbeta + sqrt3 * a * sbeta) * stheta * fM) / (sqrt27 * fER.G());
    ddist_dsigmatrial(1, 2) = -((sqrt3 * a * sbeta - 3 * a * cbeta) * stheta * fM) / (sqrt27 * fER.G());

    ddist_dsigmatrial(2, 0) = 1.;
    ddist_dsigmatrial(2, 1) = 1.;
    ddist_dsigmatrial(2, 2) = 1.;
}

void TPZYCCamClayPV::DFuncCart(STATE theta, STATE beta, STATE a, TPZFNMatrix<9, STATE> &DFunccart) const {
    const STATE stheta = sin(theta);
    const STATE ctheta = cos(theta);
    const STATE sbeta = sin(beta);
    const STATE cbeta = cos(beta);
    const REAL b = bFromTheta(theta);
    const REAL sqrt3 = sqrt(3.);

    DFunccart(0, 0) = a * b * stheta + 2. * a * cbeta * ctheta * fM / 3.;
    DFunccart(0, 1) = -2 * a * sbeta * stheta * fM / 3.;
    DFunccart(0, 2) = -1 - b * ctheta + 2. * cbeta * stheta * fM / 3.;

    DFunccart(1, 0) = a * b * stheta + (sqrt3 * a * sbeta - a * cbeta) * ctheta * fM / 3.;
    DFunccart(1, 1) = ((sqrt3 * a * cbeta + a * sbeta) * stheta * fM) / 3.;
    DFunccart(1, 2) = -1 - b * ctheta + (sqrt3 * sbeta - cbeta) * stheta * fM / 3.;

    DFunccart(2, 0) = -(sqrt3 * a * sbeta + a * cbeta) * ctheta * fM / 3. - a * b*stheta;
    DFunccart(2, 1) = ((a * sbeta - sqrt3 * a * cbeta) * stheta * fM) / 3.;
    DFunccart(2, 2) = -1 + b * ctheta + (sqrt3 * sbeta + cbeta) * stheta * fM / 3.;
}

void TPZYCCamClayPV::ProjectSigmaDep(const TPZVec<REAL> &sigma_trial_pv, const REAL aPrev, TPZVec<REAL> &sigma, REAL &aProj, TPZFMatrix<REAL> &GradSigma) const {
    TPZVec<REAL> yield(NYield);
    this->Phi(sigma_trial_pv, aPrev, yield);

    if (yield[0] <= 0.) {
        sigma = sigma_trial_pv;
        aProj = aPrev;
        GradSigma.Identity();
    } else {
        const REAL tol=1.e-5;
        bool threeEigEqual = (fabs(sigma_trial_pv[0] - sigma_trial_pv[1]) < tol && fabs(sigma_trial_pv[1] - sigma_trial_pv[2]) < tol);
//        if (threeEigEqual){
//            ProjectToSurfaceConstantBeta(sigma_trial_pv, aPrev, sigma, aProj, 1.e-5);
//            // we can compute the tangent matrix
//            STATE theta, beta;
//            SurfaceParam(sigma, aProj, theta, beta);
//            TPZFNMatrix<9, STATE> ddist_dsigmatrial(3, 3), jac(3, 3), DFunccart(3, 3);
//            GradSigmaTrial(sigma_trial_pv, theta, beta, aProj, ddist_dsigmatrial);
//            D2DistanceToSurface(sigma_trial_pv, theta, beta, aProj, jac);
//            jac.Solve_LU(&ddist_dsigmatrial);
//            DFuncCart(theta, beta, aProj, DFunccart);
//            DFunccart.Multiply(ddist_dsigmatrial, GradSigma);
//            GradSigma *= -1.;
//        } else {
            ProjectToSurface(sigma_trial_pv, aPrev, sigma, aProj, 1.e-5);
            // we can compute the tangent matrix
            STATE theta, beta;
            SurfaceParam(sigma, aProj, theta, beta);
            TPZFNMatrix<9, STATE> ddist_dsigmatrial(3, 3), jac(3, 3), DFunccart(3, 3);
            GradSigmaTrial(sigma_trial_pv, theta, beta, aProj, ddist_dsigmatrial);
            D2DistanceToSurface(sigma_trial_pv, theta, beta, aProj, jac);
            jac.Solve_LU(&ddist_dsigmatrial);
            DFuncCart(theta, beta, aProj, DFunccart);
            DFunccart.Multiply(ddist_dsigmatrial, GradSigma);
            GradSigma *= -1.;
//        }
    }
}

TPZYCCamClayPV::~TPZYCCamClayPV() {
}

