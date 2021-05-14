/* 
 * File:   TPZYCDruckerPragerPV.cpp
 * Author: quinelato
 * 
 * Created on September 12, 2017, 2:24 PM
 */

#include "pzlog.h"
#include "TPZYCDruckerPragerPV.h"
#include "TPZHWTools.h"

#ifdef PZ_LOG
static TPZLogger loggerConvTest("ConvTest");
#endif

TPZYCDruckerPragerPV::TPZYCDruckerPragerPV() : fER(fCap.fER), fM(fCap.fM), fPt(fCap.fPt), fLogHardening(fCap.fLogHardening), fLogBulkModulus(fCap.fLogBulkModulus), fA0(fCap.fA0), fE0(fCap.fE0) {
}

TPZYCDruckerPragerPV::TPZYCDruckerPragerPV(const TPZYCDruckerPragerPV& other) : fCap(other.fCap), fER(fCap.fER), fM(fCap.fM), fPt(fCap.fPt), fLogHardening(fCap.fLogHardening), fLogBulkModulus(fCap.fLogBulkModulus), fA0(fCap.fA0), fE0(fCap.fE0) {
}

TPZYCDruckerPragerPV& TPZYCDruckerPragerPV::operator=(const TPZYCDruckerPragerPV& other) {
    fCap = other.fCap;
    return *this;
}

void TPZYCDruckerPragerPV::SetUp(const TPZElasticResponse &ER, REAL gamma, REAL m, REAL pt, REAL logHardening, REAL logBulkModulus, REAL a0, REAL e0) {
    fCap.SetUp(ER, gamma, m, pt, logHardening, logBulkModulus, a0, e0);
}

void TPZYCDruckerPragerPV::SetElasticResponse(const TPZElasticResponse &ER) {
    fER = ER;
}

int TPZYCDruckerPragerPV::ClassId() const {
    return Hash("TPZYCDruckerPragerPV");
}

void TPZYCDruckerPragerPV::Read(TPZStream& buf, void* context) {
    fCap.Read(buf, context);
}

void TPZYCDruckerPragerPV::Write(TPZStream& buf, int withclassid) const {
    fCap.Write(buf, withclassid);
}

REAL TPZYCDruckerPragerPV::InitialDamage(const TPZVec<REAL> &stress_p) const{
    std::cout << "Method not implemented." << std::endl;
    DebugStop();
    return -1.0;
}

REAL TPZYCDruckerPragerPV::bFromP(const REAL p, const REAL a) const {
    return fCap.bFromP(p, a);
}

REAL TPZYCDruckerPragerPV::bFromTheta(REAL theta) const {
    return fCap.bFromTheta(theta);
}

void TPZYCDruckerPragerPV::Phi(TPZVec<REAL> sigmaPV, REAL a, TPZVec<REAL> &phi) const {
    phi.resize(NYield);

    TPZTensor<REAL> sigma;
    sigma.XX() = sigmaPV[0];
    sigma.YY() = sigmaPV[1];
    sigma.ZZ() = sigmaPV[2];

    const REAL p = sigma.I1() / 3.;
    const REAL sqrt3 = sqrt(3);

    TPZVec<REAL> phiCap;
    fCap.Phi(sigmaPV, a, phiCap);

    phi[0] = sqrt(sigma.J2()) + fM * (p - fPt) / sqrt3;
    phi[1] = phiCap[0];
}

void TPZYCDruckerPragerPV::SurfaceF1InCyl(const REAL xi, const REAL beta, TPZVec<REAL> &returnValue) const {
    const REAL sqrt3 = sqrt(3);
    REAL sqrtj2 = sqrt3 * abs(fM * (xi / sqrt3 - fPt)) / 3.;
    REAL rho = M_SQRT2*sqrtj2;
    returnValue[0] = xi; // hydrostatic component
    returnValue[1] = rho; // deviatoric component
    returnValue[2] = beta; // Lode angle
}

void TPZYCDruckerPragerPV::SurfaceF2InCyl(const REAL theta, const REAL beta, const REAL a, TPZVec<REAL> &returnValue) const {
    fCap.SurfaceInCyl(theta, beta, a, returnValue);
}

REAL TPZYCDruckerPragerPV::ResLF1(const TPZVec<STATE> &sigma_trial_pv, const TPZVec<STATE> &sigma_proj_pv, const STATE a, const STATE aPrev) const {
    STATE I1Trial = (sigma_trial_pv[0] + sigma_trial_pv[1] + sigma_trial_pv[2]);
    STATE I1Proj = (sigma_proj_pv[0] + sigma_proj_pv[1] + sigma_proj_pv[2]);
    return (I1Trial - I1Proj) - 3 * fER.K()*(PlasticVolumetricStrain(a) - PlasticVolumetricStrain(aPrev));
}

REAL TPZYCDruckerPragerPV::ResLF2(const TPZVec<STATE> &sigma_trial_pv, const STATE theta, const STATE beta, const REAL a, const REAL aPrev) const {
    return fCap.ResLFunc(sigma_trial_pv, theta, beta, a, aPrev);
}

STATE TPZYCDruckerPragerPV::DResLF1(const TPZVec<STATE> &sigma_trial_pv, const TPZVec<STATE> &sigma_proj_pv, const STATE a, const STATE aPrev) const {
    STATE log_a = log(a);
    STATE log_a0 = log(fA0);
    STATE log_a0_log_a = log_a0 - log_a;

    return 3 * fER.K()*(fLogHardening - fLogBulkModulus) / (a * (log_a0_log_a * (1 + fE0) + pow(log_a0_log_a, 2.) * fLogBulkModulus) * fLogHardening + ((1 + fE0) * log_a0_log_a * fLogBulkModulus + pow(fE0 + 1, 2)));
}

REAL TPZYCDruckerPragerPV::DistanceToSurfaceF1(const TPZVec<REAL> &sigma_trial_pv, const REAL xi, const REAL beta) const {
    TPZManVector<STATE, 3> cyl(3);
    SurfaceF1InCyl(xi, beta, cyl);
    TPZManVector<STATE, 3> cart(3);
    TPZHWTools::FromHWCylToHWCart(cyl, cart);
    TPZManVector<STATE, 3> carttrial(3);
    TPZHWTools::FromPrincipalToHWCart(sigma_trial_pv, carttrial);

    return ((1. / (3. * fER.K()))*(carttrial[0] - cart[0])*(carttrial[0] - cart[0]))
            +(1. / (2. * fER.G()))*((carttrial[1] - cart[1])*(carttrial[1] - cart[1])+(carttrial[2] - cart[2])*(carttrial[2] - cart[2]));
}

REAL TPZYCDruckerPragerPV::DistanceToSurfaceF2(const TPZVec<REAL> &sigma_trial_pv, const REAL theta, const REAL beta, const REAL a) const {

    return fCap.DistanceToSurface(sigma_trial_pv, theta, beta, a);
}

void TPZYCDruckerPragerPV::DDistanceToSurfaceF1(const TPZVec<STATE> &sigma_trial_pv, const STATE xi, const STATE beta, TPZVec<STATE> &fxn) const {

    const REAL sqrt3 = sqrt(3);
    const STATE sbeta = sin(beta);
    const STATE cbeta = cos(beta);
    const REAL K = fER.K();
    const REAL G = fER.G();
    const STATE sig1trial = sigma_trial_pv[0];
    const STATE sig2trial = sigma_trial_pv[1];
    const STATE sig3trial = sigma_trial_pv[2];
    fxn.Resize(2);
    fxn[0] = ((6 * xi - 6 * sig1trial) * G + 3 * M_SQRT2 * (cbeta * sig2trial + sbeta * sig3trial) * K * fM + (2 * xi - 2 * sqrt3 * fPt) * K * pow(fM, 2)) / (9 * G * K);
    fxn[1] = M_SQRT2 * ((sqrt3 * sbeta * fPt * sig2trial - sqrt3 * cbeta * fPt * sig3trial + (cbeta * sig3trial - sbeta * sig2trial) * xi) * fM) / (3 * G);
}

void TPZYCDruckerPragerPV::DDistanceToSurfaceF2(const TPZVec<STATE> &sigma_trial_pv, const STATE theta, const STATE beta, const REAL a, const REAL aPrev, TPZVec<STATE> &fxn) const {

    fCap.DDistanceToSurface(sigma_trial_pv, theta, beta, a, aPrev, fxn);
}

void TPZYCDruckerPragerPV::D2DistanceToSurfaceF1(const TPZVec<STATE> &sigma_trial_pv, const STATE xi, const STATE beta, TPZFNMatrix<4, STATE> &jac) const {

    const REAL sqrt3 = sqrt(3);
    const STATE sbeta = sin(beta);
    const STATE cbeta = cos(beta);
    const REAL K = fER.K();
    const REAL G = fER.G();
    const STATE sig1trial = sigma_trial_pv[0];
    const STATE sig2trial = sigma_trial_pv[1];
    const STATE sig3trial = sigma_trial_pv[2];

    jac.Resize(2, 2);
    jac(0, 0) = 2. / (3. * K) + 2. * pow(fM, 2.) / (9. * G);
    jac(0, 1) = fM * M_SQRT2 * (cbeta * sig3trial - sbeta * sig2trial) / (3. * G);

    jac(1, 0) = jac(0, 1);
    jac(1, 1) = M_SQRT2 * ((sqrt3 * cbeta * fPt * sig2trial + sqrt3 * sbeta * fPt * sig3trial - (sbeta * sig3trial + cbeta * sig2trial) * xi) * fM) / (3. * G);
}

void TPZYCDruckerPragerPV::D2DistanceToSurfaceF2(const TPZVec<STATE> &sigma_trial_pv, const STATE theta, const STATE beta, const REAL a, TPZFNMatrix<9, STATE> &jac) const {

    fCap.D2DistanceToSurface(sigma_trial_pv, theta, beta, a, jac);
}

void TPZYCDruckerPragerPV::ProjectToSurfaceF1(const TPZVec<REAL> &sigma_trial_pv, const REAL aPrev, TPZVec<REAL> &sigma, REAL &aProj, const REAL tol) const {
    const REAL sqrt3 = sqrt(3);
    REAL xi = 0.;
    REAL xi_distance = std::numeric_limits<REAL>::max();
    TPZManVector<STATE> sigma_trial_star;
    TPZHWTools::FromPrincipalToHWCart(sigma_trial_pv, sigma_trial_star);
    STATE beta = atan2(sigma_trial_star[2], sigma_trial_star[1]);
    {
        const REAL initial_xi_guess = fPt*sqrt3; // multiplying by sqrt converts from p to xi coordinates
        const REAL final_xi_guess = (fPt - (1 + fCap.fGamma) * aPrev) * sqrt3;
        const unsigned int n_steps_xi = 40;
        const REAL xi_interval = (initial_xi_guess - final_xi_guess) / n_steps_xi;
        for (unsigned int i = 0; i < n_steps_xi; ++i) {
            REAL xi_guess = initial_xi_guess + i*xi_interval;
            REAL distance = DistanceToSurfaceF1(sigma_trial_pv, xi_guess, beta);
            if (fabs(distance) < fabs(xi_distance)) {
                xi = xi_guess;
                xi_distance = distance;
            }
        }
    }

    REAL residual_norm = std::numeric_limits<REAL>::max();
    TPZFNMatrix<4, STATE> xn(2, 1, 0.), sol(2, 1, 0.);
    TPZManVector<STATE> fxn(2);
    xn(0, 0) = xi;
    xn(1, 0) = beta;
    for (unsigned int i = 0; i < 30; ++i) {
        TPZFNMatrix<4, STATE> jac(2, 2);
        D2DistanceToSurfaceF1(sigma_trial_pv, xn(0), xn(1), jac);
        DDistanceToSurfaceF1(sigma_trial_pv, xn(0), xn(1), fxn);
        for (unsigned int k = 0; k < 2; ++k) {
            sol(k, 0) = fxn[k];
        }
        residual_norm = Norm(sol);

#ifdef PZ_LOG
        if (loggerConvTest.isDebugEnabled()) {
            std::stringstream outfile; //("convergencF1.txt");
            outfile << i << " " << log(residual_norm) << '\n';
            //jac.Print(outfile);
            //outfile<< "\n xn " << " "<< fxn <<endl;
            //outfile<< "\n res " << " "<< residual_norm <<endl;
            LOGPZ_DEBUG(loggerConvTest, outfile.str());
        }
#endif

        jac.Solve_LU(&sol);
        xn = xn - sol;
        if (residual_norm < tol) break;
    }

    TPZManVector<STATE, 3> sigma_proj_cyl(3);
    SurfaceF1InCyl(xn[0], xn[1], sigma_proj_cyl);

    TPZHWTools::FromHWCylToPrincipal(sigma_proj_cyl, sigma);

    STATE aGuess = aPrev;
    STATE resl = ResLF1(sigma_trial_pv, sigma, aGuess, aPrev);
    while (resl < 0.) {
        aGuess += 1.;
        resl = ResLF1(sigma_trial_pv, sigma, aGuess, aPrev);
    }
    unsigned int count;
    for (count = 0; count < 30; ++count) {
        STATE dresl = DResLF1(sigma_trial_pv, sigma, aGuess, aPrev);

#ifdef PZDEBUG
        if (IsZero(dresl)) {
            DebugStop();
        }
#endif

        aGuess -= resl / dresl;
        resl = ResLF1(sigma_trial_pv, sigma, aGuess, aPrev);
        if (fabs(resl) < tol) break;
    }

#ifdef PZDEBUG
    {
        if (count >= 30) {

            DebugStop();
        }
    }
#endif

    aProj = aGuess;
}

void TPZYCDruckerPragerPV::ProjectToSurfaceF2(const TPZVec<REAL> &sigma_trial_pv, const REAL aPrev, TPZVec<REAL> &sigma, REAL &aProj, const REAL tol) const {

    fCap.ProjectToSurface(sigma_trial_pv, aPrev, sigma, aProj, tol);
}

void TPZYCDruckerPragerPV::ProjectSigma(const TPZVec<REAL> &sigma_trial_pv, const REAL aPrev, TPZVec<REAL> &sigma_pv, REAL &aProj, int &m_type, TPZFMatrix<REAL> * gradient) const {
    
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
    STATE I1 = sigma_trial_pv[0] + sigma_trial_pv[1] + sigma_trial_pv[2];

    if (I1 / 3. >= fPt - aPrev) {
        if (yield[0] > 0.) {
            ProjectToSurfaceF1(sigma_trial_pv, aPrev, sigma_pv, aProj, 1.e-5);
        } else {
            sigma_pv = sigma_trial_pv;
            aProj = aPrev;
        }
    } else {
        if (yield[1] > 0.) {
            ProjectToSurfaceF2(sigma_trial_pv, aPrev, sigma_pv, aProj, 1.e-5);
        } else {
            sigma_pv = sigma_trial_pv;
            aProj = aPrev;
        }
    }
}

void TPZYCDruckerPragerPV::SurfaceParamF1(const TPZVec<STATE> &sigma_pv, STATE &xi, STATE &beta, const REAL tol) const {
    TPZManVector<STATE> sigmaHWCyl(3);
    TPZHWTools::FromPrincipalToHWCyl(sigma_pv, sigmaHWCyl);
    xi = sigmaHWCyl[0];
    beta = sigmaHWCyl[2];
#ifdef PZDEBUG
    STATE dist = DistanceToSurfaceF1(sigma_pv, xi, beta);
    if (fabs(dist) > tol) {
        DebugStop();
    }
#endif
}

void TPZYCDruckerPragerPV::SurfaceParamF2(const TPZVec<STATE> &sigma_pv, const STATE a, STATE &theta, STATE &beta, const REAL tol) const {
    fCap.SurfaceParam(sigma_pv, a, theta, beta);
}

void TPZYCDruckerPragerPV::GradF1SigmaTrial(const TPZVec<REAL> &sigma_trial_pv, const REAL xi, const REAL beta, const REAL aProj, TPZFNMatrix<6, STATE> &ddist_dsigmatrial) const {
    const REAL sqrt3 = sqrt(3);
    const STATE sbeta = sin(beta);
    const STATE cbeta = cos(beta);
    const REAL G = fER.G();
    const REAL K = fER.K();

    ddist_dsigmatrial(0, 0) = -2. / (3. * K);
    ddist_dsigmatrial(0, 1) = M_SQRT2 * cbeta * fM / (3. * G);
    ddist_dsigmatrial(0, 2) = M_SQRT2 * sbeta * fM / (3. * G);

    ddist_dsigmatrial(1, 0) = 0.;
    ddist_dsigmatrial(1, 1) = -(M_SQRT2 * sbeta * xi - M_SQRT2 * sqrt3 * sbeta * fPt) * fM / (3. * G);
    ddist_dsigmatrial(1, 2) = (M_SQRT2 * cbeta * xi - M_SQRT2 * sqrt3 * cbeta * fPt) * fM / (3. * G);
}

void TPZYCDruckerPragerPV::GradF2SigmaTrial(const TPZVec<REAL> &sigma_trial_pv, const REAL theta, const REAL beta, const REAL aProj, TPZFNMatrix<9, STATE> &ddist_dsigmatrial) const {
    fCap.GradSigmaTrial(sigma_trial_pv, theta, beta, aProj, ddist_dsigmatrial);
}

void TPZYCDruckerPragerPV::DF1Cart(STATE xi, STATE beta, TPZFNMatrix<6, STATE> &DFunccart) const {
    const REAL sqrt3 = sqrt(3);
    const REAL sqrt27 = pow(3, 1.5);
    const STATE sbeta = sin(beta);
    const STATE cbeta = cos(beta);

    DFunccart(0, 0) = -(2. * sqrt3 * cbeta * fM - sqrt27) / 9.;
    DFunccart(0, 1) = (2 * sqrt3 * sbeta * xi - 6 * sbeta * fPt) * fM / 9.;

    DFunccart(1, 0) = -((sqrt3 * sbeta - cbeta) * fM - 3.) / sqrt27;
    DFunccart(1, 1) = -((-3 * sbeta - sqrt27 * cbeta) * fPt + (sqrt3 * sbeta + 3 * cbeta) * xi) * fM / 9.;

    DFunccart(2, 0) = (3. + (sqrt3 * sbeta + cbeta) * fM) / sqrt27;
    DFunccart(2, 1) = -((sqrt27 * cbeta - 3. * sbeta) * fPt + (sqrt3 * sbeta - 3 * cbeta) * xi) * fM / 9.;
}

void TPZYCDruckerPragerPV::DF2Cart(STATE theta, STATE beta, STATE a, TPZFNMatrix<9, STATE> &DFunccart) const {
    fCap.DFuncCart(theta, beta, a, DFunccart);
}

void TPZYCDruckerPragerPV::ProjectSigmaDep(const TPZVec<REAL> &sigma_trial_pv, const REAL aPrev, TPZVec<REAL> &sigma, REAL &aProj, TPZFMatrix<REAL> &GradSigma) const {
    TPZVec<REAL> yield(NYield);
    this->Phi(sigma_trial_pv, aPrev, yield);
    STATE I1 = sigma_trial_pv[0] + sigma_trial_pv[1] + sigma_trial_pv[2];

    if (I1 / 3. >= fPt - aPrev) {
        if (yield[0] > 0.) {
        ProjectToSurfaceF1(sigma_trial_pv, aPrev, sigma, aProj, 1.e-5);
        // we can compute the tangent matrix
        TPZFNMatrix<4, STATE> jac(2, 2);
        TPZFNMatrix<6, STATE> ddist_dsigmatrial(2, 3), DFunccart(3, 2);
        STATE xi, beta;
        SurfaceParamF1(sigma, aProj, xi, beta);
        GradF1SigmaTrial(sigma_trial_pv, xi, beta, aProj, ddist_dsigmatrial);
        D2DistanceToSurfaceF1(sigma_trial_pv, xi, beta, jac);
        jac.Solve_LU(&ddist_dsigmatrial);
        DF1Cart(xi, beta, DFunccart);
        DFunccart.Multiply(ddist_dsigmatrial, GradSigma);
        GradSigma *= -1.;
        } else {
            sigma = sigma_trial_pv;
            aProj = aPrev;
            GradSigma.Identity();
        }
    } else {
        if (yield[1] > 0.) {
        ProjectToSurfaceF2(sigma_trial_pv, aPrev, sigma, aProj, 1.e-5);
        // we can compute the tangent matrix
        TPZFNMatrix<9, STATE> ddist_dsigmatrial(3, 3), jac(3, 3), DFunccart(3, 3);
        STATE theta, beta;
        SurfaceParamF2(sigma, aProj, theta, beta);
        GradF2SigmaTrial(sigma_trial_pv, theta, beta, aProj, ddist_dsigmatrial);
        D2DistanceToSurfaceF2(sigma_trial_pv, theta, beta, aProj, jac);
        jac.Solve_LU(&ddist_dsigmatrial);
        DF2Cart(theta, beta, aProj, DFunccart);
        DFunccart.Multiply(ddist_dsigmatrial, GradSigma);
        GradSigma *= -1.;
        } else {
            sigma = sigma_trial_pv;
            aProj = aPrev;
            GradSigma.Identity();
        }
    }
}

STATE TPZYCDruckerPragerPV::PlasticVolumetricStrain(STATE a) const {

    return fCap.PlasticVolumetricStrain(a);
}

TPZYCDruckerPragerPV::~TPZYCDruckerPragerPV() {
}

