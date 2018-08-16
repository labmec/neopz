/**
 * @file
 */


#include "TPZPlasticStepPV.h"
#include "TPZElasticResponse.h"

#include "TPZSandlerExtended.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZElasticResponse.h"

#include "pzlog.h"
#include "TPZYCCamClayPV.h"
#include "TPZYCDruckerPragerPV.h"

//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.material.TPZPlasticStepPV"));
//#endif
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.poroelastoplastic2"));
#endif

#ifdef LOG4CXX
static LoggerPtr logger2(Logger::getLogger("plasticity.poroelastoplastic"));
#endif

#define NewTangetQ

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::InitialDamage(const TPZTensor<REAL> & sigma, REAL & k){
    TPZTensor<REAL>::TPZDecomposed sigma_p;
    sigma.EigenSystem(sigma_p);
    k = fYC.InitialDamage(sigma_p.fEigenvalues);
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> * tangent) {
    bool require_tangent_Q = (tangent != NULL);
    
#ifdef PZDEBUG
    if (require_tangent_Q) {
        // Check for required dimensions of tangent
        if (tangent->Rows() != 6 || tangent->Cols() != 6) {
            std::cerr << "Unable to compute the tangent operator. Required tangent array dimensions are 6x6." << std::endl;
            DebugStop();
        }
    }
#endif

    // Initialization and spectral decomposition for the elastic trial stress state
    TPZTensor<REAL> eps_tr = epsTotal - fN.fEpsP;
    TPZTensor<REAL> sig_tr;
    fER.Compute(eps_tr, sig_tr);
    TPZTensor<REAL>::TPZDecomposed sig_eigen_system;
    sig_tr.EigenSystem(sig_eigen_system);
    
    int m_type = 0;
    STATE nextalpha = 0;
    TPZManVector<REAL, 3> sig_projected(3, 0.);
    
    // ReturMap in the principal values
    if (require_tangent_Q) {
        TPZTensor<REAL>::TPZDecomposed eps_eigen_system;
        eps_tr.EigenSystem(eps_eigen_system);
        
        TPZFMatrix<REAL> gradient(3, 3, 0.);
        fYC.ProjectSigma(sig_eigen_system.fEigenvalues, fN.fAlpha, sig_projected, nextalpha, m_type, &gradient);
        TangentOperator(gradient, eps_eigen_system, sig_eigen_system, *tangent);
    } else{
        fYC.ProjectSigma(sig_eigen_system.fEigenvalues, fN.fAlpha, sig_projected, nextalpha, m_type);
    }
    
    fN.fAlpha = nextalpha;
    fN.fMType = m_type;
    
#ifdef LOG4CXX_KEEP
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Sig Trial " << sigtrvec << "\nSig Project " << sigprvec << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    // Reconstruction of sigmaprTensor
    sig_eigen_system.fEigenvalues = sig_projected; // Under the assumption of isotropic material eigen vectors remain unaltered
    sigma = TPZTensor<REAL>(sig_eigen_system);
    
    TPZTensor<REAL> eps_e_Np1;
    fER.ComputeDeformation(sigma, eps_e_Np1);
    fN.fEpsT = epsTotal;
    fN.fEpsP = epsTotal - eps_e_Np1;
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStressComputeStrain(const TPZTensor<REAL> &sigma, TPZTensor<REAL> &epsTotal, TPZFMatrix<REAL> * tangent)
{
    
    bool require_tangent_Q = (tangent != NULL);
    
#ifdef PZDEBUG
    if (require_tangent_Q) {
        // Check for required dimensions of tangent
        if (tangent->Rows() != 6 || tangent->Cols() != 6) {
            std::cerr << "Unable to compute the tangent operator. Required tangent array dimensions are 6x6." << std::endl;
            DebugStop();
        }
    }
#endif
    
    TPZTensor<REAL>::TPZDecomposed sig_eigen_system;
    TPZTensor<REAL> sig_tr;
    
    // Initialization and spectral decomposition for the elastic trial stress state
    TPZTensor<REAL> eps_tr, eps_p_N, eps_e_Np1;
    eps_p_N = fN.fEpsP;
    eps_tr = epsTotal;
    eps_tr -= eps_p_N;
    fER.Compute(eps_tr, sig_tr);
    sig_tr.EigenSystem(sig_eigen_system);
    
    int m_type = 0;
    STATE nextalpha = 0;
    TPZManVector<REAL, 3> sig_projected(3, 0.);
    
    // ReturMap in the principal values
    if (require_tangent_Q) {
        // Required data when tangent is needed
        TPZTensor<REAL>::TPZDecomposed eps_eigen_system;
        TPZFMatrix<REAL> gradient(3, 3, 0.);
        
        eps_tr.EigenSystem(eps_eigen_system);
        
        fYC.ProjectSigma(sig_eigen_system.fEigenvalues, fN.fAlpha, sig_projected, nextalpha, m_type, &gradient);
        TangentOperator(gradient, eps_eigen_system, sig_eigen_system, *tangent);
    }
    else{
        fYC.ProjectSigma(sig_eigen_system.fEigenvalues, fN.fAlpha, sig_projected, nextalpha, m_type);
    }
    
    
    fN.fAlpha = nextalpha;
    fN.fMType = m_type;
    
#ifdef LOG4CXX_KEEP
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Sig Trial " << sigtrvec << "\nSig Project " << sigprvec << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    // Reconstruction of sigmaprTensor
    sig_eigen_system.fEigenvalues = sig_projected; // Under the assumption of isotropic material eigen vectors remain unaltered
//    sigma = TPZTensor<REAL>(sig_eigen_system);
    
    fER.ComputeDeformation(sigma, eps_e_Np1);
    fN.fEpsT = epsTotal;
    eps_p_N = epsTotal;
    eps_p_N -= eps_e_Np1; // plastic strain update
    fN.fEpsP = eps_p_N;
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) {
    
    TPZTensor<REAL>::TPZDecomposed sig_eigen_system, eps_eigen_system;
    TPZTensor<REAL> sigtr;

    TPZTensor<REAL> epsTr, epsPN, epsElaNp1;
    epsPN = fN.fEpsP;
    epsTr = epsTotal;
    epsTr -= epsPN; // Porque soh tem implementado o operator -=

    // Compute and Decomposition of SigTrial
    fER.Compute(epsTr, sigtr); // sigma = lambda Tr(E)I + 2 mu E
    epsTr.EigenSystem(eps_eigen_system);
    epsTr.ComputeEigenvectors(eps_eigen_system);
    sigtr.EigenSystem(sig_eigen_system);
    sigtr.ComputeEigenvectors(sig_eigen_system);

    TPZManVector<REAL, 3> sigtrvec(sig_eigen_system.fEigenvalues), sigprvec(3, 0.);

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
		sig_eigen_system.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
    STATE printPlastic = fN.Alpha();
#endif

    // ReturnMap in the principal values
    STATE nextalpha = 0;
    TPZFNMatrix<9> GradSigma(3, 3, 0.);
    fYC.ProjectSigmaDep(sigtrvec, fN.fAlpha, sigprvec, nextalpha, GradSigma);
    fN.fAlpha = nextalpha;

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Sig Trial " << sigtrvec << "\nSig Project " << sigprvec << std::endl;
        GradSigma.Print("GradSigma", sout, EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    // Reconstruction of sigmaprTensor
    sig_eigen_system.fEigenvalues = sigprvec; // updating the projected values used inside TangentOperator method.
    sigma = TPZTensor<REAL>(sig_eigen_system);
    TangentOperator(GradSigma, eps_eigen_system, sig_eigen_system, Dep);

    fER.ComputeDeformation(sigma, epsElaNp1);
    fN.fEpsT = epsTotal;
    epsPN = epsTotal;
    epsPN -= epsElaNp1; // Transforma epsPN em epsPNp1
    fN.fEpsP = epsPN;


#ifdef LOG4CXX
    if (logger2->isDebugEnabled()) {
        if (fabs(printPlastic - fN.fAlpha) > 1.e-4) {
            std::stringstream sout;
            TPZVec<STATE> phi;
            TPZTensor<STATE> epsElastic(fN.fEpsT);
            epsElastic -= fN.fEpsP;
            Phi(epsElastic, phi);
            sout << " \n phi = [";
            for (int i = 0; i < phi.size(); i++) {
                sout << phi[i] << " ";
            }

            sout << " ] " << endl;

            sout << " \n eigenvalues Sigma = [";
            for (int i = 0; i < 3; i++) {
                sout << sig_eigen_system.fEigenvalues[i] << " ";
            }

            sout << " ] " << endl;



            LOGPZ_DEBUG(logger2, sout.str())
        }
    }
#endif
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::TangentOperator(TPZFMatrix<REAL> & gradient,TPZTensor<REAL>::TPZDecomposed & eps_eigen_system, TPZTensor<REAL>::TPZDecomposed & sig_eigen_system, TPZFMatrix<REAL> & Tangent){
    

    //Montando a matriz tangente
    int kival[] = {0, 0, 0, 1, 1, 2};
    int kjval[] = {0, 1, 2, 1, 2, 2};
    REAL G = fER.G();
    REAL lambda = fER.Lambda();
    
    // Coluna da matriz tangente
    for (unsigned int k = 0; k < 6; ++k) {
        const unsigned int ki = kival[k];
        const unsigned int kj = kjval[k];
        for (unsigned int i = 0; i < 3; ++i) {
            for (unsigned int j = 0; j < 3; ++j) {
                REAL temp = 2 * G * eps_eigen_system.fEigenvectors[j][kj] * eps_eigen_system.fEigenvectors[j][ki];
                if (ki == kj) {
                    temp += lambda;
                } else {
                    temp *= 2.;
                }
                for (int l = 0; l < 6; ++l) {
                    const unsigned int li = kival[l];
                    const unsigned int lj = kjval[l];
                    Tangent(l, k) += temp * gradient(i, j) * eps_eigen_system.fEigenvectors[i][li] * eps_eigen_system.fEigenvectors[i][lj];
                }/// l
            }///j
        }///i
    }///k
    
    REAL deigensig = 0., deigeneps = 0.;
    TPZFNMatrix<36> RotCorrection(6, 6, 0.);
    // Correcao do giro rigido
    for (unsigned int i = 0; i < 2; ++i) {
        for (unsigned int j = i + 1; j < 3; ++j) {
            deigeneps = eps_eigen_system.fEigenvalues[i] - eps_eigen_system.fEigenvalues[j];
            deigensig = sig_eigen_system.fEigenvalues[i] - sig_eigen_system.fEigenvalues[j];
    
            TPZFNMatrix<9, REAL> tempMat(3, 3, 0.);
            REAL factor = 0.;
            if (!IsZero(deigeneps)) {
                factor = deigensig / deigeneps;
            } else {
                factor = fER.G() * (gradient(i, i) - gradient(i, j) - gradient(j, i) + gradient(j, j)); // expression C.20
            }
            tempMat = ProdT(eps_eigen_system.fEigenvectors[i], eps_eigen_system.fEigenvectors[j]) + ProdT(eps_eigen_system.fEigenvectors[j], eps_eigen_system.fEigenvectors[i]);
            
            // expression C.14
            for (unsigned int k = 0; k < 6; ++k) {
                const unsigned int ki = kival[k];
                const unsigned int kj = kjval[k];
                TPZFNMatrix<9> ColCorr(3, 3, 0.);
                TPZFNMatrix<6> ColCorrV(6, 1, 0.);
                if (ki == kj) {
                    ColCorr = (eps_eigen_system.fEigenvectors[j][ki] * eps_eigen_system.fEigenvectors[i][kj]) * factor * tempMat;
                } else {
                    ColCorr = (eps_eigen_system.fEigenvectors[j][ki] * eps_eigen_system.fEigenvectors[i][kj] + eps_eigen_system.fEigenvectors[j][kj] * eps_eigen_system.fEigenvectors[i][ki]) * factor * tempMat;
                }
                ColCorrV = FromMatToVoight(ColCorr);
                for (int l = 0; l < 6; l++) {
                    RotCorrection(l, k) += ColCorrV(l, 0);
                }
            }
        } // j
    } // i
    
    Tangent += RotCorrection;
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::TaylorCheck(TPZTensor<REAL> &EpsIni, TPZTensor<REAL> &deps, REAL kprev, TPZVec<REAL> &conv) {
    TPZTensor<REAL> eps1, eps2, SigmaTemp, Sigma1, Sigma2;
    TPZFNMatrix <36> dSigDe(6, 6, 0.);
    TPZStack<REAL> coef;

    fN.fEpsP.Scale(0.);
    fN.fEpsT.Scale(0.);
    fN.fAlpha = kprev;
    this->ApplyStrainComputeDep(EpsIni, SigmaTemp, dSigDe);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "EpsIni " << EpsIni << "\nSigmaTemp " << SigmaTemp << "\ndSidDe " << dSigDe << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fN.fEpsP.Scale(0.);
    fN.fEpsT.Scale(0.);
    fN.fAlpha = kprev;

    REAL scale = 1.;
    REAL alphatable[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (int i = 0; i < 6; i++) {
        alphatable[i] *= scale;
    }
    for (int ia = 0; ia < 5; ia++) {
        REAL alpha1 = alphatable[0];
        REAL alpha2 = alphatable[ia + 1];
        eps1.Scale(0.);
        eps2.Scale(0.);
        eps1 = EpsIni;
        eps2 = EpsIni;
        eps1.Add(deps, alpha1);
        eps2.Add(deps, alpha2);

        fN.fEpsT = EpsIni;
        this->ApplyStrainComputeSigma(eps1, Sigma1);
        fN.fEpsP.Scale(0.);
        fN.fEpsT.Scale(0.);
        fN.fAlpha = kprev;

        fN.fEpsT = EpsIni;
        this->ApplyStrainComputeSigma(eps2, Sigma2);
        fN.fEpsP.Scale(0.);
        fN.fEpsT.Scale(0.);
        fN.fAlpha = kprev;

        TPZFNMatrix <6> deps1(6, 1, 0.), deps2(6, 1, 0.);
        TPZFNMatrix <9> depsMat(3, 3, 0.);
        depsMat = deps;
        deps1 = FromMatToVoight(depsMat);
        deps2 = FromMatToVoight(depsMat);

        TPZFNMatrix <6> tanmult1(6, 1, 0.), tanmult2(6, 1, 0.);
        dSigDe.Multiply(deps1, tanmult1);
        dSigDe.Multiply(deps2, tanmult2);

        for (int i = 0; i < 6; i++) {
            tanmult1(i, 0) *= alpha1;
            tanmult2(i, 0) *= alpha2;
        }

        TPZFNMatrix <9> SigMatTemp33(3, 3, 0.);
        TPZFNMatrix <6> sigprMat(6, 1, 0.), sigpr1Mat(6, 1, 0.), sigpr2Mat(6, 1, 0.);
        SigMatTemp33 = SigmaTemp;
        sigprMat = FromMatToVoight(SigMatTemp33);
        SigMatTemp33 = Sigma1;
        sigpr1Mat = FromMatToVoight(SigMatTemp33);
        SigMatTemp33 = Sigma2;
        sigpr2Mat = FromMatToVoight(SigMatTemp33);

        TPZFNMatrix<6> error1(6, 1, 0.), error2(6, 1, 0.);
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sigprMat.Print("sigprMat", sout);
            sigpr1Mat.Print("sigpr1Mat", sout);
            tanmult1.Print("tanmult1", sout);
            sigpr2Mat.Print("sigpr2Mat", sout);
            tanmult2.Print("tanmult2", sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        for (int i = 0; i < 6; i++) {
            error1(i, 0) = sigpr1Mat(i, 0) - sigprMat(i, 0) - tanmult1(i, 0);
            error2(i, 0) = sigpr2Mat(i, 0) - sigprMat(i, 0) - tanmult2(i, 0);
        }
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            error1.Print("error1:", sout);
            error2.Print("error2:", sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        REAL n;
        REAL norm1, norm2;
        norm1 = NormVecOfMat(error1);
        norm2 = NormVecOfMat(error2);
        n = (log(norm1) - log(norm2)) / (log(alpha1) - log(alpha2));
        coef.push_back(n);
    }
    conv = coef;
    std::cout << "coef = " << coef << std::endl;
}

template <class YC_t, class ER_t >
REAL TPZPlasticStepPV<YC_t, ER_t>::ComputeNFromTaylorCheck(REAL alpha1, REAL alpha2, TPZFMatrix<REAL> &error1Mat, TPZFMatrix<REAL> &error2Mat)
{
    REAL norm1, norm2, n;
    norm1 = NormVecOfMat(error1Mat);
    norm2 = NormVecOfMat(error2Mat);
    n = log(norm1 / norm2) / log(alpha1 / alpha2);
    return n;
}

REAL NormVecOfMat(TPZFNMatrix <9> mat)
{
    REAL norm = 0.;
    for (int i = 0; i < mat.Rows(); i++) {
        norm += mat(i, 0) * mat(i, 0);
    }
    norm = sqrt(norm);
    return norm;
}

REAL InnerVecOfMat(TPZFMatrix<REAL> &m1,TPZFMatrix<REAL> &m2)
{
    REAL dot = 0.;
    for (int i = 0; i < m1.Rows(); i++) {
        dot += m1(i, 0) * m2(i, 0);
    }
    return dot;
}

TPZFMatrix<REAL> ProdT(TPZManVector<REAL,3> &v1, TPZManVector<REAL,3> &v2) {
    TPZFMatrix<REAL> mat(3, 3, 0.);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mat(i, j) = v1[i] * v2[j];
        }
    }
    return mat;
}

TPZFNMatrix <6> FromMatToVoight(TPZFNMatrix <9> mat)
{
    TPZFNMatrix <6> voi(6, 1, 0.);
    int k = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = i; j < 3; j++) {
            voi(k++, 0) = mat(i, j);
        }
    }
    return voi;
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrain(const TPZTensor<REAL> &epsTotal)
{

    std::cout << " \n this method is not implemented in PlasticStepPV. ";
    DebugStop();

}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyLoad(const TPZTensor<REAL> & GivenStress, TPZTensor<REAL> &epsTotal) {
    //@TODO: Refactor this code
    TPZPlasticState<STATE> prevstate = GetState();
    epsTotal = prevstate.fEpsP;
    TPZTensor<STATE> GuessStress, Diff, Diff2;
    TPZFNMatrix<36, STATE> Dep(6, 6, 0.0);
    TPZFNMatrix<6, STATE> DiffFN(6, 1);

    ApplyStrainComputeSigma(epsTotal, GuessStress, &Dep);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        Dep.Print("Dep = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    Diff = GivenStress - GuessStress;

    STATE norm = Norm(Diff), normprev;
    STATE tol = 1.e-7;
    int counter = 0;

    while (norm>tol && counter<30) {
        CopyFromTensorToFMatrix(Diff, DiffFN);
        Dep.Solve_LU(&DiffFN);
        CopyFromFMatrixToTensor(DiffFN, Diff);
        TPZTensor<STATE> epsprev(epsTotal);
        normprev = norm;
        STATE scale = 1.;
        int counter2 = 0;
        do {
            epsTotal = epsprev;
            epsTotal.Add(Diff, scale);

            ApplyStrainComputeSigma(epsTotal, GuessStress, &Dep);
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                Dep.Print("Dep = ", sout, EMathematicaInput);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif

            fN = prevstate;
            Diff2 = GivenStress-GuessStress;
            CopyFromTensorToFMatrix(Diff2, DiffFN);
            Dep.Solve_LU(&DiffFN);
            norm = Norm(Diff2);
            //scale*=0.5;
            counter2++;
        } while (norm >= normprev && counter2 < 30);
        Diff = Diff2;
        counter++;
    }
    ApplyStrainComputeSigma(epsTotal, GuessStress, &Dep);
}

template <class YC_t, class ER_t >
TPZPlasticState<STATE>  TPZPlasticStepPV<YC_t, ER_t>::GetState() const
{
    return fN;
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::Phi(const TPZTensor<STATE> &eps, TPZVec<REAL> &phi) const
{
    TPZTensor<STATE> sigma;
    fER.Compute(eps, sigma);
    TPZTensor<STATE>::TPZDecomposed DecSig;
    sigma.EigenSystem(DecSig);
    TPZVec<STATE> sigvec(3);
    sigvec[0] = DecSig.fEigenvalues[0];
    sigvec[1] = DecSig.fEigenvalues[1];
    sigvec[2] = DecSig.fEigenvalues[2];
    fYC.Phi(sigvec, fN.Alpha(), phi);
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::SetState(const TPZPlasticState<REAL> &state)
{
    fN = state;
}

template<class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::Write(TPZStream& buf, int withclassid) const {
    fYC.Write(buf, withclassid);
    fER.Write(buf, withclassid);
    buf.Write(&fResTol);
    buf.Write(&fMaxNewton);
    fN.Write(buf, withclassid);
}

template<class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::Read(TPZStream& buf, void* context) {
    fYC.Read(buf, context);
    fER.Read(buf, context);
    buf.Read(&fResTol);
    buf.Read(&fMaxNewton);
    fN.Read(buf, context);
}

/** @brief Object which represents the yield criterium */
//YC_t fYC;

/** @brief Object representing the elastic response */
//ER_t fER;

/** @brief Residual tolerance accepted in the plastic loop processes */
//REAL fResTol;

/** @brief Maximum number of Newton interations allowed in the nonlinear solvers */
//int fMaxNewton;	// COLOCAR = 30 (sugestao do erick!)




/** @brief Plastic State Variables (EpsT, EpsP, Alpha) at the current time step */
//TPZPlasticState<STATE> fN;

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::CopyFromFMatrixToTensor(TPZFMatrix<STATE> FNM,TPZTensor<STATE> &copy)
{
    FNM.Resize(6, 1);
    copy.XX() = FNM(0, 0);
    copy.XY() = FNM(1, 0);
    copy.XZ() = FNM(2, 0);
    copy.YY() = FNM(3, 0);
    copy.YZ() = FNM(4, 0);
    copy.ZZ() = FNM(5, 0);
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::CopyFromTensorToFMatrix(TPZTensor<STATE> tensor,TPZFMatrix<STATE> &copy)
{

    copy(0, 0) = tensor.XX();
    copy(1, 0) = tensor.XY();
    copy(2, 0) = tensor.XZ();
    copy(3, 0) = tensor.YY();
    copy(4, 0) = tensor.YZ();
    copy(5, 0) = tensor.ZZ();
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::SetElasticResponse(TPZElasticResponse &ER)
{
    fER = ER;
    fYC.SetElasticResponse(ER);
}

template class TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>;
template class TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>;
template class TPZPlasticStepPV<TPZYCCamClayPV, TPZElasticResponse>;
template class TPZPlasticStepPV<TPZYCDruckerPragerPV, TPZElasticResponse>;

/*
 // Correcao do giro rigido
 for (int i = 0; i < 2; i++) {
 for (int j = i+1; j<3 ; j++) {
 deigeneps = DecompEps.fEigenvalues[i]  - DecompEps.fEigenvalues[j];
 deigensig = sigprvec[i] - sigprvec[j];
 TPZFNMatrix<9,REAL> tempMat(3,1,0.);
 depsMat.Multiply(epsegveFromProj[i], tempMat);
 REAL deij = InnerVecOfMat(tempMat,epsegveFromProj[j]);
 REAL factor = 0.;
 if (!IsZero(deigeneps)) {
 factor = deigensig * deij / deigeneps;
 }
 else {
 factor = fER.G() * ( GradSigma(i,i) - GradSigma(i,j) - GradSigma(j,i) + GradSigma(j,j) ) * deij;
 }
 std::cout << "factor = " << factor << std::endl;
 std::cout << "G = " << fER.G() << std::endl;
 GradSigma.Print("GradSigma");
 tempMat.Redim(3, 3);
 tempMat = ProdT(epsegveFromProj[i],epsegveFromProj[j]) + ProdT(epsegveFromProj[j],epsegveFromProj[i]);
 factorMat += tempMat * factor;
 
 }
 }
 */
