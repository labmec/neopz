/**
 * @file
 */

#include "TPZPlasticStep.h"
#include "TPZElasticResponse.h"
#include "TPZLadeNelsonElasticResponse.h"
#include "TPZThermoForceA.h"
#include "TPZYCVonMises.h"
#include "TPZYCTresca.h"
#include "tpzyctrescaregularized.h"
#include "TPZYCLadeKim.h"
#include "TPZLadeKimThermoForceA.h"
#include "TPZYCSandlerDimaggio.h"
#include "TPZYCSandlerDimaggioL.h"
#include "TPZYCSandlerDimaggioL2.h"
#include "TPZSandlerDimaggio.h"
#include "TPZSandlerDimaggioThermoForceA.h"
#include "tpzycvonmisescombtresca.h"
#include "pzfmatrix.h"
#include "pzstepsolver.h"
#include "pzvec_extras.h"
#include "TPZPlasticState.h"
#include "TPZYCDruckerPrager.h"
#include "TPZYCRankine.h"


TPZPlasticIntegrMem<REAL, 4> teste;

using namespace std;

#include "tfad.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger pointloadconfig("plasticity.loadconfig");
static TPZLogger plasticIntegrLogger("plasticity.plasticIntegr");
#endif

#ifdef PZ_LOG
static TPZLogger logger("pz.PLASTIC_STEP.main");
static TPZLogger loggerx("pz.PLASTIC_STEP.main");
static TPZLogger loggerPlasticResidual("PLASTIC_RESIDUAL");
static TPZLogger loggerDEP1("pz.PLASTIC_STEP.DEP1");
static TPZLogger loggerDEP2("pz.PLASTIC_STEP.DEP2");
#endif

template <class YC_t, class TF_t, class ER_t>
int TPZPlasticStep<YC_t, TF_t, ER_t>::IntegrationSteps() const {
    return fPlasticMem.NElements() - 2;
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::SetState_Internal(const TPZPlasticState<REAL> & state) {
    fN = state;

#ifdef PZ_LOG
    {
        std::stringstream sout1, sout2;
        sout1 << ">>> SetState_Internal ***";
        sout2 << "\nfN.m_eps_p << " << fN.m_eps_p << "\nfN.m_hardening << " << fN.m_hardening;
        LOGPZ_DEBUG(logger, sout1.str().c_str());
        LOGPZ_DEBUG(logger, sout2.str().c_str());
    }
#endif
}

template <class YC_t, class TF_t, class ER_t>
TPZPlasticState<REAL> TPZPlasticStep<YC_t, TF_t, ER_t>::GetState_Internal() const {
    return fN;
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::SetState(
        const TPZPlasticState<REAL> & state) {
    int multipl = SignCorrection();
    fN = state;
    fN.m_eps_p *= multipl;
    fN.m_eps_t *= multipl;

#ifdef PZ_LOG
    {
        std::stringstream sout1, sout2;
        sout1 << ">>> SetUp ***";
        sout2 << "\nfN.m_eps_p << " << fN.m_eps_p << "\nfN.m_hardening << " << fN.m_hardening;
        LOGPZ_DEBUG(logger, sout1.str().c_str());
        LOGPZ_DEBUG(logger, sout2.str().c_str());
    }
#endif
}

template <class YC_t, class TF_t, class ER_t>
TPZPlasticState<REAL> TPZPlasticStep<YC_t, TF_t, ER_t>::GetState() const {
    int multipl = SignCorrection();
    TPZPlasticState<REAL> N(fN);
    N.m_eps_p *= multipl;
    N.m_eps_t *= multipl;

    return N;
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::SetUp(
        const TPZTensor<REAL> & epsTotal) {
    fN.m_eps_t = epsTotal;

#ifdef PZ_LOG
    {
        std::stringstream sout1, sout2;
        sout1 << ">>> SetUp ***";
        LOGPZ_DEBUG(logger, sout1.str().c_str());
    }
#endif
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyStrain(const TPZTensor<REAL> &epsTotal) {
    TPZTensor<REAL> epsTotal_Internal(epsTotal);
    epsTotal_Internal *= SignCorrection();
    ApplyStrain_Internal(epsTotal_Internal);
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma,  TPZFMatrix<REAL> * tangent) {
    
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
    
    int multipl = SignCorrection();
    TPZTensor<REAL> epsTotal_Internal(epsTotal);
    epsTotal_Internal *= multipl;
    ApplyStrainComputeSigma_Internal(epsTotal_Internal, sigma);
    sigma *= multipl;
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) {
    int multipl = SignCorrection();
    TPZTensor<REAL> epsTotal_Internal(epsTotal);
    epsTotal_Internal *= multipl;
    ApplyStrainComputeDep_Internal(epsTotal_Internal, sigma, Dep);
    sigma *= multipl;
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal) {
    TPZTensor<REAL> sigma_Internal(sigma);
    sigma_Internal *= SignCorrection();
    ApplyLoad_Internal(sigma_Internal, epsTotal);
    epsTotal *= SignCorrection();
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const {
    TPZTensor<REAL> epsTotal_Internal(epsTotal);
    epsTotal_Internal *= SignCorrection();

    Phi_Internal(epsTotal_Internal, phi);

    return;
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::SetTensionSign(int s) {
    fInterfaceTensionSign = (s >= 0) ? 1 : -1;
}

template <class YC_t, class TF_t, class ER_t>
int TPZPlasticStep<YC_t, TF_t, ER_t>::SignCorrection() const {
    return fMaterialTensionSign * fInterfaceTensionSign;
}

template <class YC_t, class TF_t, class ER_t>
template<class T>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ComputePlasticVars(const TPZPlasticState<T> & state_T,
        TPZTensor<T> & sigma_T,
        T & A_T)const {
    // variable representing the stress and elastic deformation
    TPZTensor<T> epsE_T(state_T.EpsT());
    // subtract the value of the current plastic deformation
    epsE_T.Add(state_T.EpsP(), T(-1.));
    // compute the stress of the elastic response
    fER.ComputeStress(epsE_T, sigma_T);
    // compute the value of the thermo dynamical force for the given damage variable
    A_T = fTFA.Compute(state_T.VolHardening());
}

template <class YC_t, class TF_t, class ER_t>
bool TPZPlasticStep<YC_t, TF_t, ER_t>::IsStrainElastic(const TPZPlasticState<REAL> &state)const {
#ifdef PZ_LOG
    {
        std::stringstream sout1, sout2;
        sout1 << ">>> IsStrainElastic ***";
        sout2 << "\nstate.EpsT() << " << state.EpsT();
        LOGPZ_DEBUG(logger, sout1.str().c_str());
        LOGPZ_DEBUG(logger, sout2.str().c_str());
    }
#endif

    TPZTensor<REAL> sigma;
    REAL A;

    ComputePlasticVars<REAL>(state, sigma, A);

    // compute the value of the yield functions
    TPZManVector<REAL, 10> phi(YC_t::NYield);

    fYC.Compute(sigma, A, phi, 0);
    // verify if any yield function indicates plastification
    int i;
    for (i = 0; i < YC_t::NYield; i++) {
        if (phi[i] > 0.) break;
    }


    // if we are in the elastic range
    if (i == YC_t::NYield) {
#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "*** IsStrainElastic *** Strain yet in the elastic range - no damage variable needs update.\nExiting method ApplyStrain."
                    << "\n Phi = ";
            for (int j = 0; j < YC_t::NYield; j++)sout << phi[j] << "  ";
            LOGPZ_DEBUG(logger, sout.str().c_str());
        }
#endif
        return true;
    }
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "*** IsStrainElastic *** Strain exceeds the elastic range - damage variables need update"
                << "\n Phi = ";
        for (int j = 0; j < YC_t::NYield; j++)sout << phi[j] << "  ";
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif
    return false;

}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyStrain_Internal(const TPZTensor<REAL> &epsTotal) {

#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << ">>> ApplyStrain_Internal *** Imposed epsTotal << " << epsTotal;
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << ">>> ApplyStrain_Internal *** Imposed epsTotal << " << epsTotal;
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif

    //ProcessStrainNoSubIncrement(epsTotal);
    ProcessStrain(epsTotal);

    int n = fPlasticMem.NElements();

    // load the integrated values as the current state
    TPZPlasticStep<YC_t, TF_t, ER_t>::SetState_Internal(fPlasticMem[n - 1].m_elastoplastic_state);

#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "*** ProcessStrain *** Exiting Method.";
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif

}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ProcessStrainNoSubIncrement(const TPZTensor<REAL> &epsTotal, const EElastoPlastic ep) {

#ifdef PZ_LOG
    {
        std::stringstream sout1;
        sout1 << ">>> ProcessStrainNoSubIncrement ***";
        LOGPZ_DEBUG(logger, sout1.str().c_str());
    }
#endif

    REAL yieldMultipl = 0;

    fPlasticMem.Resize(0);
    PushPlasticMem(fN,
            0. /*k*/,
            0. /*lambda*/,
            TPZManVector<REAL, YC_t::NYield>(YC_t::NYield, 0)/*delGamma*/,
            TPZManVector<int, YC_t::NYield>(YC_t::NYield, 0)/*validEqs*/,
            0 /*forceYield*/);

    TPZPlasticState<REAL> stateAtYield(fN), Np1(fN); // Np1 state with fN guesses
    Np1.m_eps_t = epsTotal;

    bool elastic = true;

    if (ep == EForceElastic) {
#ifdef PZ_LOG
        {
            std::stringstream sout1;
            sout1 << ">>> ProcessStrainNoSubIncrement *** behaviour imposed to be Elastic";
            LOGPZ_DEBUG(logger, sout1.str().c_str());
        }
#endif
        elastic = true;
    }

    if (ep == EForcePlastic) {
#ifdef PZ_LOG
        {
            std::stringstream sout1;
            sout1 << ">>> ProcessStrainNoSubIncrement *** behaviour imposed to be Plastic";
            LOGPZ_DEBUG(logger, sout1.str().c_str());
        }
#endif
        elastic = false;
    }

    if (ep == EAuto) elastic = IsStrainElastic(Np1);

    if (elastic) {
        PushPlasticMem(Np1,
                1. /*k*/,
                0. /*lambda unused - elastic state*/,
                TPZManVector<REAL, YC_t::NYield>(YC_t::NYield, 0)/*delGamma*/,
                TPZManVector<int, YC_t::NYield>(YC_t::NYield, 0)/*validEqs*/,
                0 /*forceYield*/);
        return;
    }

    // Plastic Integration needed

    yieldMultipl = FindPointAtYield(Np1.EpsT(), stateAtYield);

    PushPlasticMem(stateAtYield,
            yieldMultipl /*k*/,
            0. /*lambda unused - elastic state*/,
            TPZManVector<REAL, YC_t::NYield>(YC_t::NYield, 0)/*delGamma*/,
            TPZManVector<int, YC_t::NYield>(YC_t::NYield, 0)/*validEqs*/,
            0 /*forceYield*/);

    REAL multipl = 0.99;
    TPZTensor<REAL> DeltaEpsP_guess = Np1.m_eps_t;
    DeltaEpsP_guess.Add(stateAtYield.m_eps_t, -1.);
    Np1.m_eps_p.Add(DeltaEpsP_guess, multipl);


    //PlasticIntegrate(stateAtYield, Np1, fIntegrTol);
    int succeeded;

    TPZManVector<REAL, YC_t::NYield> delGamma(YC_t::NYield, 0.);
    TPZManVector<int, YC_t::NYield> validEqs(YC_t::NYield, 0);

    REAL normEpsPErr = 0.;
    //    int counter = 0;
    REAL lambda = 0.;

    succeeded = PlasticLoop(stateAtYield, Np1, delGamma, normEpsPErr, lambda, validEqs);

    if (normEpsPErr < fIntegrTol && succeeded) {
        PushPlasticMem(Np1, 1., lambda, delGamma, validEqs, fYC.GetForceYield());
        return;
    }

    TPZTensor<REAL> deltaEpsTotal(Np1.EpsT());
    deltaEpsTotal.Add(stateAtYield.EpsT(), -1.);
    TPZFMatrix<REAL> residual_mat(6, 1);

    REAL resnorm = 0.;
    do {



        Np1.m_eps_t.Add(deltaEpsTotal, 1.);
        succeeded = PlasticLoop(stateAtYield, Np1, delGamma, normEpsPErr, lambda, validEqs);

        TPZTensor<REAL> res(Np1.m_eps_t), epstyield(stateAtYield.m_eps_t);
        int k;
        for (k = 0; k < 6; k++)residual_mat(k, 0) = res.fData[k] - epstyield.fData[k];
        for (k = 0; k < 6; k++)resnorm += pow(residual_mat(k, 0), 2.);
        resnorm = sqrt(resnorm);


    } while (resnorm > 1.e-3 && succeeded);

    //    succeeded = PlasticLoop(stateAtYield, Np1, delGamma, normEpsPErr, lambda, validEqs);
    PushPlasticMem(Np1, 1., lambda, delGamma, validEqs, fYC.GetForceYield());

    cout << "\n ProcessStrainNoSubIncrement  ";
    cout << "\n fIntegrTol = " << fIntegrTol;
    cout << "\n lambda = " << lambda;
    cout << "\n delGamma = " << delGamma;
    cout << "\n fIntegrTol = " << fIntegrTol;

}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ProcessStrain(const TPZTensor<REAL> &epsTotal, const EElastoPlastic ep) {

#ifdef PZ_LOG
    {
        std::stringstream sout1;
        sout1 << ">>> ProcessStrain ***";
        LOGPZ_DEBUG(logger, sout1.str().c_str());
    }
#endif

    REAL yieldMultipl = 0;

    fPlasticMem.Resize(0);
    PushPlasticMem(fN,
            0. /*k*/,
            0. /*lambda*/,
            TPZManVector<REAL, YC_t::NYield>(YC_t::NYield, 0)/*delGamma*/,
            TPZManVector<int, YC_t::NYield>(YC_t::NYield, 0)/*validEqs*/,
            0 /*forceYield*/);

    TPZPlasticState<REAL> stateAtYield(fN),
            Np1(fN); // Np1 state with fN guesses
    Np1.m_eps_t = epsTotal;

    bool elastic = true;

    if (ep == EForceElastic) {
#ifdef PZ_LOG
        {
            std::stringstream sout1;
            sout1 << ">>> ProcessStrain *** behaviour imposed to be Elastic";
            LOGPZ_DEBUG(logger, sout1.str().c_str());
        }
#endif
        elastic = true;
    }

    if (ep == EForcePlastic) {
#ifdef PZ_LOG
        {
            std::stringstream sout1;
            sout1 << ">>> ProcessStrain *** behaviour imposed to be Plastic";
            LOGPZ_DEBUG(logger, sout1.str().c_str());
        }
#endif
        elastic = false;
    }

    if (ep == EAuto) {
        bool result = IsStrainElastic(Np1);
        elastic = (result == 1);
    }

    if (elastic) {
        PushPlasticMem(Np1,
                1. /*k*/,
                0. /*lambda unused - elastic state*/,
                TPZManVector<REAL, YC_t::NYield>(YC_t::NYield, 0)/*delGamma*/,
                TPZManVector<int, YC_t::NYield>(YC_t::NYield, 0)/*validEqs*/,
                0 /*forceYield*/);
        return;
    }

    // Plastic Integration needed

    yieldMultipl = FindPointAtYield(Np1.EpsT(), stateAtYield);

    PushPlasticMem(stateAtYield,
            yieldMultipl /*k*/,
            0. /*lambda unused - elastic state*/,
            TPZManVector<REAL, YC_t::NYield>(YC_t::NYield, 0)/*delGamma*/,
            TPZManVector<int, YC_t::NYield>(YC_t::NYield, 0)/*validEqs*/,
            0 /*forceYield*/);

    REAL multipl = 0.99;
    TPZTensor<REAL> DeltaEpsP_guess = Np1.m_eps_t;
    DeltaEpsP_guess.Add(stateAtYield.m_eps_t, -1.);
    Np1.m_eps_p.Add(DeltaEpsP_guess, multipl);

    PlasticIntegrate(stateAtYield, Np1, fIntegrTol);



}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyStrainComputeDep_Internal(const TPZTensor<REAL> &epsTotal,
        TPZTensor<REAL> &sigma,
        TPZFMatrix<REAL> &Dep) {

#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << ">>> ApplyStrainComputeDep_Internal *** Imposed epsTotal << " << epsTotal;
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << ">>> ApplyStrainComputeDep_Internal *** Imposed epsTotal << " << epsTotal;
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif

    ApplyStrain_Internal(epsTotal);

#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "*** ApplyStrainComputeDep *** \n Calling ComputeDep";
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif

    ComputeDep(sigma, Dep);
    /*
    TPZTensor<STATE> sigma2;
    TPZFNMatrix<36,STATE> Dep2(6,6);
    ComputeDep2(sigma2, Dep2);
#ifdef PZ_LOG
    {
        std::stringstream sout1;
        sout1 << "Imposed epsTotal << " << epsTotal
        << "\nResulted in Sigma = " << sigma << "\nDep = \n" << Dep << std::endl;
        
        LOGPZ_DEBUG(loggerDEP1,sout1.str());
    }
#endif
#ifdef PZ_LOG
    {
        std::stringstream sout1;
        sout1 << "Imposed epsTotal << " << epsTotal
        << "\nResulted in Sigma = " << sigma2 << "\nDep = \n" << Dep2 << std::endl;
        
        LOGPZ_DEBUG(loggerDEP2,sout1.str());
    }
#endif
     */
}

/*
template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ComputeDep2(TPZTensor<REAL> & sigma, TPZFMatrix<REAL> &Dep)
{
    
    const int nyield = YC_t::NYield;
    const int nVarsResidual = 7+nyield;
    const int nVarsTensor = 6;
    
    int i,j;
    
    typedef TFad<nVarsTensor, REAL> TFAD;
    typedef TFad<nVarsResidual, TFAD> TFAD_FAD;
    
    REAL normEpsPErr, resnorm;
    
    EStatus status;
    
    TPZPlasticState<TFAD_FAD> Nk_FADFAD,
    Nkp1_FADFAD;
    
    TPZPlasticState<REAL > Nk, Nkp1;
    
    TPZManVector<TFAD_FAD, nyield>  delGamma_FADFAD(nyield);
    TPZManVector<REAL, nyield>              delGamma(nyield);
    TPZManVector<TFAD_FAD, nVarsResidual> epsRes_FADFAD(nVarsResidual);
    TPZManVector<REAL, nVarsResidual> epsRes(nVarsResidual);
    TPZDiffMatrix< TFAD > tangent_FAD(nVarsResidual,nVarsResidual),
    residual_FAD(nVarsResidual,1),
    Sol_FAD(nVarsResidual,1);
    
    TPZPlasticState<TFAD> Nk_FAD;
    TPZTensor<TFAD> sigma_FAD;
    TFAD A_FAD;
    
    TPZTensor<REAL> diffPlasticStrain;
    
    int n = fPlasticMem.NElements();
    
    if(n < 2)
    {
#ifdef PZ_LOG
        if(plasticIntegrLogger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << ">>> ComputeDep *** Insufficient Plastic Mem Entries: " << n << ".";
            LOGPZ_ERROR(plasticIntegrLogger,sout.str().c_str());
        }
#endif
#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << ">>> ComputeDep *** Insufficient Plastic Mem Entries: " << n << ".";
            LOGPZ_ERROR(logger,sout.str().c_str());
        }
#endif
        return;
    }
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << ">>> ComputeDep *** Plastic Mem Entries: " << n << ".";
        if( n == 2 ) sout << "\nTwo entries indicate pure elastic step";
        LOGPZ_DEBUG(logger,sout.str().c_str());
    }
#endif
    
    if(n==2)
    {// pure elastic step - no plastic loop necessary
        fPlasticMem[1].m_elastoplastic_state.CopyTo(Nk_FAD);
        for(i = 0; i < nVarsTensor; i++)Nk_FAD.m_eps_t.fData[i].diff(i,nVarsTensor);
        
    }else
    {// plastification occurs
        
        // plastic Loop analogue with derivative evaluation
        
        // Initializing last available elastic step
        fPlasticMem[1].m_elastoplastic_state.CopyTo(Nk_FADFAD);
        for(j = 0; j < nVarsTensor; j++)Nk_FADFAD.m_eps_t.fData[j].val().fastAccessDx(j) = fPlasticMem[1].fK;
        
        TPZManVector<STATE,10> pivots(nVarsResidual,0.);
        
        REAL disturbFactor = sqrt(fResTol); // imposing an initial guess very close to the real
        // solution in order to ensure residual drop and convergence.
        // This is very important to accumulate good FAD derivatives
        
        for(i = 2; i < n; i++)
        {
            InitializePlasticFAD(fPlasticMem[i].m_elastoplastic_state,
                                 fPlasticMem[i].fDelGamma,
                                 Nkp1_FADFAD,
                                 delGamma_FADFAD,
                                 nVarsResidual);
            fYC.SetForceYield(fPlasticMem[i].fForceYield); // imposing the same assumptions made in the plasticLoop
            
            Nkp1_FADFAD.m_hardening.val().val() -=
            ( fPlasticMem[i].m_elastoplastic_state.m_hardening - fPlasticMem[i-1].m_elastoplastic_state.m_hardening )
 * disturbFactor;
            
            for(j = 0; j < nyield; j++)delGamma_FADFAD[j].val().val() *= (1.0 - disturbFactor);
            
            for(j = 0; j < nVarsTensor; j++)
            {
                Nkp1_FADFAD.m_eps_p.fData[j].val().val() -=
                ( fPlasticMem[i].m_elastoplastic_state.m_eps_p.fData[j] - fPlasticMem[i-1].m_elastoplastic_state.m_eps_p.fData[j] )
 * disturbFactor;
                Nkp1_FADFAD.m_eps_t.fData[j].val().diff(j, nVarsTensor);
                Nkp1_FADFAD.m_eps_t.fData[j].val().fastAccessDx(j) = fPlasticMem[i].fK; // setting the derivatives of EpsT with respect to deltaEpsT
            }
            
#ifdef PZ_LOG
            {
                std::stringstream sout;
                sout << "*** ComputeDep *** Before Matrix Invertion: Plastic step number " << i-1 << " of " << n-1
                << "\n*Nk_FADFAD=\n" << setw(10) << Nk_FADFAD
                << "\n*Nkp1_FADFAD=\n" << setw(10) << Nkp1_FADFAD
                << "\n*delGamma_FADFAD=\n" << setw(10) << delGamma_FADFAD;
                LOGPZ_DEBUG(logger,sout.str().c_str());
            }
#endif
            
            int NewtonCounter = 0;
            
            do{//(resnorm > fResTol && NewtonCounter < fMaxNewton)
                
                PlasticResidual<TFAD_FAD, TFAD_FAD>(Nk_FADFAD, Nkp1_FADFAD,
                                                    delGamma_FADFAD, epsRes_FADFAD,
                                                    normEpsPErr);
                
                
//                PlasticResidualRK<TFAD_FAD, TFAD_FAD>(Nk_FADFAD, Nkp1_FADFAD,
//                 delGamma_FADFAD, epsRes_FADFAD,
//                 normEpsPErr);
                
                tangent_FAD.Reset(); // resets the LU Decomposition flag
                
                int precond = 1;
                int resetInvalidEqs = 1;
                ExtractTangent(epsRes_FADFAD, residual_FAD,
                               resnorm, tangent_FAD, // TPZFMatrix for T1=fad<real> type
                               fPlasticMem[i].fValidEqs, pivots,
                               precond, resetInvalidEqs);
#ifdef PZ_LOG
                TPZDiffMatrix< TFAD > tangentcopy(tangent_FAD);
#endif
                status = tangent_FAD.Decompose_LU();
                if(status == EZeroPivot)
                {
                    std::stringstream sout;
                    sout << "*** ComputeDep *** ### Decompose_LU error! - ZeroPivot ### No inversion will be performed";
#ifdef PZ_LOG
                    sout << "\nMatrix before decomposition\n";
                    sout << tangentcopy;
                    LOGPZ_ERROR(plasticIntegrLogger,sout.str().c_str());
#endif
#ifdef PZ_LOG
                    LOGPZ_ERROR(logger,sout.str().c_str());
#endif
                    cout << endl << sout.str().c_str();
                }
                
                Sol_FAD = residual_FAD;
                status = tangent_FAD.Substitution(&Sol_FAD);
                if(status != EOk)
                {
                    std::stringstream sout;
                    if(status == EIncompDim)sout << "*** ComputeDep *** ### LU Substitution error! - IncompatibleDimensions ### No inversion will be performed";
                    if(status == EZeroPivot)sout << "*** ComputeDep *** ### LU Substitution error! - ZeroPivot ### No inversion will be performed";
#ifdef PZ_LOG
                    LOGPZ_ERROR(plasticIntegrLogger,sout.str().c_str());
#endif
#ifdef PZ_LOG
                    LOGPZ_ERROR(logger,sout.str().c_str());
#endif
                    cout << endl << sout.str().c_str();
                }
                
                
//                 REAL lambda = UpdatePlasticVars(Nk_FADFAD, Nkp1_FADFAD,
//                 delGamma_FADFAD, epsRes_FADFAD,
//                 Sol_FAD, fPlasticMem[i].fValidEqs,
//                 0); //Do not update variables internally
                
                
                REAL lambda = 1.; // forcing unity because the guess is very close to the solution
                
                
#ifdef PZ_LOG
                {
                    std::stringstream sout;
                    sout << "*** ComputeDep *** After " << NewtonCounter
                    << "-th Matrix Invertion: Plastic step number " << i-1 << " of " << n-1
                    << "\nSol_FAD=\n " << setw(10) << Sol_FAD
                    << "\nepsRes_FADFAD=\n " << setw(10) << epsRes_FADFAD
                    << "\nNkp1_FADFAD=\n" << setw(10) <<  Nkp1_FADFAD
                    << "\ndelGamma_FADFAD=\n" << delGamma_FADFAD
                    << "\nfLambda = " << setw(10) << lambda;
                    LOGPZ_DEBUG(logger,sout.str().c_str());
                }
#endif
                
                for(j=0; j<nVarsTensor; j++) Nkp1_FADFAD.m_eps_p.fData[j].val() -= lambda*Sol_FAD(j);
                Nkp1_FADFAD.m_hardening.val() -= lambda*Sol_FAD(j++);
                for(j=0; j<YC_t::NYield; j++) delGamma_FADFAD[j].val() -= lambda*Sol_FAD(j+7);
                
                for(j=0;j<nVarsTensor;j++)
                {
                    Nk  .m_eps_p.fData[j] = Nk_FADFAD  .m_eps_p.fData[j].val().val();
                    Nk  .m_eps_t.fData[j] = Nk_FADFAD  .m_eps_t.fData[j].val().val();
                    Nkp1.m_eps_p.fData[j] = Nkp1_FADFAD.m_eps_p.fData[j].val().val();
                    Nkp1.m_eps_t.fData[j] = Nkp1_FADFAD.m_eps_t.fData[j].val().val();
                }
                Nk  .m_hardening = Nk_FADFAD.  m_hardening.val().val();
                Nkp1.m_hardening = Nkp1_FADFAD.m_hardening.val().val();
                for(j=0; j<YC_t::NYield; j++) delGamma[j] = delGamma_FADFAD[j].val().val();
                
                PlasticResidual<REAL, REAL>(Nk, Nkp1, delGamma, epsRes, normEpsPErr);
                //PlasticResidual<TFAD_FAD, TFAD_FAD>(Nk_FADFAD, Nkp1_FADFAD, delGamma_FADFAD, epsRes_FADFAD, normEpsPErr);
                
                // updating the residual
                for(j=0; j < nyield; j++)
                    if(fPlasticMem[i].fValidEqs[j] == 0)epsRes[j + 7] = 0.;
                resnorm = 0.;
                for(j=0; j < nVarsResidual; j++)
                    resnorm += pow(epsRes[j] , 2.);
                resnorm = sqrt(resnorm);
                
                NewtonCounter++;
            }while((resnorm > fResTol && NewtonCounter < fMaxNewton) || NewtonCounter < 2);
            
            for(j = 0; j < 6; j++)diffPlasticStrain.fData[j] = Nkp1_FADFAD.m_eps_p.fData[j].val().val()
                - fPlasticMem[i].m_elastoplastic_state.m_eps_p.fData[j];
            
#ifdef PZ_LOG
            {
                std::stringstream sout;
                sout << "*** ComputeDep *** substep " << i-1 << " of " << n-2
                << " solved with " << NewtonCounter << " iterations and residual = " << resnorm;
                if(resnorm > fResTol)
                {
                    sout << "\n#### Truncated Newton ####. Results are unpredictable";
                    sout << "\nDifferences in the plastic strain:" << diffPlasticStrain;
                    sout << "\nDifferences in alpha:" << Nkp1_FADFAD.m_hardening.val().val() - fPlasticMem[i].m_elastoplastic_state.m_hardening;
                    sout << "\nAlpha = " << fPlasticMem[i].m_elastoplastic_state.m_hardening;
                    LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
                }else
                {
                    LOGPZ_DEBUG(plasticIntegrLogger,sout.str().c_str());
                }
            }
#endif
            
#ifdef PZ_LOG
            {
                std::stringstream sout;
                sout << "*** ComputeDep *** substep " << i-1 << " of " << n-2
                << " solved with " << NewtonCounter << " iterations and residual = " << resnorm;
                if(resnorm > fResTol)
                {
                    sout << "\n#### Truncated Newton ####. Results are unpredictable";
                    sout << "\nDifferences in the plastic strain:" << diffPlasticStrain;
                    sout << "\nDifferences in alpha:" << Nkp1_FADFAD.m_hardening.val().val() - fPlasticMem[i].m_elastoplastic_state.m_hardening;
                    sout << "\nAlpha = " << fPlasticMem[i].m_elastoplastic_state.m_hardening;
                    LOGPZ_WARN(logger,sout.str().c_str());
                }else
                {
                    LOGPZ_DEBUG(logger,sout.str().c_str());
                }
            }
#endif
            
            // substep marching
            for(j = 0; j < nVarsTensor; j++)
            {
                Nk_FADFAD.m_eps_t.fData[j].val()  = Nkp1_FADFAD.m_eps_t.fData[j].val(); // ignoring first derivatives
                Nk_FADFAD.m_eps_p.fData[j].val()  = Nkp1_FADFAD.m_eps_p.fData[j].val(); // ignoring first derivatives
            }
            Nk_FADFAD.m_hardening.val() = Nkp1_FADFAD.m_hardening.val();// ignoring first derivatives
        }
        
        // decreasing derivative order
        for(i = 0; i < nVarsTensor; i++)
        {
            Nk_FAD.m_eps_t.fData[i]  = Nk_FADFAD.m_eps_t.fData[i].val();
            Nk_FAD.m_eps_p.fData[i]  = Nk_FADFAD.m_eps_p.fData[i].val();
        }
        Nk_FAD.m_hardening = Nk_FADFAD.m_hardening.val();
        
    }// end of plasticLoop analogue
    
    //at this point, all the residual variables should contain their real derivatives to detaEpsT
    // in the case the loading is either elastic or elastoplastic.
    ComputePlasticVars(Nk_FAD, sigma_FAD, A_FAD);
    
    for(i = 0; i < nVarsTensor; i++)for(j = 0; j < nVarsTensor; j++)Dep(i,j) = sigma_FAD.fData[i].dx(j);
    
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "*** ComputeDep *** \nsigma_FAD= \n" << sigma_FAD
        << "\nDep =\n" << Dep;
        LOGPZ_DEBUG(logger,sout.str().c_str());
    }
#endif
    
    sigma_FAD.CopyTo(sigma);
    
}
 */

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ComputeDep(TPZTensor<REAL> & sigma, TPZFMatrix<REAL> &Dep) {

    const int nyield = YC_t::NYield;
    const int nVarsResidual = 7 + nyield;
    const int nVarsTensor = 6;

    int j;

    typedef TFad<nVarsTensor, REAL> TFAD;
    typedef TFad<nVarsResidual, REAL> TFAD_RES;

    REAL normResidual, normResidual2, resnorm;

    EStatus status;

    TPZPlasticState<TFAD_RES> Nk_FADRES, Nkp1_FADRES;

    TPZPlasticState<REAL > Nk, Nkp1;

    TPZManVector<TFAD_RES, nyield> delGamma_FADRES(nyield);
    TPZManVector<REAL, nyield> delGamma(nyield);

    TPZManVector<TFAD, nyield> delGamma_FAD(nyield);

    TPZManVector<TFAD_RES, nVarsResidual> Residual_FADRES(nVarsResidual);
    TPZManVector<TFAD, nVarsResidual> Residual_FAD(nVarsResidual);
    TPZManVector<REAL, nVarsResidual> epsRes(nVarsResidual);
    TPZDiffMatrix< REAL > tangent(nVarsResidual, nVarsResidual), residual_RES(nVarsResidual, 1), Sol_RES(nVarsResidual, 1);
    TPZDiffMatrix<STATE> tangentFAD(nVarsResidual, 6), residualFAD_STATE(nVarsResidual, 1);

    TPZPlasticState<TFAD> Nk_FAD, Nkp1_FAD;
    TPZTensor<TFAD> sigma_FAD;
    TFAD A_FAD;

    TPZTensor<REAL> diffPlasticStrain;

    int n = fPlasticMem.NElements();

    if (n < 2) {
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << ">>> ComputeDep2 *** Insufficient Plastic Mem Entries: " << n << ".";
            LOGPZ_ERROR(plasticIntegrLogger, sout.str().c_str());
        }
#endif
#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << ">>> ComputeDep2 *** Insufficient Plastic Mem Entries: " << n << ".";
            LOGPZ_ERROR(logger, sout.str().c_str());
        }
#endif
        return;
    }
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << ">>> ComputeDep2 *** Plastic Mem Entries: " << n << ".";
        if (n == 2) sout << "\nTwo entries indicate pure elastic step";
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    if (n == 2) {// pure elastic step - no plastic loop necessary
        fPlasticMem[1].m_elastoplastic_state.CopyTo(Nkp1_FAD);
        for (int i = 0; i < nVarsTensor; i++)Nkp1_FAD.m_eps_t.fData[i].fastAccessDx(i) = 1.;

    } else {// plastification occurs

        // plastic Loop analogue with derivative evaluation

        // Initializing last available elastic step
        fPlasticMem[1].m_elastoplastic_state.CopyTo(Nk_FADRES);
        fPlasticMem[1].m_elastoplastic_state.CopyTo(Nkp1_FAD);
        //        for(j = 0; j < nVarsTensor; j++)Nk_FADFAD.m_eps_t.fData[j].val().fastAccessDx(j) = fPlasticMem[1].fK;

        // solution in order to ensure residual drop and convergence.
        // This is very important to accumulate good FAD derivatives
        TPZManVector<STATE, 10> pivots(nVarsResidual, 0.);


        for (int plasticstep = 2; plasticstep < n; plasticstep++) {
            InitializePlasticFAD(fPlasticMem[plasticstep].m_elastoplastic_state,
                    fPlasticMem[plasticstep].fDelGamma,
                    Nkp1_FADRES,
                    delGamma_FADRES,
                    nVarsResidual);
            fYC.SetForceYield(fPlasticMem[plasticstep].fForceYield); // imposing the same assumptions made in the plasticLoop
            Nk_FAD = Nkp1_FAD;
            fPlasticMem[plasticstep].m_elastoplastic_state.CopyTo(Nkp1_FAD);
            for (int j = 0; j < 6; j++) {
                Nkp1_FAD.m_eps_t[j].fastAccessDx(j) = fPlasticMem[plasticstep].fK;
            }
            for (int i = 0; i < YC_t::NYield; i++) {
                delGamma_FAD[i] = fPlasticMem[plasticstep].fDelGamma[i];
            }



#ifdef PZ_LOG
            {
                std::stringstream sout;
                sout << "*** ComputeDep2 *** Before Matrix Invertion: Plastic step number " << plasticstep - 1 << " of " << n - 1
                        << "\n*Nk_FADFAD = " << setw(10) << Nk_FADRES
                        << "\n*Nkp1_FADFAD = " << setw(10) << Nkp1_FADRES
                        << "\n*delGamma_FADFAD = " << setw(10) << delGamma_FADRES
                        << "\nNk_FAD= " << Nk_FAD
                        << "\nNkp1_FAD = " << Nkp1_FAD
                        << "\ndelGamma_FAD = " << delGamma_FAD;

                LOGPZ_DEBUG(logger, sout.str())
            }
#endif

            PlasticResidual<TFAD_RES, TFAD_RES>(Nk_FADRES, Nkp1_FADRES,
                    delGamma_FADRES, Residual_FADRES,
                    normResidual);
            PlasticResidual<TFAD, TFAD>(Nk_FAD, Nkp1_FAD, delGamma_FAD, Residual_FAD, normResidual2);


            tangent.Reset(); // resets the LU Decomposition flag

            int precond = 0;
            int precondtangent = 1;
            int resetInvalidEqs = 1;
            ExtractTangent(Residual_FADRES, residual_RES,
                    resnorm, tangent, // TPZFMatrix for T1=fad<real> type
                    fPlasticMem[plasticstep].fValidEqs, pivots,
                    precondtangent, resetInvalidEqs);

            ExtractTangent(Residual_FAD, residualFAD_STATE, resnorm, tangentFAD, fPlasticMem[plasticstep].fValidEqs, pivots, precond, resetInvalidEqs);

            for (int i = 0; i < tangentFAD.Rows(); i++) {
                for (int j = 0; j < tangentFAD.Cols(); j++) {
                    tangentFAD(i, j) /= pivots[i];
                }
            }
#ifdef PZ_LOG
            TPZDiffMatrix< REAL > tangentcopy(tangent);
#endif
            status = tangent.Decompose_LU();
            if (status == EZeroPivot) {
                std::stringstream sout;
                sout << "*** ComputeDep2 *** ### Decompose_LU error! - ZeroPivot ### No inversion will be performed";
#ifdef PZ_LOG
                sout << "\nMatrix before decomposition\n";
                sout << tangentcopy;
                LOGPZ_ERROR(plasticIntegrLogger, sout.str().c_str());
#endif
#ifdef PZ_LOG
                LOGPZ_ERROR(logger, sout.str().c_str());
#endif
                cout << endl << sout.str().c_str();
            }

            Sol_RES = tangentFAD;
            status = tangent.Substitution(&Sol_RES);
            if (status != EOk) {
                std::stringstream sout;
                if (status == EIncompDim)sout << "*** ComputeDep2 *** ### LU Substitution error! - IncompatibleDimensions ### No inversion will be performed";
                if (status == EZeroPivot)sout << "*** ComputeDep2 *** ### LU Substitution error! - ZeroPivot ### No inversion will be performed";
#ifdef PZ_LOG
                LOGPZ_ERROR(plasticIntegrLogger, sout.str().c_str());
#endif
#ifdef PZ_LOG
                LOGPZ_ERROR(logger, sout.str());
#endif
                cout << endl << sout.str().c_str();
            }



#ifdef PZ_LOG
            {
                std::stringstream sout;
                sout << "*** ComputeDep2 *** Plastic step number " << plasticstep - 1 << " of " << n - 1
                        << "\nSol_FAD = " << setw(10) << Sol_RES
                        << "\nResidual_FADRES = " << setw(10) << Residual_FADRES
                        << "\ntangentFAD = " << tangentFAD;

                LOGPZ_DEBUG(logger, sout.str().c_str());
            }
#endif

            for (j = 0; j < nVarsTensor; j++) for (int k = 0; k < 6; k++) Nkp1_FAD.m_eps_p.fData[j].fastAccessDx(k) = -Sol_RES(j, k);
            for (int k = 0; k < 6; k++) Nkp1_FAD.m_hardening.fastAccessDx(k) = -Sol_RES(j, k);
            j++;
            for (int k = 0; k < 6; k++) for (j = 0; j < YC_t::NYield; j++) delGamma_FAD[j].fastAccessDx(k) = -Sol_RES(j + 7, k);


#ifdef PZ_LOG
            {
                std::stringstream sout;
                sout << "*** ComputeDep2 *** substep " << plasticstep - 1 << " of " << n - 2
                        << " solved with residual = " << resnorm;
            }
#endif

        }
    } // end of condition on number of plastic steps

    //at this point, all the residual variables should contain their real derivatives to detaEpsT
    // in the case the loading is either elastic or elastoplastic.
    ComputePlasticVars(Nkp1_FAD, sigma_FAD, A_FAD);

    for (int i = 0; i < nVarsTensor; i++) for (j = 0; j < nVarsTensor; j++)Dep(i, j) = sigma_FAD.fData[i].dx(j);

#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "*** ComputeDep2 *** \nsigma_FAD= \n" << sigma_FAD
                << "\nDep =\n" << Dep;
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif

    sigma_FAD.CopyTo(sigma);

}

//template <class YC_t, class TF_t, class ER_t>
//REAL TPZPlasticStep<YC_t, TF_t, ER_t>::FindPointAtYield(
//                                                        const TPZTensor<REAL> &epsTotalNp1,
//                                                        TPZPlasticState<REAL> &stateAtYield)const
//{
//    
//#ifdef PZ_LOG
//    
//    int plasticIntegrOutput;
//    
//    {
//        std::stringstream sout;
//        sout << ">>> FindPointAtYield ***";
//        LOGPZ_DEBUG(logger,sout.str().c_str());
//    }
//#endif
//    
//    const int nVars = 1;
//    typedef TFad<nVars, REAL> TFAD_One;
//    TPZTensor<REAL> deltaEps;
//    TPZTensor<TFAD_One>     deltaEps_FAD, sigma_FAD;
//    TFAD_One A_FAD, multipl_FAD;
//    TPZManVector<TFAD_One,10> phi_FAD(YC_t::NYield);
//    REAL multiplN, minMultipl = 1.;
//    TPZPlasticState<TFAD_One> stateAtYield_FAD;
//    
//    stateAtYield.CopyTo(stateAtYield_FAD);
//    
//    
//    // delta epsilon contains epstotal_N+1 - epsTotal_N
//    epsTotalNp1.CopyTo(deltaEps);
//    deltaEps.Add(fN.EpsT(), -1.);
//    deltaEps.CopyTo(deltaEps_FAD);
//    
//    int i, nyield = YC_t::NYield;
//    for(i = 0; i < nyield; i++)
//    { // searching for the multiplier for each yield surface
//        
//        multiplN = 0.; // the plastic multiplier is likely to be zero
//        // because of possible previous plastifications,
//        // although a null value may lead to no derivative
//        // or numerical instabilities when deltaEpsT=0
//        int count = 0;
//        do
//        {
//            if(count  > 0)
//            {
//                REAL derX = phi_FAD[i].dx(0);
//                if(fabs(derX)<1.e-8) derX += 1.e-8;
//                multiplN -= phi_FAD[i].val() / derX; // avoiding division by zero
//                if(multiplN > 1.0)
//                {
//#ifdef PZ_LOG
//                    {
//                        std::stringstream sout;
//                        sout << "*** FindPointAtYield *** multiplication factor = " << multiplN << " set to 1.0";
//                        sout << "\nPlastification is known to occur within this load step and the guess for the nest Newton's step is clipped to multipl=1";
//                        LOGPZ_DEBUG(logger,sout.str().c_str());
//                    }
//#endif
//                    multiplN = 1.0; // The step is known as plastic in advance
//                    // (IsStrainElastic previously called)
//                    // The check above avoids the multiplicator to be too large when the
//                    // first derivative evaluation is very low. In very nonlinear models
//                    // and specially in those where the stiffness matrix is very low at the
//                    // null stress state this could happen and, if this statement weren't here,
//                    // this newton loop could be extremely slow to converge since the
//                    // initial guesses would be too far from the correct answer.
//                }
//                
//                if(multiplN < -1.0)
//                {
//#ifdef PZ_LOG
//                    {
//                        std::stringstream sout;
//                        sout << "*** FindPointAtYield *** multiplication factor = " << multiplN << " set to 1.0";
//                        sout << "\nPlastification is known to occur within this load step. Forcing restart with initial guess biased towards the highest range attempting to help the code reach physically meaningful results";
//                        LOGPZ_DEBUG(logger,sout.str().c_str());
//                    }
//#endif
//                    multiplN = 1.0;  // This check attempts to ensure that the
//                    // solution of this plastification step does not explode backwards to the
//                    // desired loading path. This may happen when the yield function is too
//                    // nonlinear and may present a second solution in the negative range of
//                    // the desired loading path. By setting the initial guess equal to 1.0
//                    // (that means the whole loading path is plastic - true only in perfect
//                    // plastic materials) the code attempts to bias the solver towards the
//                    // physical meaningful solution.
//                }
//            }
//            count++;
//            
//            multipl_FAD = multiplN;
//            multipl_FAD.diff(0,nVars);
//            fN.EpsT().CopyTo(stateAtYield_FAD.m_eps_t);
//            stateAtYield_FAD.m_eps_t.Add(deltaEps_FAD, multipl_FAD);
//            
//            // this method computes sigma, A(alpha) substracting E_p from epstotal
//            ComputePlasticVars(stateAtYield_FAD,
//                               sigma_FAD,
//                               A_FAD);
//            
//            // compute the value of the yield functions
//            fYC.Compute(sigma_FAD, A_FAD, phi_FAD, 0);
//            
//            
//        }while(fabs (phi_FAD[i].val()) > fResTol && count < fMaxNewton);
//        
//#ifdef PZ_LOG
//        {
//            if(count >= fMaxNewton )
//            {
//                std::stringstream sout;
//                sout << "*** FindPointAtYield *** multiplication factor= ";
//                sout << multiplN << " such that yield of function ";
//                sout << i << " = " << phi_FAD[i].val() << "; ";
//                sout << "\n#### Truncated Newton after "
//                << fMaxNewton << " steps with phi[" << i << "] = "
//                << TPZExtractVal::val(phi_FAD[i]) << "####. Results are unpredictable";
//                LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
//                //plasticIntegrOutput = 1;
//            }
//        }
//#endif
//#ifdef PZ_LOG
//        {
//            std::stringstream sout;
//            sout << "*** FindPointAtYield *** multiplication factor= ";
//            sout << multiplN << " such that yield of function ";
//            sout << i << " = " << phi_FAD[i].val();
//            if(count >= fMaxNewton )
//            {
//                sout << "\n#### Truncated Newton after "
//                << fMaxNewton << " steps with phi[" << i << "] = "
//                << TPZExtractVal::val(phi_FAD[i])
//                << "####. It appears in such load path the plastic yield solution was found in the opposite direction."
//                << "\nPlease check the other(s) yield function(s) to guarantee at least one of them is plastifying within the proposed load path.";
//                LOGPZ_WARN(logger,sout.str().c_str());
//                LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
//                
//            }else{
//                LOGPZ_DEBUG(logger,sout.str().c_str());
//            }
//        }
//#endif
//        
//        if(multiplN < minMultipl)minMultipl = multiplN;
//    }
//    
//    if(minMultipl < - fResTol)
//    {
//#ifdef PZ_LOG
//        {
//            std::stringstream sout;
//            sout << "*** FindPointAtYield *** Ignoring deltaStrainMultiplier = " << minMultipl
//            << " and setting it to zero.\n\t###### WARN: EpsTotalN isn't elastic!! ######";
//            LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
//        }
//#endif
//#ifdef PZ_LOG
//        {
//            std::stringstream sout;
//            sout << "*** FindPointAtYield *** Ignoring deltaStrainMultiplier = " << minMultipl
//            << " and setting it to zero.\n\t###### WARN: EpsTotalN isn't elastic!! ######";
//            LOGPZ_WARN(logger,sout.str().c_str());
//        }
//#endif
//    }
//    
//    
//    if(minMultipl < 0.)minMultipl = 0.; //avoiding rounding errors
//    
//    stateAtYield = fN;
//    stateAtYield.m_eps_t.Add(deltaEps, minMultipl);
//    
//#ifdef PZ_LOG
//    {
//        std::stringstream sout;
//        sout << "<<< FindPointAtYield *** Exiting Method with deltaStrainMultiplier = " << minMultipl;
//        LOGPZ_DEBUG(logger,sout.str().c_str());
//        if(plasticIntegrOutput)LOGPZ_DEBUG(plasticIntegrLogger,sout.str().c_str());
//    }
//#endif
//    
//    return minMultipl;
//}

template <class YC_t, class TF_t, class ER_t>
REAL TPZPlasticStep<YC_t, TF_t, ER_t>::FindPointAtYield(
        const TPZTensor<REAL> &epsTotalNp1,
        TPZPlasticState<REAL> &stateAtYield)const {

#ifdef PZ_LOG

    int plasticIntegrOutput;

    {
        std::stringstream sout;
        sout << ">>> FindPointAtYield ***";
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif

    const int nVars = 1;
    typedef TFad<nVars, REAL> TFAD_One;
    TPZTensor<REAL> deltaEps;
    TFAD_One A_FAD, multipl_FAD;
    TPZManVector<TFAD_One, 10> phi_FAD(YC_t::NYield);


    // this will be the result of the method (it is always a number between 0 and 1)
    REAL multiplN, minMultipl = 1.;

    TPZPlasticState<TFAD_One> stateAtYield_FAD;

    stateAtYield.CopyTo(stateAtYield_FAD);

    TPZTensor<STATE> sigma;
    STATE A;
    TPZManVector<STATE> phi(YC_t::NYield);

    // this method computes sigma, A(alpha) substracting E_p from epstotal
    ComputePlasticVars(stateAtYield, sigma, A);
    // compute the value of the yield functions
    fYC.Compute(sigma, A, phi, 0);


    // delta epsilon contains epstotal_N+1 - epsTotal_N
    epsTotalNp1.CopyTo(deltaEps);
    deltaEps.Add(fN.EpsT(), -1.);

    int i, nyield = YC_t::NYield;


    for (i = 0; i < nyield; i++) { // searching for the multiplier for each yield surface

        multiplN = 1.; // the plastic multiplier is likely to be zero
        // because of possible previous plastifications,
        // although a null value may lead to no derivative
        // or numerical instabilities when deltaEpsT=0

        TPZTensor<TFAD_One> deltaEps_FAD, sigma_FAD;
        deltaEps.CopyTo(deltaEps_FAD);

        // this will set the derivative in the value of stateatyield
        TFAD_One mult(multiplN, 0);
        deltaEps_FAD *= mult;
        stateAtYield.CopyTo(stateAtYield_FAD);
        stateAtYield_FAD.m_eps_t += deltaEps_FAD;

        // stateAtYield is at the next point

        // this method computes sigma, A(alpha) substracting E_p from epstotal
        ComputePlasticVars(stateAtYield_FAD, sigma_FAD, A_FAD);
        // compute the value of the yield functions
        fYC.Compute(sigma_FAD, A_FAD, phi_FAD, 0);

#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << "State at yield " << stateAtYield_FAD << std::endl;
            sout << "Sigma " << sigma_FAD << std::endl;
            sout << "phi " << phi_FAD << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif

        // the next value will be elastic. The multiplier cannot be diminished
        // if the value is less or equal to zero the multiplier is one
        if (phi_FAD[i].val() < fResTol) {
            continue;
        }

        // if the point was already on the plastic surface, the multiplier is zero
        if (phi[i] > -fResTol) {
            minMultipl = 0.;
            // nothing else to do, the multiplier can t be smaller
            return minMultipl;
        }

        // we have to find a scalar value which brings the point on the surface
        while (fabs(phi_FAD[i].val()) > fResTol) {

            REAL phi_prev = phi_FAD[i].val();
            deltaEps.CopyTo(deltaEps_FAD);

            // this will set the derivative in the value of stateatyield
            TFAD_One mult(multiplN, 0);
            deltaEps_FAD *= mult;
            stateAtYield.CopyTo(stateAtYield_FAD);
            stateAtYield_FAD.m_eps_t += deltaEps_FAD;

            // this method computes sigma, A(alpha) substracting E_p from epstotal
            ComputePlasticVars(stateAtYield_FAD, sigma_FAD, A_FAD);
            // compute the value of the yield functions
            fYC.Compute(sigma_FAD, A_FAD, phi_FAD, 0);

#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "State at yield " << stateAtYield_FAD << std::endl;
                sout << "Sigma " << sigma_FAD << std::endl;
                sout << "phi " << phi_FAD << std::endl;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            if (fabs(phi_FAD[i].val()) > fabs(phi_prev)) {
                DebugStop();
            }

            REAL derX = phi_FAD[i].dx(0);
            if (fabs(derX) < 1.e-8) derX += 1.e-8;
            REAL delmult = -phi_FAD[i].val() / derX; // avoiding division by zero
            if (multiplN == 1.0 && delmult > 0.) {
                // there is nothing to search further for this yield surface
                break; // The step is known as plastic in advance
            }
            multiplN += delmult;
            if (multiplN > 1.) multiplN = 1.;
        }

#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "*** FindPointAtYield *** multiplication factor= ";
            sout << multiplN << " such that yield of function ";
            sout << i << " = " << phi_FAD[i].val();
            // TODO @gustavobat: I created the variable 'count' here to be able to compile the code, please fix this.
            int count = 0;
            if (count >= fMaxNewton) {
                sout << "\n#### Truncated Newton after "
                        << fMaxNewton << " steps with phi[" << i << "] = "
                        << TPZExtractVal::val(phi_FAD[i])
                        << "####. It appears in such load path the plastic yield solution was found in the opposite direction."
                        << "\nPlease check the other(s) yield function(s) to guarantee at least one of them is plastifying within the proposed load path.";
                LOGPZ_WARN(logger, sout.str().c_str());
                LOGPZ_WARN(plasticIntegrLogger, sout.str().c_str());

            } else {
                LOGPZ_DEBUG(logger, sout.str().c_str());
            }
        }
#endif

        if (multiplN < minMultipl)minMultipl = multiplN;
    }

    if (minMultipl < fResTol) {
        minMultipl = 0.;
#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "*** FindPointAtYield *** Ignoring deltaStrainMultiplier = " << minMultipl
                    << " and setting it to zero.\n\t###### WARN: EpsTotalN isn't elastic!! ######";
            LOGPZ_WARN(plasticIntegrLogger, sout.str().c_str());
        }
#endif
#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "*** FindPointAtYield *** Ignoring deltaStrainMultiplier = " << minMultipl
                    << " and setting it to zero.\n\t###### WARN: EpsTotalN isn't elastic!! ######";
            LOGPZ_WARN(logger, sout.str().c_str());
        }
#endif
    }


    stateAtYield = fN;
    stateAtYield.m_eps_t.Add(deltaEps, minMultipl);

#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "<<< FindPointAtYield *** Exiting Method with deltaStrainMultiplier = " << minMultipl;
        LOGPZ_DEBUG(logger, sout.str().c_str());
        if (plasticIntegrOutput)LOGPZ_DEBUG(plasticIntegrLogger, sout.str().c_str());
    }
#endif

    return minMultipl;
}

/**
 * @brief Proposes an update to the plastic variables and estimates the relative error
 * comitted in this update. Neither internal variable are used nor changed.
 * In the Np1 variables, EpsT is imposed [in] and the Alpha and EpsP are evaluated.
 * It returns 1 if suceeded of 0 if tolerance not achieved.
 * @param N [in] Plastic state variables at time N
 * @param Np1 [in/out] Plastic state variables at time N+1
 * @param delGamma [in/out] plastic multipliers
 */
template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::InitialGuess(
        const TPZPlasticState<REAL> &N,
        TPZPlasticState<REAL> &Np1,
        TPZVec<REAL> &delGamma,
        TPZVec<int> &valideqs
        ) {

}

template <class YC_t, class TF_t, class ER_t>
int TPZPlasticStep<YC_t, TF_t, ER_t>::PlasticLoop(
        const TPZPlasticState<REAL> &NN,
        TPZPlasticState<REAL> &Np1,
        TPZVec<REAL> &delGamma,
        REAL &normEpsPErr,
        REAL &lambda,
        TPZVec<int> & validEqs) {

    // 6 variables for tensor and others: alpha and deltaGamma
    const int nVars = 7 + YC_t::NYield;
    int i;
    const REAL Tol = fResTol;

#ifdef PZ_LOG
    if (pointloadconfig.isDebugEnabled()) {
        std::stringstream sout;
        sout << "alpha = " << NN.m_hardening << std::endl;
        TPZFNMatrix<6, REAL> Ep(NN.m_eps_p), EpsT(Np1.m_eps_t);
        TPZTensor<REAL> sigmaT, deformElast;
        Ep.Print("Ep = ", sout, EMathematicaInput);
        EpsT.Print("Etotal = ", sout, EMathematicaInput);
        deformElast = Np1.m_eps_t;
        deformElast.Add(NN.m_eps_p, -1.);
        fER.ComputeStress(deformElast, sigmaT);
        TPZFNMatrix<6, REAL> sigma(sigmaT);
        sigma.Print("sigmaTrial = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(pointloadconfig, sout.str())
    }
#endif
    InitialGuess(NN, Np1, delGamma, validEqs);
#ifdef PZ_LOG
    {
        std::stringstream sout; //1, sout2;
        sout << " >>> PlasticLoop ***\n";
        sout << "Np1 = " << Np1
                << "\ndelGamma = " << delGamma
                << "\nlambda = " << lambda
                << "\nValid eqs = " << validEqs;
        //LOGPZ_DEBUG(logger,sout.str().c_str());
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif

    typedef TFad<nVars, REAL> TFAD;

    TPZPlasticState<TFAD> Np1_FAD;
    TPZTensor<TFAD> sigmaNp1_FAD;
    TPZManVector<TFAD, nVars> epsRes_FAD(nVars), delGamma_FAD(YC_t::NYield);
    TFAD phiRes_FAD;
    // ResVal to hold the residual vector 
    // Sol will contain the values of the unknown values
    TPZFNMatrix<nVars> ResVal(nVars, 1, 0.), Sol(nVars, 1, 0.);
    TPZFNMatrix<nVars * nVars> tangent(nVars, nVars, 0.); // Jacobian matrix
    REAL resnorm = 0.;

    TPZManVector<STATE, nVars> pivots(nVars);


    int countReset = 0;
    for (int i = 0; i < validEqs.size(); i++) {
        if (validEqs[i] != 0) {
            countReset = 1;
        }
    }
    int countNewton = 0;
    int switchvalid = 0;

    do {//while(RemoveInvalidEqs(epsRes_FAD));
        switchvalid++;
        // verifying if it is necessary to impose post-peak material behavior
        TPZTensor<REAL> sigmaGuess;
        REAL AGuess;
        ComputePlasticVars<REAL>(Np1, sigmaGuess, AGuess);
        fYC.SetYieldStatusMode(sigmaGuess, AGuess);

        InitializePlasticFAD(Np1, delGamma, Np1_FAD, delGamma_FAD);
#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "Before plastic residual\n";
            sout << "Np1_FAD\n" << Np1_FAD << std::endl;
            sout << "delGamma_FAD " << delGamma_FAD << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif

        PlasticResidual<REAL, TFAD>(NN, Np1_FAD, delGamma_FAD, epsRes_FAD, normEpsPErr);

#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "After Plastic residual\nepsRes_FAD = \n" << epsRes_FAD << std::endl;
            sout << "delGamma_FAD " << delGamma_FAD << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif

        if (countReset == 0)InitializeValidEqs(epsRes_FAD, validEqs);

#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "*** PlasticLoop *** PlasticLoop main loop with "
                    << countReset << "-th valid equation set: {"
                    << validEqs << "}.";
            LOGPZ_DEBUG(logger, sout.str().c_str());
        }
#endif


        //      REAL LineSearch(const TPZFMatrix &Wn, TPZFMatrix DeltaW, TPZFMatrix &NextW, REAL tol, int niter);

        ExtractTangent(epsRes_FAD, ResVal, resnorm, tangent, validEqs, pivots, 1/*precond*/, 1/*ResetUnvalidEqs*/);
#ifdef PZ_LOG
        {
            std::stringstream sout;
            tangent.Print("A = ", sout, EMathematicaInput);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        countNewton = 0;

        do {

            countNewton++;

            TPZFNMatrix<nVars * nVars> *matc = new TPZFNMatrix<nVars * nVars>(nVars, nVars);
            *matc = tangent;
            TPZStepSolver<REAL> st(matc);
            st.SetDirect(ELU);

            //           cout << " RESVAL " <<endl;
            //           cout << ResVal << endl;


            //  ofstream fileout("rigidez.nb");
            //  tangent.Print("Rigidez = ", fileout, EMathematicaInput);
            //       ofstream fileout2("ResVal.nb");
            // ResVal.Print("ResVal = ", fileout2, EMathematicaInput);
            //  cout << tangent << endl;
            //   cout << ResVal << endl;


            // invert the tangent matrix and put the correction in the Sol variable
            st.Solve(ResVal, Sol, 0);

            // update the independent variables (their partial derivatives remain unchanged)
            lambda = UpdatePlasticVars(NN, Np1_FAD, delGamma_FAD, epsRes_FAD, Sol, validEqs);


            if (countNewton > 8) {
                std::cout << "countNewton = " << countNewton << " resnorm " << resnorm << " tolerance " << fResTol << "\n";
            }
            // recompute the residual
            PlasticResidual<REAL, TFAD>(NN, Np1_FAD, delGamma_FAD, epsRes_FAD, normEpsPErr);

#ifdef PZ_LOG
            {
                std::stringstream sout;
                sout << "Plastic residual " << epsRes_FAD << std::endl;
                sout << "delGamma_FAD " << delGamma_FAD << std::endl;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            // extract the values of the residual vector
            ExtractTangent(epsRes_FAD, ResVal, resnorm, tangent, validEqs, pivots, 1, 1);

#ifdef PZ_LOG
            {
                std::stringstream sout1, sout2;
                sout1 << "*** PlasticLoop *** Newton's " << countNewton
                        << "-th iteration with residual norm = " << resnorm;
                //sout1 << "\nTangent Matrix:\n" << tangent
                // << "\nResidual:\n" << ResVal;
                // LOGPZ_DEBUG(logger,sout1.str().c_str());
            }
#endif

        } while (resnorm > Tol && countNewton < fMaxNewton);

#ifdef PZ_LOG
        {

            std::stringstream sout;
            sout << "*** PlasticLoop *** Exiting Newton's scheme after " << countNewton
                    << " steps and with residual norm = " << resnorm;
            if (resnorm > Tol) {
                sout << "\n###### Max Newton steps achieved - Truncating Newton ######";
                LOGPZ_WARN(plasticIntegrLogger, sout.str().c_str());
            } else {
                LOGPZ_DEBUG(plasticIntegrLogger, sout.str().c_str());
            }
        }
#endif
#ifdef PZ_LOG
        {
            std::stringstream sout1, sout2;
            sout1 << "*** PlasticLoop *** Exiting Newton's scheme after " << countNewton
                    << " steps and with residual norm = " << resnorm;
            LOGPZ_DEBUG(logger, sout1.str().c_str());
            if (resnorm > Tol) {
                sout2 << "\n###### Max Newton steps achieved - Truncating Newton ######";
                LOGPZ_WARN(logger, sout2.str().c_str());
            }
        }
#endif

        countReset++;

    } while (switchvalid < 3 && RemoveInvalidEqs(delGamma_FAD, epsRes_FAD, validEqs));

    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    // updating output variables
    //for(i=0; i<6; i++) Np1.m_eps_p.fData[i] = Np1_FAD.EpsP().fData[i].val();
    //Np1.m_hardening = Np1_FAD.VolHardening().val();
#ifdef PZDEBUG
    if (switchvalid == 3) {
        std::cout << __FUNCTION__ << " should stop switchvalid\n";
    }
#endif
    Np1_FAD.CopyTo(Np1);

    for (i = 0; i < YC_t::NYield; i++)
        delGamma[i] = delGamma_FAD[i].val();



    //    for(int i=0;i<6;i++)
    //    {
    //
    //        fOutfile << "{ {" << Np1.m_eps_t.XX() << "," << Np1.m_eps_t.XY() << "," << Np1.m_eps_t.XZ()
    //        <<"},{"<<Np1.m_eps_t.XY()<<","<<Np1.m_eps_t.YY()<<","<<Np1.m_eps_t.YZ()
    //        <<"},{"<<Np1.m_eps_t.XZ()<<","<<Np1.m_eps_t.YZ() << ","<<Np1.m_eps_t.ZZ()<<"} } ";
    //    }



#ifdef PZ_LOG

    {
        std::stringstream sout1, sout2;
        sout1 << "<<< PlasticLoop *** Exiting Method";
        sout2 << "\nNp1 << " << Np1
                << "\ndelGamma = " << delGamma
                << "\nValidEqs = {" << validEqs << "}";
        LOGPZ_DEBUG(loggerx, sout1.str().c_str());
        LOGPZ_DEBUG(loggerx, sout2.str().c_str());
    }

#endif
#ifdef PZ_LOG
    if (pointloadconfig.isDebugEnabled()) {
        std::stringstream sout;
        sout << "alphanp1 = " << Np1.m_hardening << std::endl;
        TPZFNMatrix<6, REAL> Ep(Np1.m_eps_p);
        Ep.Print("Epnp1 = ", sout, EMathematicaInput);
        sout << "delGamma = {" << delGamma << "};\n";
        LOGPZ_DEBUG(pointloadconfig, sout.str())
    }
#endif

    return resnorm <= Tol;
}

template <class YC_t, class TF_t, class ER_t>
int TPZPlasticStep<YC_t, TF_t, ER_t>::PlasticIntegrate(
        const TPZPlasticState<REAL> &N,
        TPZPlasticState<REAL> &Np1,
        const REAL &TolEpsPErr) {
    int succeeded;

    TPZManVector<REAL, YC_t::NYield> delGamma(YC_t::NYield, 0.);
    TPZManVector<int, YC_t::NYield> validEqs(YC_t::NYield, 0);

    REAL normEpsPErr = 0.;
    int counter = 0;
    REAL lambda = 0.;

    succeeded = PlasticLoop(N, Np1, delGamma, normEpsPErr, lambda, validEqs);

    if (normEpsPErr < TolEpsPErr && succeeded) {

        PushPlasticMem(Np1, 1., lambda, delGamma, validEqs, fYC.GetForceYield());

        return 1;
    }

    int discarded = 1;

    // Substepping state variables
    TPZPlasticState<REAL> Nk(N),
            Nkp1;

    TPZTensor<REAL> deltaEpsTotal(Np1.EpsT());
    deltaEpsTotal.Add(N.EpsT(), -1.);

    REAL k = 0., q = 1., kp1 = 0., multipl;

    ofstream outfile("integrate.txt");
    ofstream outfile2("integrate2.txt");

    while (k < 1.) {


        multipl = 0.95 * pow(TolEpsPErr / normEpsPErr, 0.5);

        if (multipl < 0.1) multipl = 0.1;
        if (multipl > 10.) multipl = 10.;

        if (!succeeded)multipl = 0.5;

        q *= multipl;

        if (q < fMinStepSize)q = fMinStepSize;

        kp1 = min((REAL) 1.0, k + q);
        q = kp1 - k; // needed when k+q exceeds 1
        Nkp1.m_eps_t = N.EpsT();
        Nkp1.m_eps_t.Add(deltaEpsTotal, kp1);
        Nkp1.m_hardening = Nk.VolHardening();
        Nkp1.m_eps_p = Nk.EpsP();
        for (int i = 0; i < YC_t::NYield; i++)delGamma[i] = 0.;

#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "Nkp1 = " << Nkp1 << endl;
            sout << "Nk = " << Nk << endl;
            sout << "k = " << k << endl;
            sout << "normEpsPErr = " << normEpsPErr << endl;
            sout << "delGamma = " << delGamma << endl;
            sout << "lambda = " << lambda << endl;
            LOGPZ_DEBUG(loggerx, sout.str())
        }
#endif


        succeeded = PlasticLoop(Nk, Nkp1, delGamma, normEpsPErr, lambda, validEqs);


        if ((normEpsPErr < TolEpsPErr && succeeded) || kp1 - k < fMinStepSize * 1.001) // 1.001 because of rounding errors
        {

            PushPlasticMem(Nkp1, kp1, lambda, delGamma, validEqs, fYC.GetForceYield());
            counter++;

            // the k-th integration step respects the estimated tolerance
            // proceeding with time evolution...
            Nk = Nkp1;
            k = kp1;

            outfile << "\n" << counter << " " << k;
        } else {
            discarded++;
        }// otherwise the answer isn't accepted and the next q will be
        // recomputed in order to lie within the desired tolerance.
        // If the method works fine this situation should rarely happen.
        REAL razaoerros = normEpsPErr / TolEpsPErr;
        outfile2 << "\n" << counter << " " << razaoerros;

    }

    return counter;
}

template <class YC_t, class TF_t, class ER_t>
template <class T>
int TPZPlasticStep<YC_t, TF_t, ER_t>::InitializeValidEqs(TPZVec<T> &res_T, TPZVec<int> & validEqs) {
#ifdef PZ_LOG
    {
        int i, n = res_T.NElements();
        std::stringstream sout;
        sout << ">>> InitializeValidEqs *** Res = ";
        for (i = 0; i < n; i++)sout << TPZExtractVal::val(res_T[i]) << " ";
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif
    const int nyield = YC_t::NYield;
    validEqs.Resize(nyield);
    int i, count = 0;
    for (i = 0; i < nyield; i++) {
        validEqs[i] = 0;
        if (TPZExtractVal::val(res_T[i + 7]) > 0.) {
            validEqs[i] = 1;
            count++;
        }
    }

    if (nyield == 1) // if there exists only one eqn it must be active
    {
        validEqs[0] = 1;
        count++;
    }

    return count;
}

template <class YC_t, class TF_t, class ER_t>
template <class T>
int TPZPlasticStep<YC_t, TF_t, ER_t>::RemoveInvalidEqs(TPZVec<T> & delGamma_T, TPZVec<T> & res_T, TPZVec<int> & validEqs) {
    const int nyield = YC_t::NYield;
    validEqs.Resize(nyield);

    if (nyield == 1)return 0; // remove Invalid Eqs should not invalidate the unique equation

    int i, count = 0;

#ifdef PZ_LOG
    {
        int i;
        std::stringstream sout;
        sout << ">>> RemoveInvalidEqs *** delGamma = ";
        for (i = 0; i < nyield; i++)sout << TPZExtractVal::val(delGamma_T[i]) << " ";
        sout << " phi = ";
        for (i = 0; i < nyield; i++)sout << TPZExtractVal::val(res_T[i + 7]) << " ";
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif

    int BoolDelGamma,
            BoolValidEq,
            BoolRes;

    for (i = 0; i < nyield; i++) {
        BoolDelGamma = TPZExtractVal::val(delGamma_T[i]) > 0.; // deltaGamma is valid
        BoolValidEq = validEqs[i]; // equation is already in use
        BoolRes = TPZExtractVal::val(res_T[i + 7]) > fResTol; // equation indicates plastification

        switch (BoolRes) {
            case (true):
                // the lines below unfortunately led to instabilities in the integration process
                //validEqs[i] = 1; // if the equation indicates plastification then it is necessary to involve it in the process
                // the equation was not computed and the plasticity function is positive
                if (!BoolValidEq) {
                    //                    validEqs[i] = 1;
                    //                    count++;
                }
                break;
            case (false):
                // if the equation does not indicate plastification but is
                // already involved in the plastification process and
                // has a valid deltaGamma then keep it in the integration process - do nothing

                // if it is involved in the plastification process but has an
                // invalid deltagamma - invalidate it
                if (BoolValidEq && !BoolDelGamma) {
                    validEqs[i] = 0;
                    count++;
                }

                // if it isn't invoved in the plastification process - delgamma is meaningless and no
                // testing should be performed.
                break;
        }
    }

    return count;
}

template <class YC_t, class TF_t, class ER_t>
template<class T1, class T_VECTOR, class T_MATRIX>//T1:input residual fad type (FAD FAD or FAD), T_MATRIX: output matrix &vector of type (FAD or REAL, respectvly)
void TPZPlasticStep<YC_t, TF_t, ER_t>::ExtractTangent(
        const TPZVec<T1> & epsRes_FAD,
        T_VECTOR & ResVal, //TPZFMatrix for the T1=fad<real> type
        REAL & resnorm, // REAL
        T_MATRIX & tangent, // TPZFMatrix for T1=fad<real> type
        TPZVec<int> & validEqs,
        TPZVec<REAL> &pivots,
        const int precond,
        const int resetInvalidEqs) {

    const int nyield = YC_t::NYield;
    const int nVars = 7 + nyield; // 6 stresses + alpha + 2 delGamma
    int i, j;
    int ncols = epsRes_FAD[0].size();

    // extract the partial derivatives to form the tangent matrix
    for (i = 0; i < nVars; i++) {
        for (j = 0; j < ncols; j++) {
            tangent(i, j) = epsRes_FAD[i].dx(j);
        }
        ResVal(i, 0) = epsRes_FAD[i].val();
    }

    //  cout << "epsRes_FAD[i] "<<endl;
    //  cout << epsRes_FAD <<endl;
    //
    //  cout << "ResVal"<<endl;
    //  cout << ResVal<<endl;
    //
    //  cout << " TANGENTE " <<endl;
    //  cout << tangent <<endl;

    // reseting the equations related to the invalid yield surfaces
    if (resetInvalidEqs) {
        for (i = 0; i < nyield; i++) {
            if (validEqs[i] == 0) {
                for (j = 0; j < ncols; j++) {
                    tangent(i + 7, j) = 0.;
                }
                if (i + 7 < ncols) {
                    tangent(i + 7, i + 7) = 1.;
                }

                ResVal(i + 7, 0) = 0.;
            }
        }
    }

    // updating residual norm
    // before preconditioning such that it keeps its physical meaning
    //cout << "ResVal"<<endl;
    //  cout << ResVal<<endl;
    resnorm = 0.;
    for (i = 0; i < nVars; i++) {
        resnorm += pow(TPZExtractVal::val(ResVal(i, 0)), 2.);
    }
    resnorm = sqrt(resnorm);

    if (precond && nVars != ncols) {
        DebugStop();
    }
    //cout << "resnorm"<<endl;
    //  cout << resnorm<<endl;
    // preconditioning the matrix/residual system
    if (precond) {
        for (i = 0; i < 7; i++) {
            REAL pivot = 0., elem;
            j = i;
            elem = TPZExtractVal::val(tangent(i, j));
            if (fabs(pivot) < fabs(elem))pivot = elem;
            pivots[i] = pivot;
            //pivot = fabs(TPZExtractVal::val(pivot) ) < fabs(TPZExtractVal::val(tangent(i,j) ) ) ? tangent(i,j) : pivot;

            for (j = 0; j < nVars; j++)
                tangent(i, j) = tangent(i, j) / pivot;
            ResVal(i, 0) = ResVal(i, 0) / pivot;
        }

        for (; i < nVars; i++) {
            REAL pivot = 0., elem;
            for (j = 0; j < nVars; j++) {
                elem = TPZExtractVal::val(tangent(i, j));
                if (fabs(pivot) < fabs(elem))pivot = elem;
                //pivot = fabs(TPZExtractVal::val(pivot) ) < fabs(TPZExtractVal::val(tangent(i,j) ) ) ? tangent(i,j) : pivot;
            }
            pivots[i] = pivot;
            for (j = 0; j < nVars; j++) {
                tangent(i, j) = tangent(i, j) / pivot;
            }
            ResVal(i, 0) = ResVal(i, 0) / pivot;
        }
    }
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << ">>> ExtractTangent *** Residual norm = " << resnorm
                << "; precond = " << precond
                << "; resetInvalidEqs = " << resetInvalidEqs
                << "; validEqs = {" << validEqs << "}";
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif
}

template <class YC_t, class TF_t, class ER_t>
template <class T1, class T2>
void TPZPlasticStep<YC_t, TF_t, ER_t>::PlasticResidual(
        const TPZPlasticState<T1> &N_T1,
        TPZPlasticState<T2> &Np1_T2,
        const TPZVec<T2> &delGamma_T2,
        TPZVec<T2> &res_T2,
        REAL &normEpsPErr,
        int silent)const {

    // a= 0 ->implicito
    // a= 0.5 -> ponto medio - segunda ordem
    const REAL a = 0.;

    // This function will be either called with template parameter T being REAL or FAD type
    // nyield indicates the number of yield functions
    const int nyield = YC_t::NYield;

    // variable to hold the value of the stress, elastic strain and residual corresponding to the plastic strain
    TPZTensor<T2> sigmaNp1_T2, epsResMidPt_T2;
    TPZTensor<T1> sigmaN_T1;
    // variable containing the error or deviation between solutions from Modified and
    // explicit Euler schemes
    TPZTensor<REAL> epsPErr;

    // vector holding the tensors of the N directions (usually derivative of the yield functions)
    TPZManVector<TPZTensor<T2>, nyield> NdirNp1_T2(nyield), NdirMidPt_T2(nyield);
    TPZManVector<TPZTensor<T1>, nyield> NdirN_T1(nyield);

    // vector holding the residual of the yield function equations
    TPZManVector<T2, nyield> phiRes_T2(nyield), HNp1_T2(nyield), HMidPt_T2(nyield);
    TPZManVector<T1, nyield> HN_T1(nyield);

    // residual of the equation corresponding to the damage variable
    T2 alphaRes_T2, ANp1_T2;
    T1 AN_T1;

    ComputePlasticVars<T1>(N_T1, sigmaN_T1, AN_T1);
    ComputePlasticVars<T2>(Np1_T2, sigmaNp1_T2, ANp1_T2);

    // Compute the values of the N directions (which will update the plastic strain)
    fYC.N(sigmaN_T1, AN_T1, NdirN_T1, 0);
    fYC.N(sigmaNp1_T2, ANp1_T2, NdirNp1_T2, 1);

    int i, j;
    for (i = 0; i < nyield; i++){
        for (j = 0; j < 6; j++){
            NdirMidPt_T2[i].fData[j] = ((NdirNp1_T2[i].fData[j]) * T2(1. - a)
                + (NdirN_T1[i].fData[j]) * T1(a));
        }
    }

    // Compute the value of H
    fYC.H(sigmaN_T1, AN_T1, HN_T1, 0);
    fYC.H(sigmaNp1_T2, ANp1_T2, HNp1_T2, 1);

    for (i = 0; i < nyield; i++) {
        HMidPt_T2[i] = (HNp1_T2[i] * T2(1. - a) + T2(HN_T1[i] * T1(a)));
    }

    //EpsRes = EpsilonPNp1 - m_eps_p - deltagamma * Ndir; // 6 eqs
    //for
    for (i = 0; i < nyield; i++) {
        NdirMidPt_T2[i].Multiply(delGamma_T2[i], T2(1.));
        epsResMidPt_T2.Add(NdirMidPt_T2[i], T2(-1.));

        NdirN_T1[i].Multiply(T1(TPZExtractVal::val(delGamma_T2[i])), T1(1.));
        for (j = 0; j < 6; j++) {
            epsPErr.fData[j] += TPZExtractVal::val(NdirN_T1[i]. fData[j])
                    - TPZExtractVal::val(NdirMidPt_T2[i].fData[j]);
        }
    }
    epsResMidPt_T2.Add(N_T1.EpsP(), T1(-1.));
    epsResMidPt_T2.Add(Np1_T2.EpsP(), T1(1.));



#ifdef PZ_LOG
    {
        std::stringstream sout;
        //sout << "\n NdirNp1_T2 = "<< NdirNp1_T2 <<endl;
        //sout << "\n NdirN_T1 = "<< NdirN_T1 <<endl;
        //LOGPZ_DEBUG(loggerx,sout.str().c_str());
    }
#endif

    // Explicit scheme relative error: estimated by the difference between the
    // first and second order schemes results.
    normEpsPErr = epsPErr.Norm();
    // The second order scheme error is estimated based on the explicit Euler
    // scheme. Its relative error measure coincides with the Explicit Scheme
    // error estimation.

    // alphaRes = alpha(n+1)-alpha(n)-Sum delGamma * H
    alphaRes_T2 = Np1_T2.VolHardening() - N_T1.VolHardening();

    T2 multiplier;

    fYC.AlphaMultiplier(Np1_T2.VolHardening(), multiplier);
    alphaRes_T2 *= multiplier;
    for (i = 0; i < nyield; i++) {
        alphaRes_T2 -= delGamma_T2[i] * HMidPt_T2[i]; // 1 eq
    }
    //    STATE multval = TPZExtractVal::val(multiplier);
    //    multval = sqrt(fabs(multval));
    //    alphaRes_T2 /= (T2)multval;

    // compute the values of the yield functions
    fYC.Compute(sigmaNp1_T2, ANp1_T2, phiRes_T2, 1); // nyield eq

    // transfer the values of the residuals to the residual vector
    // depending on the type T, the residual and its partial derivatives will have been computed
    for (i = 0; i < 6; i++) {
        res_T2[i] = epsResMidPt_T2.fData[i];
    }
    res_T2[i++] = alphaRes_T2;
    for (j = 0; j < nyield; j++) {
        res_T2[i++] = phiRes_T2[j];
    }

}

template <class YC_t, class TF_t, class ER_t>
template <class T1, class T2>
void TPZPlasticStep<YC_t, TF_t, ER_t>::PlasticResidualRK(
        const TPZPlasticState<T1> &N_T1,
        TPZPlasticState<T2> &Np1_T2,
        const TPZVec<T2> &delGamma_T2,
        TPZVec<T2> &res_T2,
        REAL &normEpsPErr,
        int silent)const {

    // This function will be either called with template parameter T being REAL or FAD type
    // nyield indicates the number of yield functions
    const int nyield = YC_t::NYield;

    // variable to hold the value of the stress, elastic strain and residual corresponding to the plastic strain
    TPZTensor<T2> sigmaNp1_T2, epsResMidPt_T2;
    TPZTensor<T1> sigmaN_T1;
    // variable containing the error or deviation between solutions from Modified and
    // explicit Euler schemes
    TPZTensor<REAL> epsPErr;

    // vector holding the tensors of the N directions (usually derivative of the yield functions)
    TPZManVector<TPZTensor<T2>, nyield> NdirNp1_T2(nyield), NdirMidPt_T2(nyield), K1N_T2(nyield), K2N_T2(nyield), K3N_T2(nyield), K4N_T2(nyield), KTOTALN_T2(nyield);
    TPZManVector<TPZTensor<T1>, nyield> NdirN_T1(nyield), NdirN_T1COPY(nyield);

    // vector holding the residual of the yield function equations
    TPZManVector<T2, nyield> phiRes_T2(nyield), HNp1_T2(nyield), HMidPt_T2(nyield), K1H_T2(nyield), K2H_T2(nyield), K3H_T2(nyield), K4H_T2(nyield), KTOTALH_T2(nyield);
    TPZManVector<T1, nyield> HN_T1(nyield);

    // residual of the equation corresponding to the damage variable
    T2 alphaRes_T2, ANp1_T2;
    T1 AN_T1;

    ComputePlasticVars<T1>(N_T1, sigmaN_T1, AN_T1);
    ComputePlasticVars<T2>(Np1_T2, sigmaNp1_T2, ANp1_T2);

    // Compute the values of the N directions (which will update the plastic strain)
    fYC.N(sigmaN_T1, AN_T1, NdirN_T1, 0);
    fYC.N(sigmaNp1_T2, ANp1_T2, NdirNp1_T2, 1);

    //k1 = epspN ,k2 = epspN +(N1/2)*k1, k3 = epspN+(N1/2)*k2, k3 = epspN+N1*k3
    int i, j;
    for (i = 0; i < nyield; i++) {
        for (j = 0; j < 6; j++) {
            K1N_T2[i].fData[j] = N_T1.m_eps_p.fData[j];
            K2N_T2[i].fData[j] = N_T1.m_eps_p.fData[j] + (NdirN_T1[i].fData[j]/*+NdirN_T2[i].fData[j]*/) * T2(0.5) * K1N_T2[i].fData[j];
            K3N_T2[i].fData[j] = N_T1.m_eps_p.fData[j] + (NdirN_T1[i].fData[j]/*+NdirN_T2[i].fData[j]*/) * T2(0.5) * K2N_T2[i].fData[j];
            K4N_T2[i].fData[j] = N_T1.m_eps_p.fData[j] + (NdirN_T1[i].fData[j]/*+NdirN_T2[i].fData[j]*/) * K3N_T2[i].fData[j];
            KTOTALN_T2[i].fData[j] = K1N_T2[i].fData[j] + T2(2.) * K2N_T2[i].fData[j] + T2(2.) * K3N_T2[i].fData[j] + K4N_T2[i].fData[j];

            NdirN_T1COPY[i].fData[j] = NdirN_T1[i].fData[j];

        }


    }


    // Compute the value of H
    fYC.H(sigmaN_T1, AN_T1, HN_T1, 0);
    fYC.H(sigmaNp1_T2, ANp1_T2, HNp1_T2, 1);

    for (i = 0; i < nyield; i++) {
        K1H_T2[i] = N_T1.m_hardening;
        K2H_T2[i] = N_T1.m_hardening + HN_T1[i] * T2(0.5) * K1H_T2[i];
        K3H_T2[i] = N_T1.m_hardening + HN_T1[i] * T2(0.5) * K2H_T2[i];
        K4H_T2[i] = N_T1.m_hardening + HN_T1[i] * K3H_T2[i];
        KTOTALH_T2[i] = K1H_T2[i]+(T2(2.) * K2H_T2[i])+(T2(2.) * K3H_T2[i]) + K2H_T2[i];

    }



    //EpsRes = EpsilonPNp1 - m_eps_p - ((deltagamma * Ndir)/6)*(k1+2k2+2k3+k4); // 6 eqs
    for (i = 0; i < nyield; i++) {
        //deltagamma * Ndir
        NdirN_T1COPY[i].Multiply(delGamma_T2[i], T2(1. / 6.));
        NdirN_T1COPY[i].Multiply(K1H_T2[i], T2(1.));
        epsResMidPt_T2.Add(NdirN_T1COPY[i], T2(-1.));

        NdirN_T1[i].Multiply(T1(TPZExtractVal::val(delGamma_T2[i])), T1(1.));
        for (j = 0; j < 6; j++)epsPErr.fData[j] += TPZExtractVal::val(NdirN_T1[i]. fData[j])
            - TPZExtractVal::val(NdirN_T1COPY[i].fData[j]);
    }
    //????
    epsResMidPt_T2.Add(N_T1.EpsP(), T1(-1.));
    epsResMidPt_T2.Add(Np1_T2.EpsP(), T1(1.));

    // Explicit scheme relative error: estimated by the difference between the
    // first and second order schemes results.
    normEpsPErr = epsPErr.Norm();
    // The second order scheme error is estimated based on the explicit Euler
    // scheme. Its relative error measure coincides with the Explicit Scheme
    // error estimation.

    // alphaRes = alpha(n+1)-alpha(n)-Sum delGamma * H
    alphaRes_T2 = Np1_T2.VolHardening() - N_T1.VolHardening();
    for (i = 0; i < nyield; i++) {
        alphaRes_T2 -= delGamma_T2[i] * HN_T1[i] * T2(1. / 6.) * KTOTALH_T2[i]; // 1 eq
    }

    // compute the values of the yield functions
    fYC.Compute(sigmaNp1_T2, ANp1_T2, phiRes_T2, 1); // nyield eq

    // transfer the values of the residuals to the residual vector
    // depending on the type T, the residual and its partial derivatives will have been computed
    for (i = 0; i < 6; i++) {
        res_T2[i] = epsResMidPt_T2.fData[i];
    }
    res_T2[i++] = alphaRes_T2;
    for (j = 0; j < nyield; j++) {
        res_T2[i++] = phiRes_T2[j];
    }


#ifdef PZ_LOG
    if (plasticIntegrLogger.isDebugEnabled()) {
        std::stringstream sout;
        /*  sout << "\n res_T2 = " << res_T2 <<endl;
         sout << "\n K1N_T2 = " << K1N_T2 <<endl;
         sout << "\n K2N_T2 = " << K2N_T2 <<endl;
         sout << "\n K3N_T2 = " << K3N_T2 << endl;
         sout << "\n K4N_T2 = " << K4N_T2 << endl;
         sout << "\n K1H_T2 = " << K1H_T2 <<endl;
         sout << "\n K2H_T2 = " << K2H_T2 <<endl;
         sout << "\n K3H_T2 = " << K3H_T2 << endl;
         sout << "\n K4H_T2 = " << K4H_T2 << endl;
         */
        sout << "\n epsResMidPt_T2 = " << epsResMidPt_T2 << endl;
        sout << "\n alphaRes_T2 = " << alphaRes_T2 << endl;
        sout << "\n NdirN_T1 = " << NdirN_T1 << endl;
        LOGPZ_DEBUG(plasticIntegrLogger, sout.str().c_str());
    }
#endif


}

template <class YC_t, class TF_t, class ER_t>
template<class T1, class T2, class TVECTOR>
REAL TPZPlasticStep<YC_t, TF_t, ER_t>::UpdatePlasticVars(
        const TPZPlasticState<T1> &N_T1,
        TPZPlasticState<T2> &Np1_T2,
        TPZVec<T2> &delGamma_T2,
        TPZVec<T2> &res_T2,
        TVECTOR & Sol_TVECTOR,
        TPZVec<int> & validEqs,
        int updateVars)const {
    const int nyield = YC_t::NYield;
    const int nVars = 7 + nyield;
    TPZManVector<REAL, nyield> delGamma(nyield);
    REAL sqrNormResN = 0, sqrNormResNp1, normEpsPErr;
    TPZManVector<REAL, nVars> res(nVars);
    TPZPlasticState<REAL> Np1, N;
    int i, j;

    Np1_T2.CopyTo(Np1);
    N_T1.CopyTo(N);

    for (i = 0; i < 6; i++) Np1.m_eps_p.fData[i] -= TPZExtractVal::val(Sol_TVECTOR(i));
    Np1.m_hardening -= TPZExtractVal::val(Sol_TVECTOR(i++));
    for (j = 0; j < nyield; j++)delGamma[j] = TPZExtractVal::val(delGamma_T2[j]) - TPZExtractVal::val(Sol_TVECTOR(i++));

    PlasticResidual<REAL, REAL>(N, Np1, delGamma, res, normEpsPErr, 1 /*silent*/);
    //    PlasticResidualRK<REAL, REAL>(N, Np1, delGamma, res, normEpsPErr, 1 /*silent*/);

    sqrNormResNp1 = pow(sdot(res, res), 0.5);

    const REAL k = 2.;
    REAL lambda = 1.; // ensuring that the lambda value will be ONE at the first step

    do {
        sqrNormResN = sqrNormResNp1;

        lambda /= k;

        Np1_T2.CopyTo(Np1);
        for (i = 0; i < 6; i++) Np1.m_eps_p.fData[i] -= lambda * TPZExtractVal::val(Sol_TVECTOR(i));
        Np1.m_hardening -= lambda * TPZExtractVal::val(Sol_TVECTOR(i++));
        for (j = 0; j < nyield; j++)delGamma[j] = TPZExtractVal::val(delGamma_T2[j]) - lambda * TPZExtractVal::val(Sol_TVECTOR(i++));

        PlasticResidual<REAL, REAL>(N, Np1, delGamma, res, normEpsPErr, 1 /*silent*/);
        //  PlasticResidualRK<REAL, REAL>(N, Np1, delGamma, res, normEpsPErr, 1 /*silent*/);

        // resetting invalid equations
        for (i = 0; i < nyield; i++)
            if (validEqs[i] == 0)
                res[i + 7] = 0.;

        sqrNormResNp1 = pow(sdot(res, res), 0.5);

    } while (sqrNormResNp1 < sqrNormResN && lambda >= fMinLambda); // ensuring that the step will be larger than fMinLambda
    // to avoid wandering within local minima.

    lambda *= k;

    if (lambda < 1.) {
#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "*** UpdatePlasticVars *** Line Search indicates lambda = " << lambda << " to ensure residual drop.";
            LOGPZ_DEBUG(logger, sout.str().c_str());
        }
#endif
    }

    if (!updateVars) return lambda;

    for (i = 0; i < 6; i++) Np1_T2.m_eps_p.fData[i] -= T2(lambda * Sol_TVECTOR(i));
    Np1_T2.m_hardening -= T2(lambda * Sol_TVECTOR(i++));
    for (j = 0; j < YC_t::NYield; j++) delGamma_T2[j] -= T2(lambda * Sol_TVECTOR(i++));

    return lambda;

}

template <class YC_t, class TF_t, class ER_t>
template<class T>
void TPZPlasticStep<YC_t, TF_t, ER_t>::InitializePlasticFAD(
        const TPZPlasticState<REAL> &state,
        const TPZVec<REAL> &delGamma,
        TPZPlasticState<T> &state_T,
        TPZVec<T> &delGamma_T,
        const int nVars)const {
    int i;
    int nVarsPlastic = 7 + YC_t::NYield;

    //  copying values
    state.CopyTo(state_T); // also initializing derivatives to null
    for (i = 0; i < YC_t::NYield; i++)delGamma_T[i] = delGamma[i];

    // Initialize the partial derivative values
    // the first 6 independent variables are the values of the plastic strains
    for (i = 0; i < 6; i++) state_T.m_eps_p.fData[i].diff(i, nVars);
    // the damage variable is the seventh variable

    state_T.m_hardening.diff(6, nVars);
    // the remaining variables are the yield function multipliers
    for (i = 7; i < nVarsPlastic; i++) delGamma_T[i - 7].diff(i, nVars);

}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ProcessLoad(const TPZTensor<REAL> &sigma, const EElastoPlastic ep) {
    const int nVars = 6;
    REAL resnorm;
    int i, k = 0;
    TPZFNMatrix<nVars * nVars> Dep_mat(nVars, nVars);
    TPZFNMatrix<nVars> residual_mat(nVars, 1),
            sol_mat(nVars, 1);

    TPZTensor<REAL> epsTotal(fN.m_eps_t), EEpsilon;

#ifdef PZ_LOG
    if (plasticIntegrLogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << ">>> ProcessLoad *** Evaluating Sigma to compute the resultant stress for the guess strain - Preparing starting tangent matrix for Newton's scheme";
        sout << "\n sigma << " << sigma;
        LOGPZ_DEBUG(plasticIntegrLogger, sout.str().c_str());
    }
#endif
#ifdef PZ_LOG
    {
        std::stringstream sout1, sout2;
        sout1 << ">>> ProcessLoad *** Evaluating Sigma to compute the resultant stress for the guess strain - Preparing starting tangent matrix for Newton's scheme";
        sout2 << "\n sigma << " << sigma;
        LOGPZ_DEBUG(logger, sout1.str().c_str());
        LOGPZ_DEBUG(logger, sout2.str().c_str());
    }
#endif

    //cout << "\nstarting ProcessStrain/ComputeDep";
    //cout.flush();
    /// DEBUG DEBUG
    /*
    TPZYCSandlerDimaggio *yc = (TPZYCSandlerDimaggio *)(&fYC);
    TPZYCSandlerDimaggioL *ycl = dynamic_cast<TPZYCSandlerDimaggioL *>(yc);
    TPZManVector<REAL> epsx(50),sigx(50),dsig(50),epsv(50),L(50),DepsVdL(50);
    
    for(int i=0; i<50; i++)
    {
        epsTotal.XX() = -0.025*i/50.;
        epsTotal.YY() = -0.025*i/50.;
        epsTotal.ZZ() = -0.025*i/50.;
        ProcessStrain(epsTotal, EAuto);
        ComputeDep(EEpsilon, Dep_mat);
        ApplyStrain(epsTotal);
        REAL diagdep = 0;
        for (int j=0; j< 6; j++) {
            diagdep += Dep_mat(j,0);
        }
        epsx[i] = epsTotal[0];
        sigx[i] = EEpsilon[0];
        dsig[i] = diagdep;
        if (ycl) {
            int n = this->fPlasticMem.size();
            L[i] = this->fPlasticMem[n-1].m_elastoplastic_state.m_hardening;
            yc->EpspFromL(L[i], epsv[i]);
            yc->DEpspDL(L[i], DepsVdL[i]);
        }
        else {
            REAL X;
            int n = this->fPlasticMem.size();
            epsv[i] = this->fPlasticMem[n-1].m_elastoplastic_state.m_hardening;
            yc->ComputeX(epsv[i], X);
            yc->SolveL(X, L[i]);
            yc->DEpspDL(L[i], DepsVdL[i]);
        }
        {
            std::stringstream sout;
            sout << "Sigma " << EEpsilon << std::endl;
            sout << "espTotal " << epsTotal << std::endl;
            sout << "D Sigma_x/Deps" << diagdep;
            Dep_mat.Print("depmat",sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
    }
    {
        std::stringstream sout;
        sout << "eps = { "<< epsx << "};\n";
        sout << "sig = { "<< sigx << "};\n";
        sout << "dsigdeps = { " << dsig << " };\n";
        sout << "epsv = { " << epsv << " };\n";
        sout << "L = { " << L << " };\n";
        sout << "depsvdl = { " << DepsVdL << " };\n";
        LOGPZ_DEBUG(logger, sout.str())
    }
    epsTotal = TPZTensor<REAL>();
     */
    ///
    // evaluating the plastic integration, stress tensor and jacobian
    //ProcessStrainNoSubIncrement(epsTotal,EAuto);
    ProcessStrain(epsTotal, EAuto);
    ComputeDep(EEpsilon, Dep_mat);

    //cout << "\nended ProcessStrain/ComputeDep";
    //cout.flush();

    resnorm = 0.;
    for (i = 0; i < nVars; i++)residual_mat(i, 0) = EEpsilon.fData[i] - sigma.fData[i];
    for (i = 0; i < nVars; i++)resnorm += pow(residual_mat(i, 0), 2.);
    resnorm = sqrt(resnorm);

    while (resnorm > fResTol && k < fMaxNewton) {
        k++;

        TPZFNMatrix<nVars * nVars> *matc = new TPZFNMatrix<nVars * nVars>(nVars, nVars);
        *matc = Dep_mat;

#ifdef PZ_LOG

        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            Dep_mat.Print("Derivative", sout);
            sout << "EEpsilon " << EEpsilon << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif

        //              cout << "\nNewton Method Step " << k << "\nDep=" << Dep_mat << endl << residual_mat;
        //              cout.flush();

        TPZStepSolver<REAL> st(matc);
        st.SetDirect(ELU);


        // invert the tangent matrix and put the correction in the Sol variable
        st.Solve(residual_mat, sol_mat, 0);

        //cout << "\n solve ended:" << sol_mat;
        //cout.flush();

        TPZTensor<REAL> epsTotalPrev(epsTotal);
        REAL scalefactor = 1.;
        REAL resnormprev = resnorm;

        do {
            for (i = 0; i < nVars; i++)epsTotal.fData[i] = epsTotalPrev.fData[i] - scalefactor * sol_mat(i, 0);

#ifdef PZ_LOG

            if (logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "Next epsTotal " << epsTotal << std::endl;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            //cout << "\nstarting ProcessStrain/ComputeDep";
            //cout.flush();

            // evaluating the plastic integration, stress tensor and jacobian
            //   ProcessStrainNoSubIncrement(epsTotal,ep);
            ProcessStrain(epsTotal, /*o original e ep e nao EAuto */ ep);
            ComputeDep(EEpsilon, Dep_mat);

            //cout << "\nended ProcessStrain/ComputeDep";
            //cout.flush();

            resnorm = 0.;
            for (i = 0; i < nVars; i++)residual_mat(i, 0) = EEpsilon.fData[i] - sigma.fData[i];
            for (i = 0; i < nVars; i++)resnorm += pow(residual_mat(i, 0), 2.);
            resnorm = sqrt(resnorm);

            scalefactor *= 0.5;
        } while (resnorm > resnormprev);
        //cout << "\nresidual = " << resnorm;

#ifdef PZ_LOG
        if (plasticIntegrLogger.isDebugEnabled()) {
            std::stringstream sout;
            sout << "*** ProcessLoad *** " << k << "-th iteration of Newton's scheme with residual = " << resnorm;
            LOGPZ_DEBUG(plasticIntegrLogger, sout.str().c_str());
        }
#endif
#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "*** ProcessLoad *** " << k << "-th iteration of Newton's scheme with residual = " << resnorm;
            LOGPZ_DEBUG(logger, sout.str().c_str());
        }
#endif

        if (k > fMaxNewton)cout << "\n*** ProcessLoad step " << k << " with res= " << resnorm;
    }

#ifdef PZ_LOG
    {
        std::stringstream sout;

        if (k > fMaxNewton) {
            sout << "<<< ProcessLoad *** Exiting Method with residual = " << resnorm
                    << " after " << k << " steps.";
            sout << "\n#### Truncated Newton ####. Results are unpredictable";
            LOGPZ_WARN(logger, sout.str().c_str());
            LOGPZ_WARN(plasticIntegrLogger, sout.str().c_str());
        } else {
            sout << "<<< ProcessLoad *** Exiting Method with residual = " << resnorm;
            LOGPZ_DEBUG(logger, sout.str().c_str());
            LOGPZ_DEBUG(plasticIntegrLogger, sout.str().c_str());
        }
    }
#endif

}

/**
 return the value of the yield function for the given strain
 */
template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::Phi_Internal(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const {
    TPZTensor<REAL> sigma;
    REAL A;

    TPZPlasticState<REAL> state(fN);
    state.m_eps_t = epsTotal;

    ComputePlasticVars<REAL>(state, sigma, A);

    fYC.Compute(sigma, A, phi, 0);

    return;
}

template <class YC_t, class TF_t, class ER_t>
template <int N>
void TPZPlasticStep<YC_t, TF_t, ER_t>::PushPlasticMem(
        const TPZPlasticState<REAL> & state,
        const REAL & k,
        const REAL & lambda,
        const TPZManVector<REAL, N> & delGamma,
        const TPZVec<int> & validEqs,
        const int forceYield) {
    TPZPlasticIntegrMem<REAL, YC_t::NYield> Mem(state, k, lambda, delGamma, validEqs, forceYield);
    fPlasticMem.Push(Mem);

}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyLoad_Internal(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal) {

#ifdef PZ_LOG
    if (plasticIntegrLogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << ">>> ApplyLoad_Internal ***"
                << " Imposed sigma << " << sigma;
        LOGPZ_DEBUG(plasticIntegrLogger, sout.str().c_str());
    }
#endif

#ifdef PZ_LOG
    if (pointloadconfig.isDebugEnabled()) {
        std::stringstream sout;
        sout << ">>> ApplyLoad_Internal ***"
                << " Imposed sigma << " << sigma;
        LOGPZ_DEBUG(pointloadconfig, sout.str().c_str());
    }
#endif

#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << ">>> ApplyLoad_Internal ***"
                << " Imposed sigma << " << sigma;
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif
    /*
     ProcessLoad(sigma, EAuto);
     
     int n = fPlasticMem.NElements();
     */

#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << ">>> ApplyLoad_Internal ***"
                << " Forcing Elastic behaviour";
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif

    ProcessLoad(sigma, EForceElastic);
    int n = fPlasticMem.NElements();

    if (!IsStrainElastic(fPlasticMem[n - 1].m_elastoplastic_state)) {
#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << ">>> ApplyLoad_Internal ***"
                    << " Forcing Plastic behaviour - Elastic attempt led to a final plastic state";
            LOGPZ_DEBUG(logger, sout.str().c_str());
        }
#endif
        ProcessLoad(sigma, EForcePlastic);
        n = fPlasticMem.NElements();
    }

    TPZPlasticStep<YC_t, TF_t, ER_t>::SetState_Internal(fPlasticMem[n - 1].m_elastoplastic_state);

    epsTotal = fN.m_eps_t;

#ifdef PZ_LOG
    {
        std::stringstream sout1, sout2;
        sout1 << "<<< ApplyLoad_Internal ***";
        sout2 << "\nOutput epsTotal << " << epsTotal;
        LOGPZ_DEBUG(logger, sout1.str().c_str());
        LOGPZ_DEBUG(logger, sout2.str().c_str());
    }
#endif

}

template <class YC_t, class TF_t, class ER_t>
REAL TPZPlasticStep<YC_t, TF_t, ER_t>::IntegrationOverview(TPZVec<REAL> & plastifLen) {
    int i, j, n = fPlasticMem.NElements();

    if (n <= 2) return 0;

    plastifLen.Fill(0.);
    REAL plasticLen = 0., deltaK;

    for (i = 2; i < n; i++) {
        deltaK = fPlasticMem[i].fK - fPlasticMem[i - 1].fK;
        for (j = 0; j < YC_t::NYield; j++)if (fPlasticMem[i].fValidEqs[j])plastifLen[j] += deltaK;
        plasticLen += deltaK;
    }

    return plasticLen;
}

template<class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::Write(TPZStream& buf, int withclassid) const {
    fYC.Write(buf, withclassid);
    fTFA.Write(buf, withclassid);
    fER.Write(buf, withclassid);
    buf.Write(&fResTol);
    buf.Write(&fIntegrTol);
    buf.Write(&fMaxNewton);
    buf.Write(&fMinLambda);
    buf.Write(&fMinStepSize);
    fN.Write(buf, withclassid);
    if (fPlasticMem.NElements() > 0) {
        DebugStop();
    }
    buf.Write(&fMaterialTensionSign);
    buf.Write(&fInterfaceTensionSign);
}

template<class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::Read(TPZStream& buf, void* context) {
    fYC.Read(buf, context);
    fTFA.Read(buf, context);
    fER.Read(buf, context);
    buf.Read(&fResTol);
    buf.Read(&fIntegrTol);
    buf.Read(&fMaxNewton);
    buf.Read(&fMinLambda);
    buf.Read(&fMinStepSize);
    fN.Read(buf, context);
    if (fPlasticMem.NElements() > 0) {
        DebugStop();
    }
    buf.Read(&fMaterialTensionSign);
    buf.Read(&fInterfaceTensionSign);
}

/// modify the elastic response. Needs to be reimplemented for each instantiation

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::SetElasticResponse(TPZElasticResponse &ER) {
    fER = ER;
}

template<>
void TPZPlasticStep<TPZYCLadeKim, TPZLadeKimThermoForceA, TPZLadeNelsonElasticResponse>::SetElasticResponse(TPZElasticResponse &ER) {
    DebugStop();
}

template <class YC_t, class TF_t, class ER_t>
TPZElasticResponse TPZPlasticStep<YC_t, TF_t, ER_t>::GetElasticResponse() const {
    return fER;
}

template<>
TPZElasticResponse TPZPlasticStep<TPZYCLadeKim, TPZLadeKimThermoForceA, TPZLadeNelsonElasticResponse>::GetElasticResponse() const {
    DebugStop();
    //Must return something
    TPZElasticResponse ret;
    return ret;
}



template <class YC_t, class TF_t, class ER_t>
TPZTensor<REAL> TPZPlasticStep<YC_t, TF_t, ER_t>::gRefDeform;


template class TPZPlasticStep<TPZYCVonMises, TPZThermoForceA, TPZElasticResponse>;
template class TPZPlasticStep<TPZYCTresca, TPZThermoForceA, TPZElasticResponse>;
template class TPZPlasticStep<TPZYCTrescaRegularized, TPZThermoForceA, TPZElasticResponse>;
template class TPZPlasticStep<TPZYCVonMisesCombTresca, TPZThermoForceA, TPZElasticResponse>;

template class TPZPlasticStep<TPZYCLadeKim, TPZLadeKimThermoForceA, TPZLadeNelsonElasticResponse>;

template class TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse>;
template class TPZPlasticStep<TPZYCRankine< TPZYCDruckerPrager >, TPZThermoForceA, TPZElasticResponse>;

#include "TPZYCMohrCoulomb.h"
#include  "TPZYCWillamWarnke.h"
#include  "TPZYCModifiedMohrCoulomb.h"
//#include  "TPZYCDruckerPragerBase.h"

template class TPZPlasticStep<TPZYCMohrCoulomb, TPZThermoForceA, TPZElasticResponse>;
template class TPZPlasticStep<TPZYCWillamWarnke, TPZThermoForceA, TPZElasticResponse>;
template class TPZPlasticStep<TPZYCModifiedMohrCoulomb, TPZThermoForceA, TPZElasticResponse>;
//template class TPZPlasticStep<TPZYCDruckerPragerBase< TPZYCMohrCoulomb >, TPZThermoForceA, TPZElasticResponse>;

template void TPZPlasticStep<TPZYCLadeKim, TPZLadeKimThermoForceA, TPZLadeNelsonElasticResponse>::
PlasticResidual<REAL, REAL>(TPZPlasticState<REAL> const &,
        TPZPlasticState<REAL> &,
        TPZVec<REAL> const&,
        TPZVec<REAL> &,
        REAL &, int) const;

template void TPZPlasticStep<TPZYCLadeKim, TPZLadeKimThermoForceA, TPZLadeNelsonElasticResponse>::
PlasticResidual<REAL, TFad<14, REAL> >(TPZPlasticState<REAL> const &,
        TPZPlasticState< TFad<14, REAL> > &,
        TPZVec<TFad<14, REAL> > const &,
        TPZVec<TFad<14, REAL> > &,
        REAL &, int) const;


//template void TPZPlasticStep<TPZYCLadeKim, TPZLadeKimThermoForceA, TPZLadeNelsonElasticResponse>::
//PlasticResidualRK<REAL, REAL>(TPZPlasticState<REAL> const &,
//                            TPZPlasticState<REAL> &,
//                            TPZVec<REAL> const&,
//                            TPZVec<REAL> &,
//                            REAL &, int)const;
//
//template void TPZPlasticStep<TPZYCLadeKim, TPZLadeKimThermoForceA, TPZLadeNelsonElasticResponse>::
//PlasticResidualRK<REAL, TFad<14,REAL> >(TPZPlasticState<REAL> const &,
//                                      TPZPlasticState< TFad<14,REAL> > &,
//                                      TPZVec<TFad<14,REAL> > const &,
//                                      TPZVec<TFad<14,REAL> > &,
//                                      REAL &, int)const;

/**
 * @brief Proposes an update to the plastic variables and estimates the relative error
 * comitted in this update. Neither internal variable are used nor changed.
 * In the Np1 variables, EpsT is imposed [in] and the Alpha and EpsP are evaluated.
 * It returns 1 if suceeded of 0 if tolerance not achieved.
 * @param N [in] Plastic state variables at time N
 * @param Np1 [in/out] Plastic state variables at time N+1
 * @param delGamma [in/out] plastic multipliers
 */
template <>
void TPZPlasticStep<TPZYCSandlerDimaggio, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>::InitialGuess(
        const TPZPlasticState<REAL> &N,
        TPZPlasticState<REAL> &Np1,
        TPZVec<REAL> &delGamma,
        TPZVec<int> &validEqs
        ) {
    TPZTensor<REAL> EpN = N.m_eps_p;
    TPZTensor<REAL> ETotal = Np1.m_eps_t;
    TPZTensor<REAL> ETrial = ETotal;
    TPZTensor<REAL> sigmaTrial;
    ETrial.Add(EpN, -1.);
    fER.ComputeStress(ETrial, sigmaTrial);
    TPZTensor<REAL> sigproj;
    fYC.InitialGuess(fER, N.m_hardening, sigmaTrial, Np1.m_hardening, delGamma, sigproj);
    TPZTensor<REAL> sigPlast(sigmaTrial);
    sigPlast.Add(sigproj, -1.);
    fER.ComputeStrain(sigPlast, Np1.m_eps_p);
    Np1.m_eps_p.Add(N.m_eps_p, 1.);
    validEqs.Fill(0);
    for (int i = 0; i < 2; i++) {
        if (delGamma[i] > 0.) {
            validEqs[i] = 1;
        }
    }
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "epsp next " << Np1.m_hardening << std::endl;
        sout << "delGamma " << delGamma << std::endl;
        sout << "validEqs " << validEqs << std::endl;
        TPZManVector<REAL, 2> Residual(2);
        fYC.Compute(sigmaTrial, N.m_hardening, Residual, 1);
        sout << "residual before projection" << Residual << std::endl;
        fYC.Compute(sigproj, Np1.m_hardening, Residual, 1);
        sout << "residual after projection" << Residual << std::endl;
        LOGPZ_DEBUG(loggerSM, sout.str())
    }
#endif
}

template <>
void TPZPlasticStep<TPZYCSandlerDimaggioL, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>::InitialGuess(
        const TPZPlasticState<REAL> &N,
        TPZPlasticState<REAL> &Np1,
        TPZVec<REAL> &delGamma,
        TPZVec<int> &validEqs
        ) {
    TPZTensor<REAL> EpN = N.m_eps_p;
    TPZTensor<REAL> ETotal = Np1.m_eps_t;
    TPZTensor<REAL> ETrial = ETotal;
    TPZTensor<REAL> sigmaTrial;
    ETrial.Add(EpN, -1.);
    fER.ComputeStress(ETrial, sigmaTrial);
    TPZTensor<REAL> sigproj;
    fYC.InitialGuess(fER, N.m_hardening, sigmaTrial, Np1.m_hardening, delGamma, sigproj);
    TPZTensor<REAL> sigPlast(sigmaTrial);
    sigPlast.Add(sigproj, -1.);
    fER.ComputeStrain(sigPlast, Np1.m_eps_p);
    Np1.m_eps_p.Add(N.m_eps_p, 1.);
    validEqs.Fill(0);
    for (int i = 0; i < 2; i++) {
        if (delGamma[i] > 0.) {
            validEqs[i] = 1;
        }
    }
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "epsp next " << Np1.m_hardening << std::endl;
        sout << "delGamma " << delGamma << std::endl;
        sout << "validEqs " << validEqs << std::endl;
        TPZManVector<REAL, 2> Residual(2);
        fYC.Compute(sigmaTrial, N.m_hardening, Residual, 1);
        sout << "residual before projection" << Residual << std::endl;
        fYC.Compute(sigproj, Np1.m_hardening, Residual, 1);
        sout << "residual after projection" << Residual << std::endl;
        LOGPZ_DEBUG(loggerSM, sout.str())
    }
#endif
}

template <>
void TPZPlasticStep<TPZYCSandlerDimaggioL2, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>::InitialGuess(
        const TPZPlasticState<REAL> &N,
        TPZPlasticState<REAL> &Np1,
        TPZVec<REAL> &delGamma,
        TPZVec<int> &validEqs
        ) {
    TPZTensor<REAL> EpN = N.m_eps_p;
    TPZTensor<REAL> ETotal = Np1.m_eps_t;
    TPZTensor<REAL> ETrial = ETotal;
    TPZTensor<REAL> sigmaTrial;
    ETrial.Add(EpN, -1.);
    fER.ComputeStress(ETrial, sigmaTrial);
    TPZTensor<REAL> sigproj;
    fYC.InitialGuess(fER, N.m_hardening, sigmaTrial, Np1.m_hardening, delGamma, sigproj);
    TPZTensor<REAL> sigPlast(sigmaTrial);
    sigPlast.Add(sigproj, -1.);
    fER.ComputeStrain(sigPlast, Np1.m_eps_p);
    Np1.m_eps_p.Add(N.m_eps_p, 1.);
    validEqs.Fill(0);
    for (int i = 0; i < 2; i++) {
        if (delGamma[i] > 0.) {
            validEqs[i] = 1;
        }
    }
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "epsp next " << Np1.m_hardening << std::endl;
        sout << "delGamma " << delGamma << std::endl;
        sout << "validEqs " << validEqs << std::endl;
        TPZManVector<REAL, 2> Residual(2);
        fYC.Compute(sigmaTrial, N.m_hardening, Residual, 1);
        sout << "residual before projection" << Residual << std::endl;
        fYC.Compute(sigproj, Np1.m_hardening, Residual, 1);
        sout << "residual after projection" << Residual << std::endl;
        LOGPZ_DEBUG(loggerSM, sout.str())
    }
#endif
}



template class TPZPlasticStep<TPZYCSandlerDimaggio, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>;

template class TPZPlasticStep<TPZYCSandlerDimaggioL, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>;

template class TPZPlasticStep<TPZYCSandlerDimaggioL2, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>;


template void TPZPlasticStep<TPZYCSandlerDimaggio, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>::ComputePlasticVars<REAL>(TPZPlasticState<REAL> const&, TPZTensor<REAL>&, REAL&) const;

template void TPZPlasticStep<TPZYCSandlerDimaggioL, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>::ComputePlasticVars<REAL>(TPZPlasticState<REAL> const&, TPZTensor<REAL>&, REAL&) const;

template void TPZPlasticStep<TPZYCSandlerDimaggioL2, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>::ComputePlasticVars<REAL>(TPZPlasticState<REAL> const&, TPZTensor<REAL>&, REAL&) const;
