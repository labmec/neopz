//
//  TRMPetrophysicsProperties.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMPetrophysicsProperties.h"


TRMPetrophysicsProperties::TRMPetrophysicsProperties()
{
    
    /** @brief one-phase flow */
    fIsOnePhaseQ = false;
    
    /** @brief two-phase flow */
    fIsTwoPhaseQ = false;
    
    /** @brief three-phase flow */
    fIsThreePhaseQ = false;
    
    /** @brief Definition of the flow system one - two or three phase */
    fSystemType.Resize(0);
    
    /** @brief phase alpha */
    fPhase_alpha = NULL;
    
    /** @brief phase beta */
    fPhase_beta = NULL;
    
    /** @brief phase gamma */
    fPhase_gamma = NULL;
    
    
}

TRMPetrophysicsProperties::~TRMPetrophysicsProperties()
{
    
}

/** @brief Set phase alpha */
void TRMPetrophysicsProperties::SetPhaseAlpha(TPZAutoPointer<TRMPhaseProperties> &alpha)
{
    fPhase_alpha = alpha;
}

/** @brief Get phase alpha */
TPZAutoPointer<TRMPhaseProperties> & TRMPetrophysicsProperties::AlphaProp()
{
    return fPhase_alpha;
}

/** @brief Set phase beta */
void TRMPetrophysicsProperties::SetPhaseBeta(TPZAutoPointer<TRMPhaseProperties> &beta)
{
    fPhase_beta = beta;
}

/** @brief Get phase beta */
TPZAutoPointer<TRMPhaseProperties> & TRMPetrophysicsProperties::BetaProp()
{
    return fPhase_beta;
}

/** @brief Set phase gamma */
void TRMPetrophysicsProperties::SetPhaseGamma(TPZAutoPointer<TRMPhaseProperties> &gamma)
{
    fPhase_gamma = gamma;
}

/** @brief Get phase gamma */
TPZAutoPointer<TRMPhaseProperties> & TRMPetrophysicsProperties::GammaProp()
{
    return fPhase_gamma;
}


// ------------------------------------------------------------------- //
// one phase flow case
// ------------------------------------------------------------------- //

// Capillary Pressure models

/** @brief Alpha-Alpha Capillary Pressure - Pa $P_{caa}$ */
void TRMPetrophysicsProperties::Pc(TPZManVector<STATE,10> &pc, TPZManVector<STATE,10> &x){
    
    int n = x.size() + 1;
    STATE val = 0.0;
    pc.Resize(n,0.0);
    pc[0] = val;
    
}

// Relative permeabilities models

/** @brief Alpha-Alpha phase realtive permeability $k_{ra}$ */
void TRMPetrophysicsProperties::Kr(TPZManVector<STATE,10> &kr, TPZManVector<STATE,10> &x){
    
    int n = x.size() + 1;
    STATE val = 1.0;
    kr.Resize(n,0.0);
    kr[0] = val;
    
}

// Mobilities models

/** @brief Alpha-Alpha phase moblity $l_{a}$ */
void TRMPetrophysicsProperties::l_single(TPZManVector<STATE,10> &l, TPZManVector<STATE,10> &x){
    
    int n = x.size() + 1;
    l.Resize(n,0.0);
    
    TPZManVector<STATE,10> rho,mu,kr;
    this->fPhase_alpha->Density(rho, x);
    this->fPhase_alpha->Viscosity(mu, x);
    this->Kr(kr, x);
    
    l[0] = kr[0]*rho[0]/mu[0];
    l[1] = kr[0]*((rho[1]/mu[0])-(rho[0]*mu[1]/(mu[0]*mu[0])));
    l[2] = kr[2]*rho[0]/mu[0];
    l[3] = kr[3]*rho[0]/mu[0];
    l[4] = kr[0]*((rho[4]/mu[0])-(rho[0]*mu[4]/(mu[0]*mu[0])));
    
}

// ------------------------------------------------------------------- //
// two phase flow case
// ------------------------------------------------------------------- //

// Capillary Pressure models

/** @brief Alpha-Alpha Capillary Pressure - Pa $P_{caa}$ */
void TRMPetrophysicsProperties::Pcab(TPZManVector<STATE,10> &pc, TPZManVector<STATE,10> &x){

    int n = x.size() + 1;
    STATE pmax = 0.0;
    pc.Resize(n,0.0);
    pc[0] = pmax*(1-x[1]);
    pc[2] = -pmax;
    
}

// Relative permeabilities models

/** @brief Alpha-Alpha phase realtive permeability $k_{ra}$ */
void TRMPetrophysicsProperties::Kra(TPZManVector<STATE,10> &kr, TPZManVector<STATE,10> &x){
    
    int n = x.size() + 1;
    kr.Resize(n,0.0);
    kr[0] = x[1];
    kr[2] = 1.0;
    
}

/** @brief Beta-Beta phase realtive permeability $k_{ra}$ */
void TRMPetrophysicsProperties::Krb(TPZManVector<STATE,10> &kr, TPZManVector<STATE,10> &x){

    int n = x.size() + 1;
    kr.Resize(n,0.0);
    kr[0] = (1-x[1]);
    kr[2] = -1.0;
}

// Mobilities models

/** @brief Alpha-Alpha phase moblity $l_{a}$ */
void TRMPetrophysicsProperties::la(TPZManVector<STATE,10> &l, TPZManVector<STATE,10> &x){
   
    int n = x.size() + 1;
    l.Resize(n,0.0);
    
    TPZManVector<STATE,10> rho_a,mu_a,kr_a;
    this->fPhase_alpha->Density(rho_a, x);
    this->fPhase_alpha->Viscosity(mu_a, x);
    this->Kra(kr_a, x);
    
    l[0] = kr_a[0]*rho_a[0]/mu_a[0];
    l[1] = kr_a[0]*((rho_a[1]/mu_a[0])-(rho_a[0]*mu_a[1]/(mu_a[0]*mu_a[0])));
    l[2] = kr_a[2]*rho_a[0]/mu_a[0];
    l[3] = kr_a[3]*rho_a[0]/mu_a[0];
    l[4] = kr_a[0]*((rho_a[4]/mu_a[0])-(rho_a[0]*mu_a[4]/(mu_a[0]*mu_a[0])));
    
}

/** @brief Beta-Beta phase moblity $l_{a}$ */
void TRMPetrophysicsProperties::lb(TPZManVector<STATE,10> &l, TPZManVector<STATE,10> &x){
    
    int n = x.size() + 1;
    l.Resize(n,0.0);
    
    TPZManVector<STATE,10> rho_b,mu_b,kr_b;
    this->fPhase_beta->Density(rho_b, x);
    this->fPhase_beta->Viscosity(mu_b, x);
    this->Krb(kr_b, x);
    
    l[0] = kr_b[0]*rho_b[0]/mu_b[0];
    l[1] = kr_b[0]*((rho_b[1]/mu_b[0])-(rho_b[0]*mu_b[1]/(mu_b[0]*mu_b[0])));
    l[2] = kr_b[2]*rho_b[0]/mu_b[0];
    l[3] = kr_b[3]*rho_b[0]/mu_b[0];
    l[4] = kr_b[0]*((rho_b[4]/mu_b[0])-(rho_b[0]*mu_b[4]/(mu_b[0]*mu_b[0])));
    
}

/** @brief Alpha-Alpha phase fractional flow $f_{a}$ */
void TRMPetrophysicsProperties::fa(TPZManVector<STATE,10> &f, TPZManVector<STATE,10> &x){
    
    int n = x.size() + 1;
    f.Resize(n,0.0);
    
    TPZManVector<STATE,10> l,la;
    this->la(la, x);
    this->l(l, x);
    
    f[0] = la[0]/l[0];
    f[1] = (la[1]/l[0]) - (la[0]*l[1]/(l[0]*l[0]));
    f[2] = (la[2]/l[0]) - (la[0]*l[2]/(l[0]*l[0]));
    f[3] = (la[3]/l[0]) - (la[0]*l[3]/(l[0]*l[0]));
    f[4] = (la[4]/l[0]) - (la[0]*l[4]/(l[0]*l[0]));

    
}

/** @brief Beta-Beta phase fractional flow $f_{a}$ */
void TRMPetrophysicsProperties::fb(TPZManVector<STATE,10> &f, TPZManVector<STATE,10> &x){
    
    int n = x.size() + 1;
    f.Resize(n,0.0);
    
    TPZManVector<STATE,10> l,lb;
    this->lb(lb, x);
    this->l(l, x);
    
    f[0] = lb[0]/l[0];
    f[1] = (lb[1]/l[0]) - (lb[0]*l[1]/(l[0]*l[0]));
    f[2] = (lb[2]/l[0]) - (lb[0]*l[2]/(l[0]*l[0]));
    f[3] = (lb[3]/l[0]) - (lb[0]*l[3]/(l[0]*l[0]));
    f[4] = (lb[4]/l[0]) - (lb[0]*l[4]/(l[0]*l[0]));
    
}


// ------------------------------------------------------------------- //
// three phase flow case
// ------------------------------------------------------------------- //


/** @brief Total moblity $l_{a}$ */
void TRMPetrophysicsProperties::l(TPZManVector<STATE,10> &l, TPZManVector<STATE,10> &x){
    
#ifdef PZDEBUG
    if(SystemType().size() == 0){
        DebugStop();
    }
#endif
    
    switch (fSystemType.size()) {
        case 1:
        {
            this->l_single(l, x);
        }
            break;
        case 2:
        {
            TPZManVector<STATE,10> l_a,l_b;
            this->la(l_a, x);
            this->lb(l_b, x);
            l.Resize(l_a.size(), 0.0);
            l[0] = l_a[0] + l_b[0];
            l[1] = l_a[1] + l_b[1];
            l[2] = l_a[2] + l_b[2];
            l[3] = l_a[3] + l_b[3];
            l[4] = l_a[4] + l_b[4];
            
        }
            break;
        case 3:
        {
            DebugStop();
        }
            break;
            
        default:
        {
            std::cout << "This code run just three-phasic systems" << std::endl;
            DebugStop();
        }
            break;
    }
    
}
