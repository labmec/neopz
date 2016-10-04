//
//  TRMGasPhase.cpp
//  PZ
//
//  Created by Omar on 5/23/16.
//
//

#include "TRMGasPhase.h"


/** @brief default constructor */
TRMGasPhase::TRMGasPhase() : TRMPhaseProperties() {
    
}

/** @brief default destructor */
TRMGasPhase::~TRMGasPhase(){
    
}

/** @brief Density - kg/m3  $\rho$ */
void TRMGasPhase::Density(TPZManVector<STATE,10> &rho, TPZManVector<STATE,10> &state_vars){
    
    switch (this->RhoModel()) {
        case 0:
        {
            this->Density_c(rho, state_vars);
        }
            break;
        case 1:
        {
            this->Density_l(rho, state_vars);
        }
            break;
        case 2:
        {
            DebugStop();
        }
            break;
        default:
        {
            std::cout << "Error: Model not implemented." << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Viscosity - Pa s  $\mu$ */
void TRMGasPhase::Viscosity(TPZManVector<STATE,10> &mu, TPZManVector<STATE,10> &state_vars){
    
    switch (this->MuModel()) {
        case 0:
        {
            this->Viscosity_c(mu, state_vars);
        }
            break;
        case 1:
        {
            this->Viscosity_l(mu, state_vars);
        }
            break;
        case 2:
        {
            DebugStop();
        }
            break;
        default:
        {
            std::cout << "Error: Model not implemented." << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Compressibility - 1/pa $c$ */
void TRMGasPhase::Compressibility(TPZManVector<STATE,10> &c, TPZManVector<STATE,10> &state_vars){
    
    switch (this->CModel()) {
        case 0:
        {
            this->Compressibility_c(c, state_vars);
        }
            break;
        case 1:
        {
            this->Compressibility_l(c, state_vars);
        }
            break;
        case 2:
        {
            DebugStop();
        }
            break;
        default:
        {
            std::cout << "Error: Model not implemented." << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/**
 * @defgroup Constant models
 * @{
 */

/** @brief Density - kg/m3  $\rho$ */
void TRMGasPhase::Density_c(TPZManVector<STATE,10> &rho, TPZManVector<STATE,10> &state_vars){
    
#ifdef PZDEBUG
    if (state_vars.size() == 0) {
        DebugStop();
    }
#endif
    
    int n = state_vars.size() + 1;
    STATE val = 10.0;
    rho.Resize(n,0.0);
    rho[0] = val;
    
}


/** @brief Viscosity - Pa s  $\mu$ */
void TRMGasPhase::Viscosity_c(TPZManVector<STATE,10> &mu, TPZManVector<STATE,10> &state_vars)
{
    
#ifdef PZDEBUG
    if (state_vars.size() == 0) {
        DebugStop();
    }
#endif
    
    int n = state_vars.size() + 1;
    STATE val = 1.0e-5;
    mu.Resize(n,0.0);
    mu[0] = val;
    
}

/** @brief Compressibility - 1/pa $c$ */
void TRMGasPhase::Compressibility_c(TPZManVector<STATE,10> &c, TPZManVector<STATE,10> &state_vars){
    
#ifdef PZDEBUG
    if (state_vars.size() == 0) {
        DebugStop();
    }
#endif
    
    int n = state_vars.size() + 1;
    STATE val = 1.0e-7;
    c.Resize(n,0.0);
    c[0] = val;
    
}

// @}

/**
 * @defgroup Linearized models
 * @{
 */

/** @brief Density - kg/m3  $\rho$ */
void TRMGasPhase::Density_l(TPZManVector<STATE,10> &rho, TPZManVector<STATE,10> &state_vars){
    
#ifdef PZDEBUG
    if (state_vars.size() == 0) {
        DebugStop();
    }
#endif
    
    int n = state_vars.size() + 1;
    TPZManVector<STATE,10> c(n,0.0);
    
    this->Density_c(rho, state_vars);
    this->Compressibility_c(c, state_vars);
    
    STATE val = rho[0];
    STATE p = state_vars[0];
    
    rho[0] = val*(1.0+c[0]*(p - fPRef));
    rho[1] = val*c[0];
    
}

/** @brief Viscosity - Pa s  $\mu$ */
void TRMGasPhase::Viscosity_l(TPZManVector<STATE,10> &mu, TPZManVector<STATE,10> &state_vars){
    
#ifdef PZDEBUG
    if (state_vars.size() == 0) {
        DebugStop();
    }
#endif
    
    this->Viscosity_c(mu, state_vars);
    
    STATE val = mu[0];
    STATE p = state_vars[0];
    STATE c = 1.0e-20;
    
    mu[0] = val*(1.0+c*(p - fPRef));
    mu[1] = val*c;
    
}

/** @brief Compressibility - 1/pa $c$ */
void TRMGasPhase::Compressibility_l(TPZManVector<STATE,10> &c, TPZManVector<STATE,10> &state_vars){
    
#ifdef PZDEBUG
    if (state_vars.size() == 0) {
        DebugStop();
    }
#endif
    
    this->Compressibility_c(c, state_vars);
    
    STATE val = c[0];
    STATE p = state_vars[0];
    STATE s = 1.0e-15;
    
    c[0] = val*(1.0+s*(p - fPRef));
    c[1] = val*s;
    
}