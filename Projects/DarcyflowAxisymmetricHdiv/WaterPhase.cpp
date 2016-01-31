//
//  WaterPhase.cpp
//  PZ
//
//  Created by Omar on 9/30/15.
//
//

#include "WaterPhase.h"


WaterPhase::WaterPhase() : Phase()
{
    
}

WaterPhase::~WaterPhase()
{
    
}


/** @brief Density - kg/m3  $\rho$ */
void WaterPhase::Density(TPZManVector<REAL> &rho, TPZManVector<REAL> state_vars)
{
    REAL Pw = state_vars[1];
    
    rho[0] = GetRho() * exp(  Getc() * (Pw - GetPRef() ));
    rho[2] = GetRho() * Getc() * exp(  Getc() * (Pw - GetPRef() ));
}

/** @brief viscosity - Pa s  $\mu$ */
void WaterPhase::Viscosity(TPZManVector<REAL> &mu, TPZManVector<REAL> state_vars)
{
    mu[0] = GetMu();
    mu[2] = 0.0;
}

/** @brief Compressibility - 1/pa $c$ */
void WaterPhase::Compressibility(TPZManVector<REAL> &c, TPZManVector<REAL> state_vars)
{
    c[0] = Getc();
    c[2] = 0.0;
}

/** @brief Kr - $k_{r}$ */
void WaterPhase::Kr(TPZManVector<REAL> &kr, TPZManVector<REAL> state_vars){
   
    REAL Sw = state_vars[2];
    REAL Swr = GetS_wett_r();
    REAL Sor = GetS_nwett_r();
    REAL Se = (Sw-Swr)/(1.0-Swr-Sor);

   
    if (fIsNonlinearKrQ) {
        kr[0] = Se*Se;
        kr[1] = 0.0;
        kr[2] = 0.0;
        kr[3] = 2.0*Se*(1.0)/(1.0-Swr-Sor);
        kr[4] = 0.0;
    }
    else{
        
        kr[0] = Se;
        kr[1] = 0.0;
        kr[2] = 0.0;
        kr[3] = 1.0*(1.0)/(1.0-Swr-Sor);
        kr[4] = 0.0;
        
    }
    
}

/** @brief Pc - $P_{c}$ */
void WaterPhase::Pc(TPZManVector<REAL> &pc, TPZManVector<REAL> state_vars){
    
//    REAL So = state_vars[2];
    
    pc[0] = 0.0;
    pc[1] = 0.0;
    pc[2] = 0.0;
    pc[3] = 0.0;
    pc[4] = 0.0;
    
}