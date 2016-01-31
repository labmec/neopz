//
//  OilPhase.cpp
//  PZ
//
//  Created by Omar on 9/30/15.
//
//

#include "OilPhase.h"


OilPhase::OilPhase() : Phase()
{
        
}

OilPhase::~OilPhase()
{
    
}


/** @brief Density - kg/m3  $\rho$ */
void OilPhase::Density(TPZManVector<REAL> &rho, TPZManVector<REAL> state_vars)
{
    REAL Po = state_vars[1];
    
    rho[0] = GetRho() * exp(  Getc() * (Po - GetPRef()) );
    rho[2] = GetRho() * Getc() * exp(  Getc() * (Po - GetPRef()) );
}

/** @brief viscosity - Pa s  $\mu$ */
void OilPhase::Viscosity(TPZManVector<REAL> &mu, TPZManVector<REAL> state_vars)
{
    mu[0] = GetMu();
    mu[2] = 0.0;
}

/** @brief Compressibility - 1/pa $c$ */
void OilPhase::Compressibility(TPZManVector<REAL> &c, TPZManVector<REAL> state_vars)
{
    c[0] = Getc();
    c[2] = 0.0;
}

/** @brief Kr - $k_{r}$ */
void OilPhase::Kr(TPZManVector<REAL> &kr, TPZManVector<REAL> state_vars){

    REAL So = state_vars[2];
    REAL Swr = GetS_wett_r();
    REAL Sor = GetS_nwett_r();
    REAL Se =  1.0 - (1.0 - So - Swr)/(1.0-Swr-Sor);

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
void OilPhase::Pc(TPZManVector<REAL> &pc, TPZManVector<REAL> state_vars){
    
    REAL S_alpha = state_vars[2];
    REAL Pc_max  = GetPc_max();
    
    pc[0] = (1.0-S_alpha)*Pc_max;
    pc[1] = 0.0;
    pc[2] = 0.0;
    pc[3] = (-1.0)*Pc_max;
    pc[4] = 0.0;
    
}
