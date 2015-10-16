//
//  WaterPhase.cpp
//  PZ
//
//  Created by Omar on 9/30/15.
//
//

#include "WaterPhase.h"


WaterPhase::WaterPhase() : ReducedPVT()
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
