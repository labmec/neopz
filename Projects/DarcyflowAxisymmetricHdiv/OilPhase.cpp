//
//  OilPhase.cpp
//  PZ
//
//  Created by Omar on 9/30/15.
//
//

#include "OilPhase.h"


OilPhase::OilPhase() : ReducedPVT()
{
        
}

OilPhase::~OilPhase()
{
    
}


/** @brief Density - kg/m3  $\rho$ */
void OilPhase::Density(TPZManVector<REAL> &rho, TPZManVector<REAL> state_vars)
{
    REAL Po = state_vars[1];
    
    rho[0] = GetRho() * ( 1.0 + Getc() * (Po - GetPRef() ));
    rho[2] = GetRho() * Getc();
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