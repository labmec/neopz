//
//  GasPhase.cpp
//  PZ
//
//  Created by Omar on 9/30/15.
//
//

#include "GasPhase.h"


GasPhase::GasPhase() : ReducedPVT()
{

}

GasPhase::~GasPhase()
{
    
}


/** @brief Density - kg/m3  $\rho$ */
void GasPhase::Density(TPZManVector<REAL> &rho, TPZManVector<REAL> state_vars)
{
    rho[0] = 0.0;
    rho[2] = 0.0;
}

/** @brief viscosity - Pa s  $\mu$ */
void GasPhase::Viscosity(TPZManVector<REAL> &mu, TPZManVector<REAL> state_vars)
{
    mu[0] = 0.0;
    mu[2] = 0.0;
}

/** @brief Compressibility - 1/pa $c$ */
void GasPhase::Compressibility(TPZManVector<REAL> &c, TPZManVector<REAL> state_vars)
{
    c[0] = 0.0;
    c[2] = 0.0;
}