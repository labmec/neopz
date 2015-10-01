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
    
    /** @brief Temperature @ reservoir conditions  - F */
    fReservoirTemperature = 180.0;
    
    /** @brief Characteristic Pressure - Pa */
    fPRef = 1.0;
    
    /** @brief Characteristic Density - kg/m3 */
    fRhoRef = 1.0;
    
    /** @brief Characteristic viscosity - Pa s */
    fMuRef = 1.0;
    
    /** @brief Density - kg/m3  $\rho$ */
    fRho = 0.0;
    
    /** @brief viscosity - Pa s  $\mu$ */
    fMu = 0.0;
    
    /** @brief Compressibility - 1/pa $c$ */
    fc = 0.0;
    
}

WaterPhase::~WaterPhase()
{
    
}


/** @brief Density - kg/m3  $\rho$ */
void WaterPhase::Density(TPZVec<REAL> &rho, TPZVec<REAL> state_vars)
{
    rho[0] = 0.0;
    rho[2] = 0.0;
}

/** @brief viscosity - Pa s  $\mu$ */
void WaterPhase::Viscosity(TPZVec<REAL> &mu, TPZVec<REAL> state_vars)
{
    mu[0] = 0.0;
    mu[2] = 0.0;
}

/** @brief Compressibility - 1/pa $c$ */
void WaterPhase::Compressibility(TPZVec<REAL> &c, TPZVec<REAL> state_vars)
{
    c[0] = 0.0;
    c[2] = 0.0;
}