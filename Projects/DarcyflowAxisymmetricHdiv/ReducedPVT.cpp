//
//  ReducedPVT.cpp
//  PZ
//
//  Created by Omar on 4/24/15.
//
//

#include "ReducedPVT.h"


ReducedPVT::ReducedPVT()
{
    
    /** @brief Temperature @ reservoir conditions  - F */
    fReservoirTemperature = 180.0;
    
    /** @brief Characteristic Pressure - Pa */
    fPref = 1.0;
    
    /** @brief Characteristic Density - kg/m3 */
    fRhoref = 1.0;
    
    /** @brief Characteristic viscosity - Pa s */
    fMuref = 1.0;
    
    
    // Water properties
    
    /** @brief Water Density - kg/m3  $\rho_{w}$ */
    fRhoWater = 0.0;
    
    /** @brief Water viscosity - Pa s  $\mu_{w}$ */
    fMuWater = 0.0;
    
    /** @brief Water Compressibility - 1/pa $c_{w}$ */
    fcWater  = 0.0;
    
    
    // Oil properties
    
    /** @brief Water Density - kg/m3  $\rho_{o}$ */
    fRhoOil = 0.0;
    
    /** @brief Water viscosity - Pa s  $\mu_{o}$ */
    fMuOil = 0.0;
    
    /** @brief Water Compressibility - 1/pa $c_{o}$ */
    fcOil = 0.0;
    
    
    // Gas properties
    
    /** @brief Water Density - kg/m3  $\rho_{g}$ */
    fRhoGas = 0.0;
    
    /** @brief Water viscosity - Pa s  $\mu_{g}$ */
    fMuGas = 0.0;
    
    /** @brief Water Compressibility - 1/pa $c_{g}$ */
    fcGas = 0.0;

    
}

ReducedPVT::~ReducedPVT()
{
    
}

// Water properties

/** @brief Water Density - kg/m3  $\rho_{w}$ */
void ReducedPVT::WaterDensity(REAL P, REAL &WaterRho, REAL &dWaterRhodP)
{
    // Linear model
    WaterRho = this->fRhoWater * (1.0 + fcWater * (P - fPref) );
    dWaterRhodP = this->fRhoWater * fcWater;
}

/** @brief Water viscosity - Pa s  $\mu_{w}$ */
void ReducedPVT::WaterViscosity(REAL P, REAL &WaterMu, REAL &dWaterMudP)
{
    WaterMu = fMuWater;
    dWaterMudP = 0.0;
}

/** @brief Water Compressibility - 1/pa $c_{w}$ */
void ReducedPVT::WaterCompressibility(REAL P, REAL &cWater, REAL &dcWaterdP)
{
    cWater = fcWater;
    dcWaterdP = 0.0;
}


// Oil properties

/** @brief Water Density - kg/m3  $\rho_{o}$ */
void ReducedPVT::OilDensity(REAL P, REAL &OilRho, REAL &dOilRhodP)
{
    // Linear model
    OilRho = this->fRhoOil * (1.0 + fcOil * (P - fPref) );
    dOilRhodP = this->fRhoOil * fcOil;
}

/** @brief Water viscosity - Pa s  $\mu_{o}$ */
void ReducedPVT::OilViscosity(REAL P, REAL &OilMu, REAL &dOilMudP)
{
    OilMu = fMuOil;
    dOilMudP = 0.0;
}

/** @brief Water Compressibility - 1/pa $c_{o}$ */
void ReducedPVT::OilCompressibility(REAL P, REAL &cOil, REAL &dcOildP)
{
    cOil = fcOil;
    dcOildP = 0.0;
}


// Gas properties

/** @brief Water Density - kg/m3  $\rho_{g}$ */
void ReducedPVT::GasDensity(REAL P, REAL &GasRho, REAL &dGasRhodP)
{
    GasRho = 0.0;
    dGasRhodP = 0.0;
}

/** @brief Water viscosity - Pa s  $\mu_{g}$ */
void ReducedPVT::GasViscosity(REAL P, REAL &GasMu, REAL &dGasMudP)
{
    GasMu = 0.0;
    dGasMudP = 0.0;
}

/** @brief Water Compressibility - 1/pa $c_{g}$ */
void ReducedPVT::GasCompressibility(REAL P, REAL &cGas, REAL &dcGasdP)
{
    cGas = 0.0;
    dcGasdP = 0.0;
}

