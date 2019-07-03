#ifndef TPZFluidDATAH
#define TPZFluidDATAH
/*
 *  ReservoirData.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include <math.h>

class FluidData {
    
public:
    
    /**
     * @ingroup Characteristic Parameters
     * @brief Define characteristic parameters for Darcy linear flow.
     * @since December 08, 2014
     */
    
    /** @brief Temperature @ reservoir conditions  - F */
    REAL fReservoirTemperature;
    
    /** @brief Characteristic Pressure - Pa */
    REAL fPref;
    
    /** @brief Characteristic Density - kg/m3 */
    REAL fRhoref;
    
    /** @brief Characteristic viscosity - Pa s */
    REAL fMuref;
    
    FluidData();
    
    ~FluidData();
    
    
    /** @brief Set the temperature @ reservoir conditions  - F */
    void SetResT(REAL ResT) {fReservoirTemperature = ResT;}
    
    /** @brief Get the temperature @ reservoir conditions  - F */
    REAL ResT() {return fReservoirTemperature;}
    
    /** @brief Set the characteristic Pressure - Pa */
    void SetPref(REAL Pref) {fPref = Pref;}
    
    /** @brief Characteristic Pressure - Pa */
    REAL Pref() {return fPref;}
    
    /** @brief Set the characteristic Density - kg/m3 */
    void Rhoref(REAL Rhoref) {fRhoref = Rhoref;}
    
    /** @brief Characteristic Density - kg/m3 */
    REAL Rhoref() {return fRhoref;}
    
    /** @brief Set the characteristic viscosity - Pa s */
    void SetMuref(REAL Muref) {fMuref = Muref;}
    
    /** @brief Characteristic viscosity - Pa s */
    REAL Muref() {return fMuref;}
    
    
    // Water properties
    
    /** @brief Water Density - kg/m3  $\rho_{w}$ */
    virtual void WaterDensity(REAL P, REAL &WaterRho, REAL &dWaterRhodP) = 0;
    
    /** @brief Water viscosity - Pa s  $\mu_{w}$ */
    virtual void WaterViscosity(REAL P, REAL &WaterMu, REAL &dWaterMudP) = 0;
    
    /** @brief Water Compressibility - 1/pa $c_{w}$ */
    virtual void WaterCompressibility(REAL P, REAL &cWater, REAL &dcWaterdP) = 0;
    
    
    // Oil properties
    
    /** @brief Water Density - kg/m3  $\rho_{o}$ */
    virtual void OilDensity(REAL P, REAL &OilRho, REAL &dOilRhodP) = 0;
    
    /** @brief Water viscosity - Pa s  $\mu_{o}$ */
    virtual void OilViscosity(REAL P, REAL &OilMu, REAL &dOilMudP) = 0;
    
    /** @brief Water Compressibility - 1/pa $c_{o}$ */
    virtual void OilCompressibility(REAL P, REAL &cOil, REAL &dcOildP) = 0;
    
    
    // Gas properties
    
    /** @brief Water Density - kg/m3  $\rho_{g}$ */
    virtual void GasDensity(REAL P, REAL &GasRho, REAL &dGasRhodP) = 0;
    
    /** @brief Water viscosity - Pa s  $\mu_{g}$ */
    virtual void GasViscosity(REAL P, REAL &GasMu, REAL &dGasMudP) = 0;
    
    /** @brief Water Compressibility - 1/pa $c_{g}$ */
    virtual void GasCompressibility(REAL P, REAL &cGas, REAL &dcGasdP) = 0;
    
    
    
};


#endif