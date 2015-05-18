//
//  ReducedPVT.h
//  PZ
//
//  Created by Omar on 4/24/15.
//
//

#ifndef __PZ__ReducedPVT__
#define __PZ__ReducedPVT__

#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include <math.h>


class ReducedPVT

{
    
public:
    
    
    /** @brief Temperature @ reservoir conditions  - F */
    REAL fReservoirTemperature;
    
    /** @brief Characteristic Pressure - Pa */
    REAL fPref;
    
    /** @brief Characteristic Density - kg/m3 */
    REAL fRhoref;
    
    /** @brief Characteristic viscosity - Pa s */
    REAL fMuref;
    
    
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
    REAL fRhoWater;

    /** @brief Water viscosity - Pa s  $\mu_{w}$ */
    REAL fMuWater;

    /** @brief Water Compressibility - 1/pa $c_{w}$ */
    REAL fcWater;


    // Oil properties

    /** @brief Water Density - kg/m3  $\rho_{o}$ */
    REAL fRhoOil;

    /** @brief Water viscosity - Pa s  $\mu_{o}$ */
    REAL fMuOil;

    /** @brief Water Compressibility - 1/pa $c_{o}$ */
    REAL fcOil;


    // Gas properties

    /** @brief Water Density - kg/m3  $\rho_{g}$ */
    REAL fRhoGas;

    /** @brief Water viscosity - Pa s  $\mu_{g}$ */
    REAL fMuGas;

    /** @brief Water Compressibility - 1/pa $c_{g}$ */
    REAL fcGas;
    
    ReducedPVT();
    
    ~ReducedPVT();
    
    
    // Water properties
    
    /** @brief Water Density - kg/m3  $\rho_{w}$ */
    void WaterDensity(REAL P, REAL &WaterRho, REAL &dWaterRhodP);
    
    /** @brief Water viscosity - Pa s  $\mu_{w}$ */
    void WaterViscosity(REAL P, REAL &WaterMu, REAL &dWaterMudP);
    
    /** @brief Water Compressibility - 1/pa $c_{w}$ */
    void WaterCompressibility(REAL P, REAL &cWater, REAL &dcWaterdP);
    
    
    // Oil properties
    
    /** @brief Water Density - kg/m3  $\rho_{o}$ */
    void OilDensity(REAL P, REAL &OilRho, REAL &dOilRhodP);
    
    /** @brief Water viscosity - Pa s  $\mu_{o}$ */
    void OilViscosity(REAL P, REAL &OilMu, REAL &dOilMudP);
    
    /** @brief Water Compressibility - 1/pa $c_{o}$ */
    void OilCompressibility(REAL P, REAL &cOil, REAL &dcOildP);
    
    
    // Gas properties
    
    /** @brief Water Density - kg/m3  $\rho_{g}$ */
    void GasDensity(REAL P, REAL &GasRho, REAL &dGasRhodP);
    
    /** @brief Water viscosity - Pa s  $\mu_{g}$ */
    void GasViscosity(REAL P, REAL &GasMu, REAL &dGasMudP);
    
    /** @brief Water Compressibility - 1/pa $c_{g}$ */
    void GasCompressibility(REAL P, REAL &cGas, REAL &dcGasdP);
    
    
    // Water properties
    
    /** @brief Set Water Density - kg/m3  $\rho_{w}$ */
    void SetRhoWater(REAL RhoWater){this->fRhoWater = RhoWater;}
    
    /** @brief Get Water Density - kg/m3  $\rho_{w}$ */
    REAL GetRhoWater(REAL RhoWater){return this->fRhoWater ;}
    
    /** @brief Set Water viscosity - pa s  $\mu_{w}$ */
    void SetMuWater(REAL MuWater){this->fMuWater = MuWater;}
    
    /** @brief Get Water viscosity - pa s  $\mu_{w}$ */
    REAL GetMuWater(REAL MuWater){return this->fMuWater ;}
    
    /** @brief Set Water Compressibility - 1/pa   $C_{w}$ */
    void SetcWater(REAL cWater){this->fcWater = cWater;}
    
    /** @brief Get Water Compressibility - 1/pa   $C_{w}$ */
    REAL GetcWater(REAL cWater){return this->fcWater ;}
    
    
    
    // Oil properties

    
    /** @brief Set Oil Density - kg/m3  $\rho_{o}$ */
    void SetRhoOil(REAL RhoOil){this->fRhoOil = RhoOil;}
    
    /** @brief Get Oil Density - kg/m3  $\rho_{o}$ */
    REAL GetRhoOil(REAL RhoOil){return this->fRhoOil ;}
    
    /** @brief Set Oil viscosity - pa s  $\mu_{o}$ */
    void SetMuOil(REAL MuOil){this->fMuOil = MuOil;}
    
    /** @brief Get Oil viscosity - pa s  $\mu_{o}$ */
    REAL GetMuOil(REAL MuOil){return this->fMuOil ;}
    
    /** @brief Set Oil Compressibility - 1/pa   $C_{o}$ */
    void SetcOil(REAL cOil){this->fcOil = cOil;}
    
    /** @brief Get Oil Compressibility - 1/pa   $C_{o}$ */
    REAL GetcOil(REAL cOil){return this->fcOil ;}
    
    
    
    // Gas properties

    
    /** @brief Set Gas Density - kg/m3  $\rho_{g}$ */
    void SetRhoGas(REAL RhoGas){this->fRhoGas = RhoGas;}
    
    /** @brief Get Gas Density - kg/m3  $\rho_{g}$ */
    REAL GetRhoGas(REAL RhoGas){return this->fRhoGas ;}
    
    /** @brief Set Gas viscosity - pa s  $\mu_{g}$ */
    void SetMuGas(REAL MuGas){this->fMuGas = MuGas;}
    
    /** @brief Get Gas viscosity - pa s  $\mu_{g}$ */
    REAL GetMuGas(REAL MuGas){return this->fMuGas ;}
    
    /** @brief Set Gas Compressibility - 1/pa   $C_{g}$ */
    void SetcGas(REAL cGas){this->fcGas = cGas;}
    
    /** @brief Get Gas Compressibility - 1/pa   $C_{g}$ */
    REAL GetcGas(REAL cGas){return this->fcGas ;}

    
    
    
};


#endif /* defined(__PZ__ReducedPVT__) */