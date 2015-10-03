//
//  OilPhase.h
//  PZ
//
//  Created by Omar on 9/30/15.
//
//

#ifndef __PZ__OilPhase__
#define __PZ__OilPhase__

#include <stdio.h>
#include "ReducedPVT.h"

#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include <math.h>


class OilPhase : public ReducedPVT

{
    
private:
    
//    /** @brief Temperature @ reservoir conditions  - F */
//    REAL fReservoirTemperature;
//    
//    /** @brief Characteristic Pressure - Pa */
//    REAL fPRef;
//    
//    /** @brief Characteristic Density - kg/m3 */
//    REAL fRhoRef;
//    
//    /** @brief Characteristic viscosity - Pa s */
//    REAL fMuRef;
//    
//    /** @brief Density - kg/m3  $\rho_{g}$ */
//    REAL fRho;
//    
//    /** @brief viscosity - Pa s  $\mu_{g}$ */
//    REAL fMu;
//    
//    /** @brief Compressibility - 1/pa $c_{g}$ */
//    REAL fc;
    
    
public:
    
    /** @brief Default constructor $ */
    OilPhase();
    
    /** @brief Default desconstructor $ */
    ~OilPhase();
    
//    /** @brief Set the temperature @ reservoir conditions  - F */
//    void SetResT(REAL ResT) {fReservoirTemperature = ResT;}
//    
//    /** @brief Get the temperature @ reservoir conditions  - F */
//    REAL ResT() {return fReservoirTemperature;}
//    
//    /** @brief Set the characteristic Pressure - Pa */
//    void SetPRef(REAL Pref) {fPRef = Pref;}
//    
//    /** @brief Characteristic Pressure - Pa */
//    REAL PRef() {return fPRef;}
//    
//    /** @brief Set the characteristic Density - kg/m3 */
//    void RhoRef(REAL Rhoref) {fRhoRef = Rhoref;}
//    
//    /** @brief Characteristic Density - kg/m3 */
//    REAL RhoRef() {return fRhoRef;}
//    
//    /** @brief Set the characteristic viscosity - Pa s */
//    void SetMuRef(REAL Muref) {fMuRef = Muref;}
//    
//    /** @brief Characteristic viscosity - Pa s */
//    REAL MuRef() {return fMuRef;}
    
    /** @brief Density - kg/m3  $\rho$ */
    void Density(TPZVec<REAL> &rho, TPZVec<REAL> state_vars);
    
    /** @brief viscosity - Pa s  $\mu$ */
    void Viscosity(TPZVec<REAL> &mu, TPZVec<REAL> state_vars);
    
    /** @brief Compressibility - 1/pa $c$ */
    void Compressibility(TPZVec<REAL> &c, TPZVec<REAL> state_vars);
    
//    /** @brief Set Density - kg/m3  $\rho$ */
//    void SetRho(REAL Rho){this->fRho = Rho;}
//    
//    /** @brief Get Density - kg/m3  $\rho$ */
//    REAL GetRho(){return this->fRho ;}
//    
//    /** @brief Set viscosity - pa s  $\mu$ */
//    void SetMu(REAL Mu){this->fMu = Mu;}
//    
//    /** @brief Get viscosity - pa s  $\mu$ */
//    REAL GetMu(){return this->fMu ;}
//    
//    /** @brief Set Compressibility - 1/pa   $C$ */
//    void Setc(REAL c){this->fc = c;}
//    
//    /** @brief Get Compressibility - 1/pa   $C_$ */
//    REAL Getc(){return this->fc ;}
    
};


#endif /* defined(__PZ__OilPhase__) */
