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
    
private:
    
    /** @brief Temperature @ reservoir conditions  - F */
    REAL fReservoirTemperature;
    
    /** @brief Characteristic Pressure - Pa */
    REAL fPRef;
    
    /** @brief Characteristic Density - kg/m3 */
    REAL fRhoRef;
    
    /** @brief Characteristic viscosity - Pa s */
    REAL fMuRef;
    
    /** @brief Density - kg/m3  $\rho_{g}$ */
    REAL fRho;
    
    /** @brief viscosity - Pa s  $\mu_{g}$ */
    REAL fMu;
    
    /** @brief Compressibility - 1/pa $c_{g}$ */
    REAL fc;
    
    
public:
    
    /** @brief Default constructor $ */
    ReducedPVT();
    
    /** @brief Default desconstructor $ */
    ~ReducedPVT();

    
    /** @brief Density - kg/m3  $\rho$ */
    virtual void Density(TPZVec<REAL> &rho, TPZVec<REAL> state_vars) = 0;
    
    /** @brief viscosity - Pa s  $\mu$ */
    virtual void Viscosity(TPZVec<REAL> &mu, TPZVec<REAL> state_vars) = 0;
    
    /** @brief Compressibility - 1/pa $c$ */
    virtual void Compressibility(TPZVec<REAL> &c, TPZVec<REAL> state_vars) = 0;
        
    /** @brief Set Density - kg/m3  $\rho$ */
    void SetPRef(REAL PRef){this->fPRef = PRef;}
    
    /** @brief Get Density - kg/m3  $\rho$ */
    REAL GetPRef(){return this->fPRef ;}
    
    /** @brief Set Density - kg/m3  $\rho$ */
    void SetRho(REAL Rho){this->fRho = Rho;}
    
    /** @brief Get Density - kg/m3  $\rho$ */
    REAL GetRho(){return this->fRho ;}
    
    /** @brief Set viscosity - pa s  $\mu$ */
    void SetMu(REAL Mu){this->fMu = Mu;}
    
    /** @brief Get viscosity - pa s  $\mu$ */
    REAL GetMu(){return this->fMu ;}
    
    /** @brief Set Compressibility - 1/pa   $C$ */
    void Setc(REAL c){this->fc = c;}
    
    /** @brief Get Compressibility - 1/pa   $C_$ */
    REAL Getc(){return this->fc ;}
    
};


#endif /* defined(__PZ__ReducedPVT__) */