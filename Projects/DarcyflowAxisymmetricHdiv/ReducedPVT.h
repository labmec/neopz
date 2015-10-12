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
    REAL fTRes;
    
    /** @brief Temperature for references values [K] */
    REAL fTRef;
    
    /** @brief Pressure for references values - Pa */
    REAL fPRef;
    
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
    
    /** @brief Copy constructor $ */
    ReducedPVT(const ReducedPVT& other)
    {
        fTRes       = other.fTRes;
        fPRef       = other.fPRef;
        fRho        = other.fRho;
        fMu         = other.fMu;
        fc          = other.fc;
        
    }
    
    /** @brief Copy assignemnt operator $ */
    ReducedPVT& operator = (const ReducedPVT& other)
    {
        if (this != & other) // prevent self-assignment
        {
            
            fTRes       = other.fTRes;
            fPRef       = other.fPRef;
            fRho        = other.fRho;
            fMu         = other.fMu;
            fc          = other.fc;
            
        }
        return *this;
    }

    /** @brief Density - kg/m3  $\rho$ */
    virtual void Density(TPZManVector<REAL> &rho, TPZManVector<REAL> state_vars) = 0;
    
    /** @brief viscosity - Pa s  $\mu$ */
    virtual void Viscosity(TPZManVector<REAL> &mu, TPZManVector<REAL> state_vars) = 0;
    
    /** @brief Compressibility - 1/pa $c$ */
    virtual void Compressibility(TPZManVector<REAL> &c, TPZManVector<REAL> state_vars) = 0;
    
    /** @brief Set Reservoir T - K  $T_{res}$ */
    void SetTRes(REAL TRes){this->fTRes = TRes;}
    
    /** @brief Get Reservoir T - K  $T_{res}$ */
    REAL GetTRes(){return this->fTRes ;}
    
    /** @brief Set Reference Temperature - K  $T$ */
    void SetTRef(REAL TRef){this->fTRef = TRef;}
    
    /** @brief Get Reference Temperature - K  $T$ */
    REAL GetTRef(){return this->fTRef ;}
    
    /** @brief Set Reference Pressure - Pa  $P_{ref}$ */
    void SetPRef(REAL PRef){this->fPRef = PRef;}
    
    /** @brief Get Reference Pressure - Pa  $P_{ref}$ */
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