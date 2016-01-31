//
//  Phase.h
//  PZ
//
//  Created by Omar on 4/24/15.
//
//

#ifndef __PZ__Phase__
#define __PZ__Phase__

#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include <math.h>


class Phase

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
    
    /** @brief Irreducible Saturation of the wetting phase - */
    REAL fS_wett_r;
    
    /** @brief Irreducible Saturation of the no wetting phase - */
    REAL fS_nwett_r;
    
    /** @brief Maximum capillary pressure value - */
    REAL fPc_max;
    
protected:
    
    /** @brief Relative permeability model */
    bool fIsNonlinearKrQ;
    
    
public:
    
    /** @brief Default constructor $ */
    Phase();
    
    /** @brief Default desconstructor $ */
    ~Phase();
    
    /** @brief Copy constructor $ */
    Phase(const Phase& other)
    {
        fTRes       = other.fTRes;
        fPRef       = other.fPRef;
        fRho        = other.fRho;
        fMu         = other.fMu;
        fc          = other.fc;
        
    }
    
    /** @brief Copy assignemnt operator $ */
    Phase& operator = (const Phase& other)
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
    
    /** @brief Relative permeability - $k_{r}$ */
    virtual void Kr(TPZManVector<REAL> &kr, TPZManVector<REAL> state_vars) = 0;
    
    /** @brief Capillar pressure - $k_{r}$ */
    virtual void Pc(TPZManVector<REAL> &pc, TPZManVector<REAL> state_vars) = 0;
    
    /** @brief Set Reservoir Temperature - Kelvin  $T_{res}$ */
    void SetTRes(REAL TRes){this->fTRes = TRes;}
    
    /** @brief Get Reservoir Temperature - Kelvin  $T_{res}$ */
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
    
    /** @brief Set the irreducible Saturation of the wetting phase - */
    void SetS_wett_r(REAL S_w_r) {fS_wett_r = S_w_r;}
    
    /** @brief Get the irreducible Saturation of the wetting phase - */
    REAL GetS_wett_r() {return fS_wett_r;}
    
    /** @brief Set the irreducible Saturation of the no wetting phase - */
    void SetS_nwett_r(REAL S_nw_r) {fS_nwett_r = S_nw_r;}
    
    /** @brief Get the irreducible Saturation of the no wetting phase - */
    REAL GetS_nwett_r() {return fS_nwett_r;}
    
    /** @brief Set maximum capillary pressure value - Pa  $P_{ref}$ */
    void SetPc_max(REAL Pc){this->fPc_max = Pc;}
    
    /** @brief Get maximum capillary pressure value - Pa  $P_{ref}$ */
    REAL GetPc_max(){return this->fPc_max ;}
    
    /** @brief Set the relative permeability model */
    void SetNonlinearKr(bool isnonlinearQ) {fIsNonlinearKrQ = isnonlinearQ; }
    
    /** @brief Get the relative permeability model */
    bool IsNonlinearKrQ() {return fIsNonlinearKrQ; }
    
    
};


#endif /* defined(__PZ__Phase__) */