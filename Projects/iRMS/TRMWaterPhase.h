//
//  TRMWaterPhase.h
//  PZ
//
//  Created by Omar on 5/23/16.
//
//

#ifndef TRMWaterPhase_h
#define TRMWaterPhase_h

#include <stdio.h>
#include "TRMPhaseProperties.h"

class TRMWaterPhase : public TRMPhaseProperties {
    
public:
    
    /** @brief default constructor */
    TRMWaterPhase();
    
    /** @brief default constructor */
    TRMWaterPhase(const TRMWaterPhase &copy)
    {
        DebugStop();
    }
    
    /** @brief default constructor */
    TRMWaterPhase &operator=(const TRMWaterPhase &copy)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief default destructor */
    ~TRMWaterPhase();
    
    /** @brief Density - kg/m3  $\rho$ */
    void Density(TPZManVector<STATE,10> &rho, TPZManVector<STATE,10> &state_vars);
    
    /** @brief Viscosity - Pa s  $\mu$ */
    void Viscosity(TPZManVector<STATE,10> &mu, TPZManVector<STATE,10> &state_vars);
    
    /** @brief Compressibility - 1/pa $c$ */
    void Compressibility(TPZManVector<STATE,10> &c, TPZManVector<STATE,10> &state_vars);
    
    /**
     * @defgroup Constant models
     * @{
     */
    
    /** @brief Density - kg/m3  $\rho$ */
    void Density_c(TPZManVector<STATE,10> &rho, TPZManVector<STATE,10> &state_vars);
    
    /** @brief Viscosity - Pa s  $\mu$ */
    void Viscosity_c(TPZManVector<STATE,10> &mu, TPZManVector<STATE,10> &state_vars);
    
    /** @brief Compressibility - 1/pa $c$ */
    void Compressibility_c(TPZManVector<STATE,10> &c, TPZManVector<STATE,10> &state_vars);
    
    // @}
    
    /**
     * @defgroup Linearized models
     * @{
     */
    
    /** @brief Density - kg/m3  $\rho$ */
    void Density_l(TPZManVector<STATE,10> &rho, TPZManVector<STATE,10> &state_vars);
    
    /** @brief Viscosity - Pa s  $\mu$ */
    void Viscosity_l(TPZManVector<STATE,10> &mu, TPZManVector<STATE,10> &state_vars);
    
    /** @brief Compressibility - 1/pa $c$ */
    void Compressibility_l(TPZManVector<STATE,10> &c, TPZManVector<STATE,10> &state_vars);
    
    // @}
    
};

#endif /* TRMWaterPhase_h */
