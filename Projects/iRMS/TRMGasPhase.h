//
//  TRMGasPhase.hpp
//  PZ
//
//  Created by Omar on 5/23/16.
//
//

#ifndef TRMGasPhase_hpp
#define TRMGasPhase_hpp

#include <stdio.h>
#include "TRMPhaseProperties.h"

class TRMGasPhase : public TRMPhaseProperties {
    
public:
    
    /** @brief default constructor */
    TRMGasPhase();
    
    /** @brief default destructor */
    ~TRMGasPhase();
    
    TRMGasPhase(const TRMGasPhase &copy)
    {
        DebugStop();
    }
    
    TRMGasPhase &operator=(const TRMGasPhase &copy)
    {
        DebugStop();
        return *this;
    }
    
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

#endif /* TRMGasPhase_hpp */
