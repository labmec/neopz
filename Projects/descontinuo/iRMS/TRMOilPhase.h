//
//  TRMOilPhase.hpp
//  PZ
//
//  Created by Omar on 5/23/16.
//
//

#ifndef TRMOilPhase_hpp
#define TRMOilPhase_hpp

#include <stdio.h>
#include "TRMPhaseProperties.h"

class TRMOilPhase : public TRMPhaseProperties {
    
public:
    
    /** @brief default constructor */
    TRMOilPhase();
    
    /** @brief default constructor */
    TRMOilPhase(const TRMOilPhase &copy)
    {
        DebugStop();
    }
    
    /** @brief default constructor */
    TRMOilPhase &operator=(const TRMOilPhase &copy)
    {
        DebugStop();
        return *this;
    }
    

    /** @brief default destructor */
    ~TRMOilPhase();
    
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

#endif /* TRMOilPhase_hpp */
