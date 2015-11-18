//
//  GasPhase.h
//  PZ
//
//  Created by Omar on 9/30/15.
//
//

#ifndef __PZ__GasPhase__
#define __PZ__GasPhase__

#include <stdio.h>
#include "Phase.h"

#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include <math.h>


class GasPhase : public Phase

{
    
private:
    
    
    /** @brief Specific gravity of Gas fraction $ */
    REAL fGas_gamma;
    
    /** @brief Carbon dioxide content yCO2 molar fraction $ */
    REAL fyCO2;
    
    /** @brief Acid sulfhidric content yH2S molar fraction $ */
    REAL fyH2S;
    
    /** @brief Nitrogen content yN2 molar fraction $ */
    REAL fyN2;
    
    /** @brief Pressure at standard contidions [Pa] $ */
    REAL fPstd;
    
    /** @brief Temperature at standard contidions [K] $ */
    REAL fTstd;
    
    /** @brief Air mass density at standard contidions [kg/m3] $ */
    REAL fRho_air_std;
    
    
public:
    
    /** @brief Default constructor $ */
    GasPhase();
    
    /** @brief Default desconstructor $ */
    ~GasPhase();
    
    /** @brief Density - kg/m3  $\rho$ */
    void Density(TPZManVector<REAL> &rho, TPZManVector<REAL> state_vars);
    
    /** @brief viscosity - Pa s  $\mu$ */
    void Viscosity(TPZManVector<REAL> &mu, TPZManVector<REAL> state_vars);
    
    /** @brief viscosity standar - Pa s  $\mu$ */
    REAL Viscosity_std();
    
    /** @brief Compressibility - 1/pa $c$ */
    void Compressibility(TPZManVector<REAL> &c, TPZManVector<REAL> state_vars);
    
    /** @brief Relative permeability - $k_{r}$ */
    void Kr(TPZManVector<REAL> &kr, TPZManVector<REAL> state_vars);
    
    /** @brief Capillar pressure - $P_{c}$ */
    void Pc(TPZManVector<REAL> &pc, TPZManVector<REAL> state_vars);
    
    /** @brief Computes the pseudo critic pressure of Gas $ */
    REAL Ppc();
    
    /** @brief Computes the pseudo critic temperature of Gas $ */
    REAL Tpc();
    
    /** @brief Computes the compressibility factor using Beggs and Brill correlation (1974)  $ */
    void Z(TPZManVector<REAL> &z, TPZManVector<REAL> state_vars);
    
    /** @brief Computes the gas formation volume factor  $ */
    void Bg(TPZManVector<REAL> &bg, TPZManVector<REAL> state_vars);
    
    
};


#endif /* defined(__PZ__GasPhase__) */





