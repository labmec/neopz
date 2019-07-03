//
//  ReducedPVT.cpp
//  PZ
//
//  Created by Omar on 4/24/15.
//
//

#include "Phase.h"


Phase::Phase()
{
    
    /** @brief Temperature @ reservoir conditions 355.3722 [K] - 180 [F] */
    fTRes = 355.3722;
    
    /** @brief Temperature for references values [K] */
    fTRef = 1.0;
    
    /** @brief Pressure for references values - Pa */
    fPRef = 1.0;
    
    /** @brief Density - kg/m3  $\rho$ */
    fRho = 0.0;
    
    /** @brief viscosity - Pa s  $\mu$ */
    fMu = 0.0;
    
    /** @brief Compressibility - 1/pa $c$ */
    fc = 0.0;
    
    /** @brief Irreducible Saturation of the wetting phase - */
    fS_wett_r = 0.0;
    
    /** @brief Irreducible Saturation of the no wetting phase - */
    fS_nwett_r = 0.0;
    
    /** @brief Relative permeability model */
    fIsNonlinearKrQ = true;
    
}

Phase::~Phase()
{
    
}