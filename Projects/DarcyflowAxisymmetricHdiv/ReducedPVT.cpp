//
//  ReducedPVT.cpp
//  PZ
//
//  Created by Omar on 4/24/15.
//
//

#include "ReducedPVT.h"


ReducedPVT::ReducedPVT()
{
    
    /** @brief Temperature @ reservoir conditions  - F */
    fTRes = 180.0;
    
    /** @brief Pressure for references values - Pa */
    fPRef = 1.0;
    
    /** @brief Density - kg/m3  $\rho$ */
    fRho = 0.0;
    
    /** @brief viscosity - Pa s  $\mu$ */
    fMu = 0.0;
    
    /** @brief Compressibility - 1/pa $c$ */
    fc = 0.0;
    
}

ReducedPVT::~ReducedPVT()
{
    
}