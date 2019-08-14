//
//  TRMPhaseProperties.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMPhaseProperties.h"


/** @brief default constructor */
TRMPhaseProperties::TRMPhaseProperties(){
    
    /** @brief Temperature @ reservoir conditions 355.3722 [K] - 180 [F] */
    fTRes = 355.3722;
    
    /** @brief Temperature for references values [K] - 20 [C] */
    fTRef = 293.15;
    
    /** @brief Pressure for references values - Pa - 1 [atm] */
    fPRef = 101325.0;
    
    /** @brief Density - kg/m3  $\rho$ */
    frho.resize(0);
    
    /** @brief viscosity - Pa s  $\mu$ */
    fmu.resize(0);
    
    /** @brief Compressibility - 1/pa $c$ */
    fc.resize(0);
    
    /** @brief Density model = {0,1,2}  */
    frho_model = 0;
    
    /** @brief Viscosity model = {0,1,2} */
    fmu_model = 0;
    
    /** @brief Compressibility model = {0,1,2} */
    fc_model = 0;
    
}

/** @brief default destructor */
TRMPhaseProperties::~TRMPhaseProperties(){
    
}
