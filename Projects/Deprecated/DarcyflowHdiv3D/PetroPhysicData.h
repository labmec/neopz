#ifndef TPZPetroPhysicDATAH
#define TPZPetroPhysicDATAH
/*
 *  ReservoirData.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include <math.h>

class PetroPhysicData {
    
public:
    
    /**
     * @ingroup Characteristic Parameters
     * @brief Define characteristic parameters for Darcy linear flow.
     * @since December 08, 2014
     */
    
    
    PetroPhysicData();
    
    ~PetroPhysicData();
    
    // Capillary Pressure
    
    /** @brief Oil-Water Capillary Pressure - Pa $P_{cow}$ */
    void Pcow(REAL Sw, REAL &Pcow, REAL &dPcowdSw);
    
    /** @brief Gas-Oil Capillary Pressure - Pa $P_{cgo}$ */
    void Pcgo(REAL So, REAL &Pcgo, REAL &dPcgodSo);
    
    /** @brief Gas-Oil Capillary Pressure - Pa $P_{cgw}$ */
    void Pcgw(REAL Sw, REAL &Pcgw, REAL &dPcgwdSw);
    
    
    // Relative permeabilities
    
    /** @brief Oil-Water Capillary Pressure - Pa $K_{rw}$ */
    void Krw(REAL Sw, REAL &krw, REAL &dkrwdSw);
    
    /** @brief Gas-Oil Capillary Pressure - Pa $K_{ro}$ */
    void Kro(REAL So, REAL &kro, REAL &dkrodSo);
    
    /** @brief Gas-Oil Capillary Pressure - Pa $K_{rg}$ */
    void Krg(REAL Sg, REAL &krg, REAL &dkrgdSw);
    
};


#endif