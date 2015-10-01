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
    void Pcwo(TPZVec<REAL> &pc_wo, TPZVec<REAL> state_vars);
    
    /** @brief Gas-Oil Capillary Pressure - Pa $P_{cgo}$ */
    void Pcog(TPZVec<REAL> &pc_og, TPZVec<REAL> state_vars);
    
    /** @brief Gas-Oil Capillary Pressure - Pa $P_{cgw}$ */
    void Pcwg(TPZVec<REAL> &pc_wg, TPZVec<REAL> state_vars);
    
    
    // Relative permeabilities
    
    /** @brief Oil-Water Capillary Pressure - Pa $K_{rw}$ */
    void Krw(TPZVec<REAL> &kr_w, TPZVec<REAL> state_vars);
    
    /** @brief Gas-Oil Capillary Pressure - Pa $K_{ro}$ */
    void Kro(TPZVec<REAL> &kr_o, TPZVec<REAL> state_vars);
    
    /** @brief Gas-Oil Capillary Pressure - Pa $K_{rg}$ */
    void Krg(TPZVec<REAL> &kr_g, TPZVec<REAL> state_vars);
    
};


#endif