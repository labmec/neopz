/*
 *  ReservoirData.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "PetroPhysicData.h"

PetroPhysicData::PetroPhysicData()
{
 
    
}

PetroPhysicData::~PetroPhysicData()
{
    
}

// Capillary Pressure

/** @brief Oil-Water Capillary Pressure - Pa $P_{cow}$ */
void PetroPhysicData::Pcwo(TPZVec<REAL> &pc_wo, TPZVec<REAL> state_vars)
{
    
//    REAL so = state_vars[2];
//    REAL sw = state_vars[3];
//    REAL sg = 1.0 - state_vars[2] - state_vars[3];
    
    pc_wo[0] = 0.0;
    pc_wo[3] = 0.0;
    pc_wo[4] = 0.0;
}

/** @brief Gas-Oil Capillary Pressure - Pa $P_{cgo}$ */
void PetroPhysicData::Pcog(TPZVec<REAL> &pc_og, TPZVec<REAL> state_vars)
{
    
//    REAL so = state_vars[2];
//    REAL sw = state_vars[3];
//    REAL sg = 1.0 - state_vars[2] - state_vars[3];
    
    pc_og[0] = 0.0;
    pc_og[3] = 0.0;
    pc_og[4] = 0.0;
}

/** @brief Gas-Oil Capillary Pressure - Pa $P_{cgw}$ */
void PetroPhysicData::Pcwg(TPZVec<REAL> &pc_wg, TPZVec<REAL> state_vars)
{
    
//    REAL so = state_vars[2];
//    REAL sw = state_vars[3];
//    REAL sg = 1.0 - state_vars[2] - state_vars[3];
    
    pc_wg[0] = 0.0;// or Pcgo(So) + Pcow(Sw)
    pc_wg[3] = 0.0;
    pc_wg[4] = 0.0;

}


// Relative permeabilities

/** @brief Water Relative permeabilities  $K_{rw}$ */
void PetroPhysicData::Krw(TPZVec<REAL> &kr_w, TPZVec<REAL> state_vars)
{
    
//    REAL so = state_vars[2];
    REAL sw = state_vars[3];
//    REAL sg = 1.0 - state_vars[2] - state_vars[3];
    
    kr_w[0] = sw;
    kr_w[3] = 0.0;
    kr_w[4] = 1.0;
  
}

/** @brief Oil Relative permeabilities  $K_{ro}$ */
void PetroPhysicData::Kro(TPZVec<REAL> &kr_o, TPZVec<REAL> state_vars)
{
    
    REAL so = state_vars[2];
//    REAL sw = state_vars[3];
//    REAL sg = 1.0 - state_vars[2] - state_vars[3];
    
    kr_o[0] = so;
    kr_o[3] = 1.0;
    kr_o[4] = 0.0;
  
}

/** @brief Gas Relative permeabilities  $K_{rg}$ */
void PetroPhysicData::Krg(TPZVec<REAL> &kr_g, TPZVec<REAL> state_vars)
{
    REAL so = state_vars[2];
    REAL sw = state_vars[3];
//    REAL sg = 1.0 - state_vars[2] - state_vars[3];
    
    kr_g[0] = 1-so-sw;
    kr_g[3] = -1.0;
    kr_g[4] = -1.0;
}
