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
void PetroPhysicData::Pcow(REAL Sw, REAL &Pcow, REAL &dPcowdSw)
{
    Pcow = 0.0;
    dPcowdSw = 0.0;
}

/** @brief Gas-Oil Capillary Pressure - Pa $P_{cgo}$ */
void PetroPhysicData::Pcgo(REAL So, REAL &Pcgo, REAL &dPcgodSo)
{
    Pcgo = 0.0;
    dPcgodSo = 0.0;
}

/** @brief Gas-Oil Capillary Pressure - Pa $P_{cgw}$ */
void PetroPhysicData::Pcgw(REAL Sw, REAL &Pcgw, REAL &dPcgwdSw)
{
    Pcgw = 0.0; // or Pcgo(So) + Pcow(Sw)
    dPcgwdSw = 0.0;
}


// Relative permeabilities

/** @brief Water Relative permeabilities  $K_{rw}$ */
void PetroPhysicData::Krw(REAL Sw, REAL &krw, REAL &dkrwdSw)
{
    krw = Sw;
    dkrwdSw = 1.0;
}

/** @brief Oil Relative permeabilities  $K_{ro}$ */
void PetroPhysicData::Kro(REAL So, REAL &kro, REAL &dkrodSo)
{
    kro = So;
    dkrodSo = 1.0;
}

/** @brief Gas Relative permeabilities  $K_{rg}$ */
void PetroPhysicData::Krg(REAL Sg, REAL &krg, REAL &dkrgdSw)
{
    krg = 0.0;
    dkrgdSw = 0.0;
}
