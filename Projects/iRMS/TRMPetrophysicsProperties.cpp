//
//  TRMPetrophysicsProperties.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMPetrophysicsProperties.h"


TRMPetrophysicsProperties::TRMPetrophysicsProperties()
{
    
    
}

TRMPetrophysicsProperties::~TRMPetrophysicsProperties()
{
    
}


// Capillary Pressure

/** @brief Oil-Water Capillary Pressure - Pa $P_{cow}$ */
void TRMPetrophysicsProperties::Pcow(REAL Sw, REAL &Pcow, REAL &dPcowdSw)
{
    Pcow = 0.0;
    dPcowdSw = 0.0;
}

/** @brief Gas-Oil Capillary Pressure - Pa $P_{cgo}$ */
void TRMPetrophysicsProperties::Pcgo(REAL So, REAL &Pcgo, REAL &dPcgodSo)
{
    Pcgo = 0.0;
    dPcgodSo = 0.0;
}

/** @brief Gas-Oil Capillary Pressure - Pa $P_{cgw}$ */
void TRMPetrophysicsProperties::Pcgw(REAL Sw, REAL &Pcgw, REAL &dPcgwdSw)
{
    Pcgw = 0.0; // or Pcgo(So) + Pcow(Sw)
    dPcgwdSw = 0.0;
}


// Relative permeabilities

/** @brief Water Relative permeabilities  $K_{rw}$ */
void TRMPetrophysicsProperties::Krw(REAL Sw, REAL &krw, REAL &dkrwdSw)
{
    krw = Sw;
    dkrwdSw = 1;
    
}

/** @brief Oil Relative permeabilities  $K_{ro}$ */
void TRMPetrophysicsProperties::Kro(REAL So, REAL &kro, REAL &dkrodSo)
{
    kro = So;
    dkrodSo = 1;
    
}

/** @brief Gas Relative permeabilities  $K_{rg}$ */
void TRMPetrophysicsProperties::Krg(REAL Sg, REAL &krg, REAL &dkrgdSw)
{
    krg = 0.0;
    dkrgdSw = 0.0;
}
