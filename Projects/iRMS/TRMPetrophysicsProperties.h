//
//  TRMPetrophysicsProperties.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMPetrophysicsProperties__
#define __PZ__TRMPetrophysicsProperties__

#include <stdio.h>
#include "pzreal.h"

class TRMPetrophysicsProperties{
    
    /**
     * @ingroup Petrophysics Properties
     * @brief Capillary Pressure and relative permeabilities models
     * @since June 09, 2015
     */
    
public:
    
    TRMPetrophysicsProperties();
    
    ~TRMPetrophysicsProperties();
    
    TRMPetrophysicsProperties(const TRMPetrophysicsProperties &copy)
    {
        DebugStop();
    }
    
    TRMPetrophysicsProperties &operator=(const TRMPetrophysicsProperties &copy)
    {
        DebugStop();
        return *this;
    }
    
    // Capillary Pressure models
    
    /** @brief Oil-Water Capillary Pressure - Pa $P_{cow}$ */
    void Pcow(REAL Sw, REAL &Pcow, REAL &dPcowdSw);
    
    /** @brief Gas-Oil Capillary Pressure - Pa $P_{cgo}$ */
    void Pcgo(REAL So, REAL &Pcgo, REAL &dPcgodSo);
    
    /** @brief Gas-Oil Capillary Pressure - Pa $P_{cgw}$ */
    void Pcgw(REAL Sw, REAL &Pcgw, REAL &dPcgwdSw);
    
    
    // Relative permeabilities models
    
    /** @brief Oil-Water Capillary Pressure - Pa $K_{rw}$ */
    void Krw(REAL Sw, REAL &krw, REAL &dkrwdSw);
    
    /** @brief Gas-Oil Capillary Pressure - Pa $K_{ro}$ */
    void Kro(REAL So, REAL &kro, REAL &dkrodSo);
    
    /** @brief Gas-Oil Capillary Pressure - Pa $K_{rg}$ */
    void Krg(REAL Sg, REAL &krg, REAL &dkrgdSw);
    
};

#endif /* defined(__PZ__TRMPetrophysicsProperties__) */
