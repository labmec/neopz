//
//  TRMPhaseTransport.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMPhaseTransport__
#define __PZ__TRMPhaseTransport__

#include <stdio.h>
#include "pzmatwithmem.h"
#include "TRMPhaseMemory.h"
#include "TRMPhaseTransport.h"
#include "pzdiscgal.h"

class TRMPhaseTransport : public TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> {
    
    /// Implementar principalemente ContributeInterface, porque a saturacao no interior do elemento é zero (por enquanto)
};

#endif /* defined(__PZ__TRMPhaseTransport__) */
