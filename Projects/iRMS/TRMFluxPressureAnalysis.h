//
//  TRMFluxPressureAnalysis.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMFluxPressureAnalysis__
#define __PZ__TRMFluxPressureAnalysis__

#include <stdio.h>
#include "pzanalysis.h"
#include "pzstepsolver.h"

class TRMTranportAnalysis;

class TRMFluxPressureAnalysis : public TPZAnalysis {
    
    // Setting up data
    // Create a struture for nonlinear analysis
    // required a initialization process of the fPreconditioner with the decomposition fo the first linear system
    // it will simply converges the residuum in one time step considering fixed saturation (inverse for transport analysis)

    TPZStepSolver<STATE> fPreconditioner;
    
    void UpdateTransientMemory(); // updates the transient memory
    
    void UpdateSolution(TRMTranportAnalysis &transportAnalysis);
    
};

#endif /* defined(__PZ__TRMFluxPressureAnalysis__) */
