//
//  TRMTransportAnalysis.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMTransportAnalysis__
#define __PZ__TRMTransportAnalysis__

#include <stdio.h>

class TRMFluxPressureAnalysis;

class TRMTransportAnalysis{
    
public:
    void UpdateSolution(TRMFluxPressureAnalysis &fluxPressureAnalysis); // no pior dos casos pega o ponto de integracao chama o solution e atualiza o ponto
    
};

#endif /* defined(__PZ__TRMTransportAnalysis__) */