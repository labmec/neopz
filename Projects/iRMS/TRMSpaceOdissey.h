//
//  TRMSpaceOdissey.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMSpaceOdissey__
#define __PZ__TRMSpaceOdissey__

#include <stdio.h>
#include "tpzautopointer.h"
#include "TRMSimulationData.h"

class TRMSpaceOdissey{
    
    // Create the computational meshes
    // Initializate the TRMSimulationData
    
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
public:
    
    void InitializeSimulationData(TRMRawData &rawdata);
    
    void CreateFluxPressureMesh();
    
    void CreateTransportMesh();
    
    void CreateGeoMechanicMesh();
    
    
    
};

#endif /* defined(__PZ__TRMSpaceOdissey__) */
