//
//  TRMOrchestra.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMOrchestra__
#define __PZ__TRMOrchestra__

#include <stdio.h>
#include "pzgmesh.h"
#include "TRMSimulationData.h"
#include "TRMFluxPressureAnalysis.h"
#include "TRMTransportAnalysis.h"

class TRMRawData;

class TRMOrchestra{
    
private:
    // Defines the global geometry
    TPZGeoMesh fgmesh;
    
    TRMFluxPressureAnalysis fFluxPressureAnalysis;
    TRMTransportAnalysis fTransportAnalysis;
    
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
protected:
    
    void InitializePressure(); //  Solve the initial conditions for pressure using a l2 projection
    
public:
    
    void SetGMesh(TPZGeoMesh &gmesh){
        fgmesh = gmesh;
    }
    
    
    void CreateCompMeshes(TRMRawData &rawdata); //will create comp meshes using space odissey. RawData may or may not be needed
    
    void RunSimulation(); // Run the time steps set in the simulation data
    
    void ExecuteOneTimeStep();
    
    void PostProcess();
    
};

#endif /* defined(__PZ__TRMOrchestra__) */
