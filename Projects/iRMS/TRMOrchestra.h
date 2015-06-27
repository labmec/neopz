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
#include "tpzautopointer.h"

// Type of structural matrices
#include "TPZSkylineNSymStructMatrix.h"

#include "TRMSimulationData.h"
#include "TRMFluxPressureAnalysis.h"
#include "TRMTransportAnalysis.h"
#include "TRMSpaceOdissey.h"

class TRMRawData;

class TRMOrchestra{
    
private:

    /** @brief Define the global geometry being used */
    TPZAutoPointer<TPZGeoMesh > fgmesh;

    /** @brief Define the mixed system analysis */
    TRMFluxPressureAnalysis fFluxPressureAnalysis;
    
    /** @brief Define the transpor equation analysis */
    TRMTransportAnalysis fTransportAnalysis;
    
    /** @brief Define simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
protected:
    
    /** @brief Solve the initial conditions for pressure using a l2 projection */
    void InitializePressure();
    
public:

    /** @brief Default constructor */
    TRMOrchestra();

    /** @brief Default desconstructor */
    ~TRMOrchestra();
    
    void SetGMesh(TPZAutoPointer<TPZGeoMesh > gmesh)
    {
        fgmesh = gmesh;
    }
    
    /** @brief Create computational meshes using space odissey */
    void CreateCompMeshes(TRMRawData &rawdata);
    
    /** @brief Create a primal analysis using space odissey */
    void CreateAnalysisPrimal(TRMSpaceOdissey spacegenerator);
    
    /** @brief Create a dual analysis using space odissey */
    void CreateAnalysisDual(TRMSpaceOdissey spacegenerator);
    
    /** @brief Run the time steps set in the simulation data */
    void RunSimulation();
    
    /** @brief Run a single time step */
    void ExecuteOneTimeStep();
    
    /** @brief Computes the post processed results */
    void PostProcess();
    
};

#endif /* defined(__PZ__TRMOrchestra__) */
