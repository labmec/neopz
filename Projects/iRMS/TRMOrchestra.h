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
    
    /** @brief The space generator */
    TRMSpaceOdissey fSpaceGenerator;
    
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
    void CreateAnalysisPrimal();

    /** @brief Create a dual analysis using space odissey */
    void CreateAnalysisDual();
    
    /** @brief Create a dual analysis on box geometry using space odissey */
    void CreateAnalysisDualonBox();
    
    /** @brief Run the time steps set in the simulation data */
    void RunSimulation();
    
    /// Transfer the flux solution to the saturation mesh
    void TransferToSaturationMesh();
    
    /** @brief Run a single time step */
    void ExecuteOneTimeStep();
    
    /** @brief Computes the post processed results */
    void PostProcess();
    
    /** @brief Project an exact solution */
    void ProjectExactSolution();
    
    /** @brief Compute the production rate of the reservoir */
    void ComputeProductionRate(std::map<REAL,STATE> &RatebyPosition, STATE &Total);
    
    /** @brief exact pressure */
    static void ExactPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &pressure);
    
    /** @brief exact flux */
    static void ExactFlux(const TPZVec<REAL> &pt, TPZVec<STATE> &flux);
    
    /** @brief exact laplacian */
    static void ExactLaplacian(const TPZVec<REAL> &pt, TPZVec<STATE> &pressure);
    
    /** @brief Compute the system of equations using transfer matrixces */
    TPZFMatrix<STATE> IntegrateResidue(TPZAutoPointer<TPZCompMesh> cmesh_multiphysics, TPZAutoPointer< TPZCompMesh> cmesh_flux, TPZAutoPointer< TPZCompMesh> cmesh_pressure, TPZAutoPointer<TRMBuildTransfers> transfer);
    
    /** @brief Compute gradient of the system of equations using transfer matrixces */
    void IntegrateGradientOfResidue(TPZAutoPointer< TPZCompMesh> cmesh_flux, TPZAutoPointer<TRMBuildTransfers> transfer);
    
};

#endif /* defined(__PZ__TRMOrchestra__) */
