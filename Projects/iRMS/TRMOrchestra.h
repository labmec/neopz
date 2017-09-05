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
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"

#include "TRMSpaceOdissey.h"
#include "TRMSimulationData.h"
#include "TRMPrimalMultiphaseAnalysis.h"
#include "TRMMonolithicMultiphaseAnalysis.h"
#include "TRMSegregatedAnalysis.h"

#include "TRMFluxPressureAnalysis.h"
#include "TRMTransportAnalysis.h"
#include "TPZCompMeshTools.h"


class TRMRawData;

class TRMOrchestra{
    
private:

    /** @brief Define the global geometry being used */
    TPZGeoMesh * fgmesh;

    /** @brief Define the space generator */
    TRMSpaceOdissey * fSpaceGenerator;
    
    /** @brief Define simulation data */
    TRMSimulationData * fSimulationData;
    
    /** @brief Define the primal with global post-processing analysis */
    TRMPrimalMultiphaseAnalysis * fPrimalMultiphaseAnalysis;

    /** @brief Define the initial monolithic multiphase analysis */
    TRMMonolithicMultiphaseAnalysis * fMonolithicMultiphaseAnalysis_I;
    
    /** @brief Define the monolithic multiphase analysis */
    TRMMonolithicMultiphaseAnalysis * fMonolithicMultiphaseAnalysis;
    
//    /** @brief Define the analysis mixed system  */
//    TPZAutoPointer<TRMFluxPressureAnalysis> fparabolic;
//    
//    /** @brief Define the analysis for transport phases */
//    TPZAutoPointer<TRMTransportAnalysis> fhyperbolic;
    
    /** @brief Define initial multiphase segregated analysis for transport phases */
    TRMSegregatedAnalysis * fSegregatedAnalysis_I;
    
    /** @brief Define multiphase segregated analysis for transport phases */
    TRMSegregatedAnalysis * fSegregatedAnalysis;
    
    
    /** @brief Define the use of monolithic analysis  */
    bool fIsMonolithicQ;
    
    /** @brief Define the use of segregated analysis */
    bool fIsSegregatedQ;

    /** @brief Define the use of segregated with GC plus projection analysis */
    bool fIsSegregatedwithCGQ;
    
protected:
    
    /** @brief Solve the initial conditions for pressure using a l2 projection */
    void InitializePressure();
    
public:

    /** @brief Default constructor */
    TRMOrchestra();

    /** @brief Default desconstructor */
    ~TRMOrchestra();
    
    /** @brief Default constructor */
    TRMOrchestra(const TRMOrchestra &copy)
    {
        DebugStop();
    }
    
    /** @brief Default constructor */
    TRMOrchestra &operator=(const TRMOrchestra &copy)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief Set autopointer of the global geometry being used */
    void SetGMesh(TPZGeoMesh * gmesh)
    {
        fgmesh = gmesh;
    }
    /** @brief Get autopointer of the global geometry being used */
    TPZAutoPointer<TPZGeoMesh > GMesh()
    {
        return fgmesh;
    }
    
    /** @brief Set the space generator */
    void SetSpaceGenerator(TRMSpaceOdissey * SpaceGenerator)
    {
        fSpaceGenerator = SpaceGenerator;
    }
    /** @brief Get the space generator */
    TRMSpaceOdissey * SpaceGenerator()
    {
        return fSpaceGenerator;
    }
    
    /** @brief Set autopointer of the simulation data */
    void SetSimulationData(TRMSimulationData * SimulationData)
    {
        fSimulationData = SimulationData;
        fSpaceGenerator->SetSimulationData(SimulationData);
    }
    /** @brief Get autopointer of the simulation data */
    TRMSimulationData * SimulationData()
    {
        return fSimulationData;
    }
    
    /** @brief Set autopointer of the primal with global post-processing analysis */
    void SetPrimalMultiphaseAnalysis(TRMPrimalMultiphaseAnalysis * PrimalMultiphaseAnalysis)
    {
        fPrimalMultiphaseAnalysis = PrimalMultiphaseAnalysis;
    }
    /** @brief Get autopointer of the primal with global post-processing analysis */
    TRMPrimalMultiphaseAnalysis * PrimalMultiphaseAnalysis()
    {
        return fPrimalMultiphaseAnalysis;
    }

    /** @brief Set autopointer of the initial monolithic multiphase analysis */
    void SetMonolithicMultiphaseAnalysis_I(TRMMonolithicMultiphaseAnalysis * MonolithicMultiphaseAnalysis_I)
    {
        fMonolithicMultiphaseAnalysis_I = MonolithicMultiphaseAnalysis_I;
    }
    /** @brief Get autopointer of the initial monolithic multiphase analysis */
    TRMMonolithicMultiphaseAnalysis * MonolithicMultiphaseAnalysis_I()
    {
        return fMonolithicMultiphaseAnalysis_I;
    }
    
    /** @brief Set autopointer of the monolithic multiphase analysis */
    void SetMonolithicMultiphaseAnalysis(TRMMonolithicMultiphaseAnalysis * MonolithicMultiphaseAnalysis)
    {
        fMonolithicMultiphaseAnalysis = MonolithicMultiphaseAnalysis;
    }
    /** @brief Get autopointer of the monolithic multiphase analysis */
    TRMMonolithicMultiphaseAnalysis * MonolithicMultiphaseAnalysis()
    {
        return fMonolithicMultiphaseAnalysis;
    }
    
//    /** @brief Set autopointer of the mixed system analysis */
//    void SetFluxPressureAnalysis(TPZAutoPointer<TRMFluxPressureAnalysis > &FluxPressureAnalysis)
//    {
//        fFluxPressureAnalysis = FluxPressureAnalysis;
//    }
//    /** @brief Get autopointer of the mixed system analysis */
//    TPZAutoPointer<TRMFluxPressureAnalysis > FluxPressureAnalysis()
//    {
//        return fFluxPressureAnalysis;
//    }
//    
//    /** @brief Set autopointer of the transpor equations analysis */
//    void SetTransportAnalysis(TPZAutoPointer<TRMTransportAnalysis > &TransportAnalysis)
//    {
//        fTransportAnalysis = TransportAnalysis;
//    }
//    /** @brief Get autopointer of the water transpor equation analysis */
//    TPZAutoPointer<TRMTransportAnalysis > TransportAnalysis()
//    {
//        return fTransportAnalysis;
//    }
    
    
    /** @brief Set the use of monolithic analysis  */
    void SetMonolithicQ(bool IsMonolithicQ){
        fIsMonolithicQ = IsMonolithicQ;
    }
    
    /** @brief Get the use of monolithic analysis  */
    bool IsMonolithicQ(){
        return fIsMonolithicQ;
    }
    
    /** @brief Set the use of segregated analysis  */
    void SetSegregatedQ(bool IsSegregatedQ){
        fIsSegregatedQ = IsSegregatedQ;
    }
    
    /** @brief Get the use of segregated analysis  */
    bool IsSegregatedQ(){
        return fIsSegregatedQ;
    }

    /** @brief Set the use of segregated with GC plus projection analysis  */
    void SetSegregatedwithCGQ(bool IsSegregatedwithCGQ){
        fIsSegregatedwithCGQ = IsSegregatedwithCGQ;
    }
    
    /** @brief Get the use of segregated with GC plus projection analysis  */
    bool IsSegregatedwithCGQ(){
        return fIsSegregatedwithCGQ;
    }
    
    /** @brief Create geometric mesh being used by space odissey */
    void BuildGeometry(bool Is3DGeometryQ);
    
    /** @brief Create geometric mesh being used by space odissey */
    void BuildGeometry2();
    
    /** @brief Create computational meshes using space odissey */
    void CreateCompMeshes();
    
    /** @brief Create a primal analysis using space odissey */
    void CreateAnalysisPrimal();

    /** @brief Create a dual analysis using space odissey */
    void CreateAnalysisDual();
    
    /** @brief Create a dual analysis on box geometry using space odissey */
    void CreateAnalysisDualonBox(bool IsInitialQ);
    
    /** @brief Create a monolithic dual analysis on box geometry using space odissey */
    void CreateMonolithicAnalysis(bool IsInitialQ);

    /** @brief Run the static problem for a single step with a large time step */
    void RunStaticProblem();
    
    /** @brief Run the evolutionary problem for all steps set in the simulation data */
    void RunEvolutionaryProblem();
    
    /** @brief Run a single time step */
    void ExecuteOneTimeStep();
    
    /** @brief Must report time */
    bool MustResporTimeQ(REAL time, bool & draw_mixed_mapQ);
    
    /** @brief Computes the post processed results */
    void PostProcess();

    /** @brief Computes the post processed results */
    void ComputationalMeshUniformRefinement(TPZCompMesh  *cmesh, int ndiv);
    
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
    
    // @}    
    
};

#endif /* defined(__PZ__TRMOrchestra__) */
