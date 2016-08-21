//
//  TRMSpaceOdissey.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
// This class defines the proper approximation space in accordance to the model. Steady state, parabolic, MHM, MHM++.

#ifndef __PZ__TRMSpaceOdissey__
#define __PZ__TRMSpaceOdissey__

#include <stdio.h>
#include "tpzautopointer.h"
#include "TRMSimulationData.h"
#include "TRMRawData.h"
#include "TRMBuildTransfers.h"
#include "pzgmesh.h"
#include "pzcmesh.h"

#include "pzbuildmultiphysicsmesh.h"
#include "pzcondensedcompel.h"

#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"



class TRMSpaceOdissey{
    
public:
    
    /** @brief Define the type of geometry being used */
    enum MGeoMeshType {ENone = 0,EBox = 1, EReservoir = 2, EGIDReservoir = 3};
    
    /** @brief Define the type of geometry being used */
    MGeoMeshType fMeshType;
    
private:
    
    /** @brief order of approximation flux and pressure */
    int fPOrder;
    
    /** @brief order of approximation transport */
    int fSOrder;
    
    /** @brief Autopointer of the Geometric mesh shared with all the classes involved */
    TPZGeoMesh * fGeoMesh;
    
    /** @brief Autopointer of Simulation data */
    TRMSimulationData * fSimulationData;
    
    /** @brief H1 computational mesh for primal approach */
    TPZCompMesh * fH1Cmesh;
    
    /** @brief Hdiv computational mesh, conservative vector field */
    TPZCompMesh * fFluxCmesh;
    
    /** @brief L2 computational mesh the restriction equation */
    TPZCompMesh * fPressureCmesh;
    
    /** @brief L2 computational mesh alpha saturation equations */
    TPZCompMesh * fAlphaSaturationMesh;

    /** @brief L2 computational mesh beta saturation equations */
    TPZCompMesh * fBetaSaturationMesh;
    
    /** @brief H1 computational mesh for Maurice Biot linear poroelasticity */
    TPZCompMesh * fGeoMechanicsCmesh;
    
    /** @brief L2 computational multiphysics mesh transport equations */
    TPZCompMesh * fTransportMesh;
    
    /** @brief Mixed computational mesh for a dual analysis */
    TPZCompMesh * fMixedFluxPressureCmesh;
    
    /** @brief H1-L2 for computational mesh fora primal analysis with Global postprocessing of fluxes */
    TPZCompMesh * fPressureSaturationCmesh;
    
    /** @brief Computational mesh for multiphase monolithic approach */
    TPZCompMesh * fMonolithicMultiphaseCmesh;
    
    void ModifyElementOrders(std::map<long,int> &elorders);

public:
    
    /** @brief Default constructor */
    TRMSpaceOdissey();
    
    /** @brief Default constructor */
    TRMSpaceOdissey(const TRMSpaceOdissey &copy)
    {
        DebugStop();
    }
    
    /** @brief Default constructor */
    TRMSpaceOdissey &operator=(const TRMSpaceOdissey &copy)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief Default desconstructor */
    ~TRMSpaceOdissey();
    
    /** @brief Initialize the simulation data */
    void InitializeSimulationData(TRMRawData &rawdata);
    
    /** @brief Set the order of approximation flux and pressure */
    void SetDefaultPOrder(int porder)
    {
        fPOrder = porder;
    }
    
    /** @brief Set the order of approximation transport */
    void SetDefaultSOrder(int sorder)
    {
        fSOrder = sorder;
    }
    
    /** @brief Create a H1 computational mesh */
    void CreateH1Cmesh();
    
    /** @brief Create a Hdiv computational mesh */
    void CreateFluxCmesh();

    /** @brief Create a Discontinuous computational mesh L2 */
    void CreatePressureCmesh();
    
    /** @brief Create a Mixed computational mesh Hdiv-L2 */
    void CreateMixedCmesh();
    
    /** @brief Create a Mixed-Transport multiphase computational mesh Hdiv-L2-L2-L2 */
    void CreateMultiphaseCmesh();
    
    /** @brief Statically condense the internal equations of the elements */
    void StaticallyCondenseEquations();
    
    /** @brief Configure the boundary conditions of a well with reservoir boundary conditions */
    void ConfigureWellConstantPressure(STATE wellpressure, STATE farfieldpressure);
    
    /** @brief Create a computational mesh L2 of alpha phase */
    void CreateAlphaTransportMesh();
    
    /** @brief Create a computational mesh L2 of beta phase  */
    void CreateBetaTransportMesh();
    
    /** @brief Create a multiphysics computational mesh L2 */
    void CreateTransportMesh();
    
    /** @brief Create a H1 computational mesh for Maurice Biot Linear Poroelasticity */
    void CreateGeoMechanicMesh();
    
    /** @brief Create the reservoir geometry */
    void CreateGeometricReservoirMesh();
    
    /** @brief Print the reservoir geometry */
    void PrintGeometry();
    
    /** @brief Create a reservoir-box geometry */
    void CreateGeometricGIDMesh(std::string &grid);
    
    /** @brief Create a reservoir-box geometry */
    void CreateGeometricBoxMesh2D(TPZManVector<REAL,2> dx, TPZManVector<REAL,2> dy);
    
    /** @brief Create a reservoir-box geometry */
    void CreateGeometricBoxMesh(TPZManVector<REAL,2> dx, TPZManVector<REAL,2> dy, TPZManVector<REAL,2> dz);
    
    /** @brief Create a reservoir-box geometry */
    void CreateGeometricExtrudedGIDMesh(std::string &grid, TPZManVector<REAL,2> dz);
    
    /** @brief Parametric function that computes elements in the x direction */
    static  void ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /** @brief Parametric function that computes elements in the y direction */
    static  void ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /** @brief Parametric function that computes elements in the z direction */
    static  void ParametricfunctionZ(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /** @brief Set autopointer of the Geometric mesh */
    void SetGmesh(TPZGeoMesh *GeoMesh){
        fGeoMesh = GeoMesh;
    }
    
    /** @brief Get autopointer of the Geometric mesh */
    TPZGeoMesh * Gmesh(){
        return fGeoMesh;
    }
    
    /** @brief Apply uniform refinement on the Geometric mesh */
    void UniformRefinement(int n_ref);
    
    /** @brief Set autopointer of Simulation data */
    void SetSimulationData(TRMSimulationData * SimulationData){
        fSimulationData = SimulationData;
    }
    
    /** @brief Get autopointer of Simulation data */
    TRMSimulationData *  SimulationData(){
        return fSimulationData;
    }
    
    /** @brief Set autopointer of H1 computational mesh for primal approach */
    void SetH1CMesh(TPZCompMesh * H1Cmesh){
        fH1Cmesh = H1Cmesh;
    }
    
    /** @brief Get autopointer of H1 computational mesh for primal approach */
    TPZCompMesh *  H1CMesh(){
        return fH1Cmesh;
    }
    
    /** @brief Set autopointer of Hdiv computational mesh conservative vector field */
    void SetFluxCmesh(TPZCompMesh * FluxCmesh){
        fFluxCmesh = FluxCmesh;
    }
    
    /** @brief Get autopointer of Hdiv computational mesh conservative vector field */
    TPZCompMesh * FluxCmesh(){
        return fFluxCmesh;
    }
    
    /** @brief Set L2 computational multiphysics mesh transport equations */
    void SetTransportMesh(TPZCompMesh * FluxCmesh_Int){
        fTransportMesh = FluxCmesh_Int;
    }
    
    /** @brief Get L2 computational multiphysics mesh transport equations */
    TPZCompMesh * TransportMesh(){
        return fTransportMesh;
    }
    
    /** @brief Set autopointer of L2 computational mesh the restriction equation */
    void SetPressureCmesh(TPZCompMesh * PressureCmesh){
        fPressureCmesh = PressureCmesh;
    }
    
    /** @brief Get autopointer of L2 computational mesh the restriction equation */
    TPZCompMesh * PressureCmesh(){
        return fPressureCmesh;
    }
    
    /** @brief Set autopointer of L2 computational mesh water saturation equations */
    void SetWaterSaturationMesh(TPZCompMesh * AlphaSaturationMesh){
        fAlphaSaturationMesh = AlphaSaturationMesh;
    }
    
    /** @brief Get autopointer of L2 computational mesh water saturation equations */
    TPZCompMesh * AlphaSaturationMesh(){
        return fAlphaSaturationMesh;
    }
    
    /** @brief Set autopointer of L2 computational mesh oil saturation equations */
    void SetBetaSaturationMesh(TPZCompMesh * BetaSaturationMesh){
        fBetaSaturationMesh = BetaSaturationMesh;
    }
    
    /** @brief Get autopointer of L2 computational mesh oil saturation equations */
    TPZCompMesh * BetaSaturationMesh(){
        return fBetaSaturationMesh;
    }
    
    /** @brief Set autopointer of H1 computational mesh for Maurice Biot linear poroelasticity */
    void SetGeoMechanicsCmesh(TPZCompMesh * GeoMechanicsCmesh){
        fGeoMechanicsCmesh = GeoMechanicsCmesh;
    }
    
    /** @brief Get autopointer of H1 computational mesh for Maurice Biot linear poroelasticity */
    TPZCompMesh * GeoMechanicsCmesh(){
        return fGeoMechanicsCmesh;
    }
    
    /** @brief Set autopointer of Mixed computational mesh for a dual analysis */
    void SetMixedFluxPressureCmesh(TPZCompMesh * MixedFluxPressureCmesh){
        fMixedFluxPressureCmesh = MixedFluxPressureCmesh;
    }
    
    /** @brief Get autopointer of Mixed computational mesh for a dual analysis */
    TPZCompMesh * MixedFluxPressureCmesh(){
        return fMixedFluxPressureCmesh;
    }
    
    /** @brief Set autopointer of H1-L2 for computational mesh fora primal analysis with Global postprocessing of fluxes */
    void SetPressureSaturationCmesh(TPZCompMesh * PressureSaturationCmesh){
        fPressureSaturationCmesh = PressureSaturationCmesh;
    }
    
    /** @brief Get autopointer of H1-L2 for computational mesh fora primal analysis with Global postprocessing of fluxes */
    TPZCompMesh * PressureSaturationCmesh(){
        return fPressureSaturationCmesh;
    }
    
    /** @brief Set autopointer of Computational mesh for multiphase monolithic approach */
    void SetMonolithicMultiphaseCmesh(TPZCompMesh * MonolithicMultiphaseCmesh){
        fMonolithicMultiphaseCmesh = MonolithicMultiphaseCmesh;
    }
    
    /** @brief Get autopointer of Computational mesh for multiphase monolithic approach */
    TPZCompMesh * MonolithicMultiphaseCmesh(){
        return fMonolithicMultiphaseCmesh;
    }
    
    /** @brief Create computational interfaces inside for flux computations  */
    void CreateInterfacesInside(TPZCompMesh * cmesh);
    
    /** @brief Adjust the polinomial order of the elements */
    void IncreaseOrderAroundWell(int numlayers);
    
    
    
};

#endif /* defined(__PZ__TRMSpaceOdissey__) */
