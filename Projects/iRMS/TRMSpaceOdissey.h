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
#include "TRMGmshReader.h"
#include "TPZVTKGeoMesh.h"


#include "pzpoisson3d.h"
#include "TRMMultiphase.h"
#include "TRMMixedDarcy.h"
#include "TRMBiotPoroelasticity.h"
#include "TPZMatLaplacian.h"
#include "pzbndcond.h"
#include "TRMPhaseTransport.h"
#include "TRMPhaseInterfaceTransport.h"

#include "pzelementgroup.h"
#include "TPZCompMeshTools.h"
#include "pzsubcmesh.h"
#include "pzcheckgeom.h"

#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternTools.h"
#include "tpzhierarquicalgrid.h"
#include "pzgeopoint.h"
#include "TRMSimworxMeshGenerator.h"
#include "TPZCompMeshTools.h"
#include "pzelchdivbound2.h"
#include "pzl2projection.h"
#include "pzshapequad.h"

#include "pzbndcond.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzcompelwithmem.h"

#include "pzcreateapproxspace.h"

class TRMSpaceOdissey{
    
public:
    
    /** @brief Define the type of geometry being used */
    enum MGeoMeshType {ENone = 0,EBox = 1, EReservoir = 2, EGIDReservoir = 3};
    
    /** @brief Define the type of geometry being used */
    MGeoMeshType fMeshType;
    
private:
    
    /** @brief order of approximation displacements */
    int fUOrder;
    
    /** @brief order of approximation flux and pressure */
    int fPOrder;
    
    /** @brief order of approximation transport */
    int fSOrder;
    
    /** @brief Autopointer of the Geometric mesh shared with all the classes involved */
    TPZGeoMesh * fGeoMesh;
    
    /** @brief Autopointer of Simulation data */
    TRMSimulationData * fSimulationData;
    
    /** @brief H1 computational mesh for biot linear poroelasticity approach */
    TPZCompMesh * fBiotCmesh;
    
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
    
    /** @brief Mixed computational mhm mesh for a dual analysis */
    TPZCompMesh * fMixedFluxPressureCmeshMHM;
    
    /** @brief H1-L2 for computational mesh for a primal analysis with Global postprocessing of fluxes */
    TPZCompMesh * fPressureSaturationCmesh;
    
    /** @brief Computational mesh for multiphase monolithic approach */
    TPZCompMesh * fMonolithicMultiphaseCmesh;
    
    void ModifyElementOrders(std::map<int64_t,int> &elorders);

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
    
    /** @brief Set the order of approximation displacements */
    void SetDefaultUOrder(int porder)
    {
        fUOrder = porder;
    }
    
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
    
    /** @brief Create a Biot H1 computational mesh */
    void CreateBiotCmesh();
    
    /** @brief Create a H1 computational mesh */
    void CreateH1Cmesh();
    
    /** @brief Create a Hdiv computational mesh */
    void CreateFluxCmesh();

    /** @brief Create a Discontinuous computational mesh L2 */
    void CreatePressureCmesh();
    
    /** @brief Create a Mixed computational mesh Hdiv-L2 */
    void CreateMixedCmesh();
    
    /** @brief Create a Mixed computational mesh Hdiv-L2 */
    void CreateMixedCmeshMHM();
    
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
    
    /** @brief Create a reservoir-box geometry with cylindrical wells */
    void CreateGeometricGmshMesh(std::string &grid);
    
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
    
    /** @brief Apply uniform refinement on the given father index mesh */
    void UniformRefinement_at_Father(int n_ref, int father_index);
    
    /** @brief Apply uniform refinement at specific material id */
    void UniformRefinement_at_MaterialId(int n_ref, int mat_id);
    
    /** @brief Apply uniform refinement around at specific material id */
    void UniformRefinement_Around_MaterialId(int n_ref, int mat_id);
    
    /** @brief Set autopointer of Simulation data */
    void SetSimulationData(TRMSimulationData * SimulationData){
        fSimulationData = SimulationData;
    }
    
    /** @brief Get autopointer of Simulation data */
    TRMSimulationData *  SimulationData(){
        return fSimulationData;
    }
    
    
    /** @brief Set autopointer of Biot H1 computational mesh for primal approach */
    void SetBiotCMesh(TPZCompMesh * BiotCmesh){
        fBiotCmesh = BiotCmesh;
    }
    
    /** @brief Get autopointer of Biot H1 computational mesh for primal approach */
    TPZCompMesh *  BiotCMesh(){
        return fBiotCmesh;
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
    
    /** @brief Set autopointer of Mixed computational mesh for a dual analysis */
    void SetMixedFluxPressureCmeshMHM(TPZCompMesh * MixedFluxPressureCmeshMHM){
        fMixedFluxPressureCmeshMHM = MixedFluxPressureCmeshMHM;
    }
    
    /** @brief Get autopointer of Mixed computational mesh for a dual analysis */
    TPZCompMesh * MixedFluxPressureCmeshMHM(){
        return fMixedFluxPressureCmeshMHM;
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
    
    /** @brief Build MHM form the current hdvi mesh */
    void BuildMixed_Mesh();
    
    /** @brief Create computational interfaces inside for flux computations  */
    void CreateInterfacesInside(TPZCompMesh * cmesh);
    
    /** @brief Adjust the polinomial order of the elements */
    void IncreaseOrderAroundWell(int numlayers);
    
    // MHM
    
    /** @brief Build MHM form the current hdvi mesh */
    void BuildMHM_Mesh();
    
    /** @brief Build MHM form the current hdvi mesh */
    void InsertSkeletonInterfaces(int skeleton_id = 0);
    
    /** @brief Sparated connects by given selected skeleton ids */
    void SeparateConnectsBySkeletonIds(TPZVec<int64_t> skeleton_ids);
    
    /** @brief Sparated connects by hdiv connect neighborhood */
    void SeparateConnectsByNeighborhood();
    
    /** @brief Construc computational macro elements */
    void BuildMacroElements();
    
    /** @brief Destruct computational macro elements */
    void UnwrapMacroElements();
    
    /** @brief Computational element refinement by index */
    void CElemtentRefinement(TPZCompMesh  *cmesh, int element_index);
    
    /** @brief Computational mesh uniform refinement */
    void CMeshRefinement(TPZCompMesh  *cmesh, int ndiv);
    
};

#endif /* defined(__PZ__TRMSpaceOdissey__) */
