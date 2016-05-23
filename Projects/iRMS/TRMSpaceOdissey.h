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

#include "TPZVTKGeoMesh.h"



class TRMSpaceOdissey{
    
public:
    
    /** @brief Define the type of geometry being used */
    enum MGeoMeshType {ENone = 0,EBox = 1, EReservoir = 2, EGIDReservoir = 3};
    
    /** @brief Define the type of geometry being used */
    MGeoMeshType fMeshType;
    
private:
    
    /** @brief order of approximation */
    int fPOrder;
    
    /** @brief Autopointer of the Geometric mesh shared with all the classes involved */
    TPZAutoPointer<TPZGeoMesh> fGeoMesh;
    
    /** @brief Autopointer of Simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
    /** @brief H1 computational mesh for primal approach */
    TPZAutoPointer<TPZCompMesh> fH1Cmesh;
    
    /** @brief Hdiv computational mesh conservative vector field */
    TPZAutoPointer<TPZCompMesh> fFluxCmesh;
    
    /** @brief L2 computational mesh the restriction equation */
    TPZAutoPointer<TPZCompMesh> fPressureCmesh;
    
    /** @brief L2 computational mesh water saturation equations */
    TPZAutoPointer<TPZCompMesh> fWaterSaturationMesh;

    /** @brief L2 computational mesh oil saturation equations */
    TPZAutoPointer<TPZCompMesh> fOilSaturationMesh;
    
    /** @brief H1 computational mesh for Maurice Biot linear poroelasticity */
    TPZAutoPointer<TPZCompMesh> fGeoMechanicsCmesh;
    
    /** @brief Mixed computational mesh for a dual analysis */
    TPZAutoPointer<TPZCompMesh> fMixedFluxPressureCmesh;
    
    /** @brief H1-L2 for computational mesh fora primal analysis with Global postprocessing of fluxes */
    TPZAutoPointer<TPZCompMesh> fPressureSaturationCmesh;
    
    /** @brief Computational mesh for multiphase monolithic approach */
    TPZAutoPointer<TPZCompMesh> fMonolithicMultiphaseCmesh;
    
    /** @brief Object that generates and contains all the sparse matrices to transfer informations */
    TPZAutoPointer<TRMBuildTransfers> fTransferGenerator;
    
    void ModifyElementOrders(std::map<long,int> &elorders);

public:
    
    /** @brief Default constructor */
    TRMSpaceOdissey();
    
    /** @brief Default desconstructor */
    ~TRMSpaceOdissey();
    
    /** @brief Initialize the simulation data */
    void InitializeSimulationData(TRMRawData &rawdata);
    
    /// Change the default polynomial order
    void SetDefaultPOrder(int porder)
    {
        fPOrder = porder;
    }
    
    /** @brief Create a H1 computational mesh */
    void CreateH1Cmesh();
    
    /** @brief Create a Hdiv computational mesh Hdiv */
    void CreateFluxCmesh();

    /** @brief Create a Discontinuous computational mesh L2 */
    void CreatePressureCmesh();
    
    /** @brief Create a Mixed computational mesh Hdiv-L2 */
    void CreateMixedCmesh();
    
    /** @brief Create a Mixed-Transport muliphase computational mesh Hdiv-L2-L2-L2 */
    void CreateMultiphaseCmesh();
    
    /** @brief Statically condense the internal equations of the elements */
    void StaticallyCondenseEquations();
    
    /** @brief Configure the boundary conditions of a well with reservoir boundary conditions */
    void ConfigureWellConstantPressure(STATE wellpressure, STATE farfieldpressure);
    
    /** @brief Create a computational mesh L2 */
    void CreateWaterTransportMesh();
    
    /** @brief Create a computational mesh L2 */
    void CreateOilTransportMesh();
    
    /** @brief Create a H1 computational mesh for Maurice Biot Linear Poroelasticity */
    void CreateGeoMechanicMesh();
    
    /** @brief Create the reservoir geometry */
    void CreateGeometricReservoirMesh();
    
    /** @brief Print the reservoir geometry */
    void PrintGeometry();
    
    /** @brief Create a reservoir-box geometry */
    void CreateGeometricBoxMesh(TPZManVector<REAL,2> dx, TPZManVector<REAL,2> dy, TPZManVector<REAL,2> dz);
    
    /** @brief Parametric function that computes elements in the x direction */
    static  void ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /** @brief Parametric function that computes elements in the y direction */
    static  void ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /** @brief Parametric function that computes elements in the z direction */
    static  void ParametricfunctionZ(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    
    
    /**
     * @defgroup Access Methods
     * @brief    Implements Access methods:
     * @{
     */
    
    /** @brief Set autopointer of the Geometric mesh */
    void SetGmesh(TPZAutoPointer<TPZGeoMesh> &GeoMesh){
        fGeoMesh = GeoMesh;
    }
    
    /** @brief Get autopointer of the Geometric mesh */
    TPZAutoPointer<TPZGeoMesh>  &Gmesh(){
        return fGeoMesh;
    }
    
    /** @brief Set autopointer of Simulation data */
    void SetSimulationData(TPZAutoPointer<TRMSimulationData> &SimulationData){
        fSimulationData = SimulationData;
    }
    
    /** @brief Get autopointer of Simulation data */
    TPZAutoPointer<TRMSimulationData>  SimulationData(){
        return fSimulationData;
    }
    
    /** @brief Set autopointer of H1 computational mesh for primal approach */
    void SetH1CMesh(TPZAutoPointer<TPZCompMesh> &H1Cmesh){
        fH1Cmesh = H1Cmesh;
    }
    
    /** @brief Get autopointer of H1 computational mesh for primal approach */
    TPZAutoPointer<TPZCompMesh>  H1CMesh(){
        return fH1Cmesh;
    }
    
    /** @brief Set autopointer of Hdiv computational mesh conservative vector field */
    void SetFluxCmesh(TPZAutoPointer<TPZCompMesh> &FluxCmesh){
        fFluxCmesh = FluxCmesh;
    }
    
    /** @brief Get autopointer of Hdiv computational mesh conservative vector field */
    TPZAutoPointer<TPZCompMesh>  FluxCmesh(){
        return fFluxCmesh;
    }
    
    /** @brief Set autopointer of L2 computational mesh the restriction equation */
    void SetPressureCmesh(TPZAutoPointer<TPZCompMesh> &PressureCmesh){
        fPressureCmesh = PressureCmesh;
    }
    
    /** @brief Get autopointer of L2 computational mesh the restriction equation */
    TPZAutoPointer<TPZCompMesh>  PressureCmesh(){
        return fPressureCmesh;
    }
    
    /** @brief Set autopointer of L2 computational mesh water saturation equations */
    void SetWaterSaturationMesh(TPZAutoPointer<TPZCompMesh> &WaterSaturationMesh){
        fWaterSaturationMesh = WaterSaturationMesh;
    }
    
    /** @brief Get autopointer of L2 computational mesh water saturation equations */
    TPZAutoPointer<TPZCompMesh>  WaterSaturationMesh(){
        return fWaterSaturationMesh;
    }
    
    /** @brief Set autopointer of L2 computational mesh oil saturation equations */
    void SetOilSaturationMesh(TPZAutoPointer<TPZCompMesh> &OilSaturationMesh){
        fOilSaturationMesh = OilSaturationMesh;
    }
    
    /** @brief Get autopointer of L2 computational mesh oil saturation equations */
    TPZAutoPointer<TPZCompMesh>  OilSaturationMesh(){
        return fOilSaturationMesh;
    }
    
    /** @brief Set autopointer of H1 computational mesh for Maurice Biot linear poroelasticity */
    void SetGeoMechanicsCmesh(TPZAutoPointer<TPZCompMesh> &GeoMechanicsCmesh){
        fGeoMechanicsCmesh = GeoMechanicsCmesh;
    }
    
    /** @brief Get autopointer of H1 computational mesh for Maurice Biot linear poroelasticity */
    TPZAutoPointer<TPZCompMesh>  GeoMechanicsCmesh(){
        return fGeoMechanicsCmesh;
    }
    
    /** @brief Set autopointer of Mixed computational mesh for a dual analysis */
    void SetMixedFluxPressureCmesh(TPZAutoPointer<TPZCompMesh> &MixedFluxPressureCmesh){
        fMixedFluxPressureCmesh = MixedFluxPressureCmesh;
    }
    
    /** @brief Get autopointer of Mixed computational mesh for a dual analysis */
    TPZAutoPointer<TPZCompMesh>  MixedFluxPressureCmesh(){
        return fMixedFluxPressureCmesh;
    }
    
    /** @brief Set autopointer of H1-L2 for computational mesh fora primal analysis with Global postprocessing of fluxes */
    void SetPressureSaturationCmesh(TPZAutoPointer<TPZCompMesh> &PressureSaturationCmesh){
        fPressureSaturationCmesh = PressureSaturationCmesh;
    }
    
    /** @brief Get autopointer of H1-L2 for computational mesh fora primal analysis with Global postprocessing of fluxes */
    TPZAutoPointer<TPZCompMesh>  PressureSaturationCmesh(){
        return fPressureSaturationCmesh;
    }
    
    /** @brief Set autopointer of Computational mesh for multiphase monolithic approach */
    void SetMonolithicMultiphaseCmesh(TPZAutoPointer<TPZCompMesh> &MonolithicMultiphaseCmesh){
        fMonolithicMultiphaseCmesh = MonolithicMultiphaseCmesh;
    }
    
    /** @brief Get autopointer of Computational mesh for multiphase monolithic approach */
    TPZAutoPointer<TPZCompMesh>  MonolithicMultiphaseCmesh(){
        return fMonolithicMultiphaseCmesh;
    }
    
    /** @brief Set autopointer of Object that generates and contains all the sparse matrices to transfer informations */
    void SetTransferGenerator(TPZAutoPointer<TRMBuildTransfers> &TransferGenerator){
        fTransferGenerator = TransferGenerator;
    }
    
    /** @brief Get autopointer of Object that generates and contains all the sparse matrices to transfer informations */
    TPZAutoPointer<TRMBuildTransfers>  TransferGenerator(){
        return fTransferGenerator;
    }
    
    
    // @}
    
    /**
     * @defgroup Utilities and other Methods
     * @brief    Implements Methods for h and p refinement:
     * @{
     */
    
    /** @brief Adjust the polinomial order of the elements */
    void IncreaseOrderAroundWell(int numlayers);
    
    // @}
    
    
};

#endif /* defined(__PZ__TRMSpaceOdissey__) */
