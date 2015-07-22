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

#include "pzgmesh.h"
#include "pzcmesh.h"

#include "pzbuildmultiphysicsmesh.h"
#include "pzcondensedcompel.h"

#include "TPZVTKGeoMesh.h"



class TRMSpaceOdissey{
    
public:
    
    /** @brief Define the type of geometry being used */
    enum MGeoMeshType {ENone = 0,EBox = 1, EReservoir = 2};
    
private:
    
    /** @brief order of approximation */
    int fPOrder;
    
    /** @brief Define the type of geometry being used */
    MGeoMeshType fMeshType;
    
    /** @brief Autopointer of the Geometric mesh shared with all the classes involved */
    TPZAutoPointer<TPZGeoMesh> fGeoMesh;
    
    /** @brief H1 computational mesh for validation */
    TPZAutoPointer<TPZCompMesh> fH1Cmesh;
    
    /** @brief Autopointer of Simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
    /** @brief Hdiv computational mesh conservative vector field */
    TPZAutoPointer<TPZCompMesh> fFluxCmesh;
    
    /** @brief L2 computational mesh the restriction equation */
    TPZAutoPointer<TPZCompMesh> fPressureCmesh;
    
    /** @brief Mixed computational mesh for a dual analysis */
    TPZAutoPointer<TPZCompMesh> fMixedFluxPressureCmesh;
    
    /** @brief H1 computational mesh for Maurice Biot Linear Poroelasticity */
    TPZAutoPointer<TPZCompMesh> fGeoMechanicsCmesh;
    
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
    
    /** @brief Statically condense the internal equations of the elements */
    void StaticallyCondenseEquations();
    
    /** @brief Configure the boundary conditions of a well with reservoir boundary conditions */
    void ConfigureWellConstantPressure(STATE wellpressure, STATE farfieldpressure);
    
    /** @brief Create a computational mesh L2 */
    void CreateTransportMesh();
    
    /** @brief Create a H1 computational mesh for Maurice Biot Linear Poroelasticity */
    void CreateGeoMechanicMesh();
    
    /** @brief Create the reservoir geometry */
    void CreateGeometricReservoirMesh();
    
    /** @brief Print the reservoir geometry */
    void PrintGeometry();
    
    /** @brief Create a reservoir-box geometry */
    void CreateGeometricBoxMesh(TPZManVector<int,2> dx, TPZManVector<int,2> dy, TPZManVector<int,2> dz);
    
    /** @brief Parametric function that computes elements in the x direction */
    static  void ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /** @brief Parametric function that computes elements in the y direction */
    static  void ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /** @brief Parametric function that computes elements in the z direction */
    static  void ParametricfunctionZ(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /**
     * @ingroup Acces methods
     * @brief Set and Get fucntion attributes
     * @since June 09, 2015
     */
    
    /** @brief Autopointer of the Geometric mesh shared with all the classes involved */
    void SetGmesh(TPZAutoPointer<TPZGeoMesh> GeoMesh){
        fGeoMesh = GeoMesh;
    }
    TPZAutoPointer<TPZGeoMesh>  GetGmesh(){
        return fGeoMesh;
    }
    
    /** @brief H1 computational mesh for validation */
    void SetH1Cmesh(TPZAutoPointer<TPZCompMesh> H1Cmesh){
        fH1Cmesh = H1Cmesh;
    }
    TPZAutoPointer<TPZCompMesh>  GetH1Cmesh(){
        return fH1Cmesh;
    }

    /// Access method for the flux mesh
    TPZAutoPointer<TPZCompMesh>  GetFluxCmesh(){
        return fFluxCmesh;
    }
    
    /// Access method for the pressure mesh
    TPZAutoPointer<TPZCompMesh>  GetPressureMesh(){
        return fPressureCmesh;
    }
    

    /** @brief Mixed computational mesh for validation */
    void SetMixedCmesh(TPZAutoPointer<TPZCompMesh> MixedCmesh){
        fMixedFluxPressureCmesh = MixedCmesh;
    }
    TPZAutoPointer<TPZCompMesh>  GetMixedCmesh(){
        return fMixedFluxPressureCmesh;
    }
    
    /// Adjust the polinomial order of the elements
    void IncreaseOrderAroundWell(int numlayers);
};

#endif /* defined(__PZ__TRMSpaceOdissey__) */
