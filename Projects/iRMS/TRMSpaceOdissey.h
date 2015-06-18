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

#include "TPZVTKGeoMesh.h"

class TRMSpaceOdissey{
    
public:
    
    /** @brief Define the type of geometry being used */
    enum MGeoMeshType {ENone = 0,EBox = 1, EReservoir = 2};
    
private:
    
    /** @brief Define the type of geometry being used */
    MGeoMeshType fMeshType;
    
    /** @brief Autopointer of the Geometric mesh shared with all the classes involved */
    TPZAutoPointer<TPZGeoMesh> fGeoMesh;
    
    /** @brief H1 computational mesh for validation */
    TPZAutoPointer<TPZCompMesh> fH1Mesh;
    
    /** @brief Autopointer of Simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
    /** @brief Hdiv computational mesh conservative vector field */
    TPZAutoPointer<TPZCompMesh> fFluxMesh;
    
    /** @brief L2 computational mesh the restriction equation */
    TPZAutoPointer<TPZCompMesh> fPressureMesh;
    
    /** @brief H1 computational mesh for Maurice Biot Linear Poroelasticity */
    TPZAutoPointer<TPZCompMesh> fGeoMechanicsMesh;
    

public:
    
    /** @brief Default constructor */
    TRMSpaceOdissey();
    
    /** @brief Initialize the simulation data */
    void InitializeSimulationData(TRMRawData &rawdata);
    
    /** @brief Create a H1 computational mesh */
    void CreateH1Mesh();
    
    /** @brief Create a Mixed computational mesh Hdiv-L2 */
    void CreateMixedMesh();
    
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
    
    
};

#endif /* defined(__PZ__TRMSpaceOdissey__) */
