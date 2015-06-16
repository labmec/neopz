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

#include "pzgmesh.h"
#include "pzcmesh.h"

/// Create the computational meshes
class TRMSpaceOdissey{
    
public:
    
    enum MGeoMeshType {ENone = 0,EBox = 1, EReservoir = 2};
    
private:
    /// Type of geometric mesh to be generated
    MGeoMeshType fMeshType;
    
    /// Geometric mesh shared by everybody
    TPZAutoPointer<TPZGeoMesh> fGeoMesh;
    
    /// H1 Mesh for initial validation
    TPZAutoPointer<TPZCompMesh> fH1Mesh;
    
    /// SimulationData, unique object for the execution
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
    /// HDiv approximation
    TPZAutoPointer<TPZCompMesh> fHDivMesh;
    
    /// Pressure mesh
    TPZAutoPointer<TPZCompMesh> fPressureMesh;
    
    /// Elastic deformation mesh
    TPZAutoPointer<TPZCompMesh> fGeoMechanicsMesh;
    
    /// Compound mesh
    TPZAutoPointer<TPZCompMesh> fFluxAndPressureMesh;

public:
    
    /// Default constructor
    TRMSpaceOdissey() : fMeshType(EBox)
    {
        
    }
    
    /// Initialize the TRMSimulationData
    void InitializeSimulationData(TRMRawData &rawdata);
    
    /// Create a H1 approximation mesh
    void CreateH1Mesh();
    
    /// Create a flux and pressure multiphysics mesh
    void CreateFluxPressureMesh();
    
    /// Create a mesh for saturation transport
    void CreateTransportMesh();
    
    /// Create a computational mesh for elastic deformation
    void CreateGeoMechanicMesh();
    
    /// Create a geometric mesh of a reservoir simulation
    void CreateGeometricReservoirMesh();
    
    /// Create a box mesh for initial testing
    void CreateGeometricBoxMesh();
    
    
};

#endif /* defined(__PZ__TRMSpaceOdissey__) */
