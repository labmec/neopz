//
//  TPZMultiphysicsCompMesh.h
//  pz
//
//  Created by Omar Durán on 3/21/19.
//

#ifndef TPZMultiphysicsCompMesh_h
#define TPZMultiphysicsCompMesh_h

#include <stdio.h>
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include "TPZMultiphysicsInterfaceEl.h"

class TPZMultiphysicsCompMesh : public TPZCompMesh {
    
    /// Vector of active physics
    TPZManVector<int,5> m_active_approx_spaces;
    
    /// Vector of computational meshes
    TPZVec<TPZCompMesh * > m_mesh_vector;
    
public:
    
    /// Default constructor
    TPZMultiphysicsCompMesh();
    
    /// Constructor based on TPZGeoMesh pointer and vector of meshes
    TPZMultiphysicsCompMesh(TPZGeoMesh * gmesh);
    
    /// Copy constructor
    TPZMultiphysicsCompMesh(const TPZMultiphysicsCompMesh &other);
    
    /// Assignement constructor
    TPZMultiphysicsCompMesh & operator=(const TPZMultiphysicsCompMesh &other);
    
    /// Automatic builder for the computational mesh structure
    void AutoBuild();
    
    /// Set active approximation spaces
    void SetActiveApproxSpaces(TPZManVector<int,5> & active_approx_spaces, TPZVec<TPZCompMesh * > & mesh_vector);
    
    void LoadSolutionFromMeshes();
    
    void LoadSolutionFromMultiPhysics();
    
private:
    
    void AddElements();
    
    void AddConnects();
};

#endif /* TPZMultiphysicsCompMesh_h */
