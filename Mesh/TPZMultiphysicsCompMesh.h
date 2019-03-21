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
#include "pzbuildmultiphysicsmesh.h"

class TPZMultiphysicsCompMesh : public TPZCompMesh {
    
    /// Vector of active physics
    TPZVec<int> m_active_approx_spaces;
    
    /// Vector of computational meshes
    TPZVec<TPZCompMesh * > m_mesh_vector;
    
public:
    
    /// Default constructor
    TPZMultiphysicsCompMesh();
    
    /// Constructor based on TPZGeoMesh pointer and vector of meshes
    TPZMultiphysicsCompMesh(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > & mesh_vector);
    
    /// Copy constructor
    TPZMultiphysicsCompMesh(const TPZMultiphysicsCompMesh &other);
    
    /// Assignement constructor
    TPZMultiphysicsCompMesh & operator=(const TPZMultiphysicsCompMesh &other);
    
    /// Automatic builder for the computational mesh structure
    void AutoBuild();
    
    /// Set active approximation spaces
    void SetActiveApproxSpaces(TPZVec<int> & active_approx_spaces);
};

#endif /* TPZMultiphysicsCompMesh_h */
