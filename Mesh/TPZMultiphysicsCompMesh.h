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
    
    /// Vector of active physics: index vector
    /// Define wich space will be active in order to generate equations. Should be defined for each space that you want to use by 0: no active or 1: active
    ///The size have to be the same as the m_mesh_vector
    TPZManVector<int,5> m_active_approx_spaces;
    
    /// Vector of computational meshes
    TPZManVector<TPZCompMesh * , 3> m_mesh_vector;
    
public:
    
    /// Default constructor
    TPZMultiphysicsCompMesh();
    
    /// Constructor based on TPZGeoMesh pointer and vector of meshes
    TPZMultiphysicsCompMesh(TPZGeoMesh * gmesh);
    
    /// Copy constructor
    TPZMultiphysicsCompMesh(const TPZMultiphysicsCompMesh &other);
    
    /// Assignement constructor
    TPZMultiphysicsCompMesh & operator=(const TPZMultiphysicsCompMesh &other);
    
    /// Destructor
    ~TPZMultiphysicsCompMesh();
    
    /// Automatic builder for the computational mesh structure
    void AutoBuild();
    
    /// Set active approximation spaces
    void BuildMultiphysicsSpace(TPZVec<int> & active_approx_spaces, TPZVec<TPZCompMesh * > & mesh_vector);
    
    void LoadSolutionFromMeshes();
    
    void LoadSolutionFromMultiPhysics();
    
    /// Get the vector of computational meshes
    TPZVec<TPZCompMesh *> & MeshVector();
    
    /// Get the vector of active physics
    TPZVec<int> & GetActiveApproximationSpaces();
    
private:
    
    void AddElements();
    
    void AddConnects();
};

#endif /* TPZMultiphysicsCompMesh_h */
