//
//  TPZMultiphysicsCompMesh.cpp
//  pz
//
//  Created by Omar Durán on 3/21/19.
//

#include "TPZMultiphysicsCompMesh.h"


TPZMultiphysicsCompMesh::TPZMultiphysicsCompMesh(){
    
    m_active_approx_spaces.Resize(0);
    m_mesh_vector.Resize(0);
    
}

TPZMultiphysicsCompMesh::TPZMultiphysicsCompMesh(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > & mesh_vector) : TPZCompMesh(gmesh){
    m_mesh_vector = mesh_vector;
}


void TPZMultiphysicsCompMesh::AutoBuild(){
    
    TPZCompMesh::AutoBuild();
    
    if (m_mesh_vector.size() != m_active_approx_spaces.size()) {
        DebugStop();
    }
    
    TPZBuildMultiphysicsMesh::AddElements(m_mesh_vector, this);
    TPZBuildMultiphysicsMesh::AddConnects(m_mesh_vector, this);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(m_mesh_vector, this);
}

void TPZMultiphysicsCompMesh::SetActiveApproxSpaces(TPZVec<int> & active_approx_spaces){
    m_active_approx_spaces = active_approx_spaces;
}

TPZMultiphysicsCompMesh::TPZMultiphysicsCompMesh(const TPZMultiphysicsCompMesh &other) : TPZCompMesh(other) {
    
}

TPZMultiphysicsCompMesh & TPZMultiphysicsCompMesh::operator=(const TPZMultiphysicsCompMesh &other){
    
    if (this != & other) // prevent self-assignment
    {
        TPZCompMesh::operator=(other);
        m_active_approx_spaces  = other.m_active_approx_spaces;
        m_mesh_vector           = other.m_mesh_vector;
    }
    return *this;
}
