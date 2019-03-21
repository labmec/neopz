//
//  TPZMultiphysicsMesh.hpp
//  pz
//
//  Created by Omar Durán on 3/20/19.
//

#ifndef TPZMultiphysicsMesh_h
#define TPZMultiphysicsMesh_h

#include <stdio.h>
#include "pzvec.h"
#include "pzcmesh.h"

class TPZMultiphysicsMesh : public TPZCompMesh {

private:
    
    TPZVector<int> m_active_meshes;
    
    TPZVector<int> m_inert_meshes;

public:
    
};

#endif /* TPZMultiphysicsMesh_h */
