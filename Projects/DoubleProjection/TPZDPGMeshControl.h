//
//  TPZDPGMeshControl.h
//  PZ
//
//  Created by Agnaldo Farias on 26/08/14.
//
//

#ifndef __PZ__TPZDPGMeshControl__
#define __PZ__TPZDPGMeshControl__

#include <iostream>

#include "pzcmesh.h"
#include "pzgmesh.h"
#include "TPZMHMeshControl.h"

///class  oriented towards the creation of Discontiunos Petrov-Galerkin (DPG) meshes
class TPZDPGMeshControl
{
    /// geometric mesh used to create the computational mesh
    TPZAutoPointer<TPZGeoMesh> fGMesh;
    
     /// computational mesh to contain the pressure elements in the coarse mesh
    TPZAutoPointer<TPZCompMesh> fPressureCoarseMesh;
    
    ///pointers to MHM mesh
    TPZMHMeshControl * fMHMControl;
    
    /// interpolation order of the internal elements at coarse mesh
    int fPOrderCoarseInternal;
    
public:
    
    TPZDPGMeshControl(){
    };
    
    TPZDPGMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<long> &coarseindices);
    
    ~TPZDPGMeshControl(){
    }
    
    TPZDPGMeshControl(const TPZDPGMeshControl &copy);
    
    TPZDPGMeshControl &operator=(const TPZDPGMeshControl &cp);
    
    void SetInternalCoarsePOrder(int order)
    {
        fPOrderCoarseInternal = order;
        fPressureCoarseMesh->SetDefaultOrder(order);
    }
    
    TPZMHMeshControl * MHMMControl(){
        return fMHMControl;
    }
};

#endif /* defined(__PZ__TPZDPGMeshControl__) */