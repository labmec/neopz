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

///class oriented towards the creation of Discontiunos Petrov-Galerkin (DPG) meshes
class TPZDPGMeshControl
{
    /// geometric mesh used to create the computational mesh
    TPZAutoPointer<TPZGeoMesh> fGMesh;
    
     /// computational mesh to contain the pressure elements in the coarse mesh
    TPZCompMesh fPressureCoarseMesh;
    
    ///pointers to MHM mesh
    TPZMHMeshControl fMHMControl;
    
    /// interpolation order of the internal elements at coarse mesh
    int fPOrderCoarseInternal;
    
    ///Materials Id of fine mesh
    int fFinerMatId;
    
    ///Materials Id of coarse mesh
    int fCoarseMatId;
    
    ///Materials Id of skeleton mesh. Elements with dimension equal to dim(mesh)-1
    int fSkeletMatId;
    
public:
    
    TPZDPGMeshControl() : fPOrderCoarseInternal(-1), fFinerMatId(0),fCoarseMatId(0), fSkeletMatId(0)
    {
        
    }
    
    TPZDPGMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<int64_t> &coarseindices);
    //~TPZDPGMeshControl();

    
    TPZDPGMeshControl(const TPZDPGMeshControl &copy);
    
    TPZDPGMeshControl &operator=(const TPZDPGMeshControl &cp);
    
    void SetPOrderMeshes(int porderfinermesh, int pordercoarsemesh, int porderskeleton)
    {
        fPOrderCoarseInternal = pordercoarsemesh;
        fPressureCoarseMesh.SetDefaultOrder(pordercoarsemesh);
        
        fMHMControl.SetInternalPOrder(porderfinermesh);
        fMHMControl.SetSkeletonPOrder(porderskeleton);
    }
    
    TPZMHMeshControl & MHMControl(){
        return fMHMControl;
    }
    
    TPZCompMesh & PressureCoarseMesh(){
        return fPressureCoarseMesh;
    }
    
    void BuildComputationalMesh();
    
    void SetMatIds(int finermat, int coarsemat, int skeletonmat){
   
        if((finermat == skeletonmat) || (coarsemat==skeletonmat)) DebugStop();
        fFinerMatId = finermat;
        fCoarseMatId = coarsemat;
        fSkeletMatId = skeletonmat;
    }
    
    /// Put the pointers to the meshes in a vector
    void GetMeshVec(TPZVec<TPZCompMesh *> &meshvec)
    {
        TPZVec<TPZCompMesh *> meshvecold;
        fMHMControl.GetMeshVec(meshvecold);
        meshvec.Resize(4);
        meshvec[0] = meshvecold[0];
        meshvec[1] = meshvecold[1];
        meshvec[2] = meshvecold[2];
        meshvec[3] = &fPressureCoarseMesh;
    }
};

#endif /* defined(__PZ__TPZDPGMeshControl__) */