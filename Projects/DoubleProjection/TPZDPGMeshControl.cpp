//
//  TPZDPGMeshControl.cpp
//  PZ
//
//  Created by Agnaldo Farias on 26/08/14.
//
//

#include "TPZDPGMeshControl.h"


TPZDPGMeshControl::TPZDPGMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<long> &coarseindices) : fGMesh(gmesh), fPOrderCoarseInternal(-1)
{
 
    fPressureCoarseMesh = new TPZCompMesh(fGMesh);
    
    fMHMControl  = new TPZMHMeshControl(fGMesh,coarseindices);
}

TPZDPGMeshControl::TPZDPGMeshControl(const TPZDPGMeshControl &copy){
    
      this->operator=(copy);
}

TPZDPGMeshControl & TPZDPGMeshControl::operator=(const TPZDPGMeshControl &cp){
    
    fGMesh = cp.fGMesh;
    fMHMControl = cp.fMHMControl;
    fPOrderCoarseInternal = cp.fPOrderCoarseInternal;
}
