//
//  TPZMHMixedMeshControl.hpp
//  PZ
//


#ifndef TPZMHMixedMeshChannelControl_hpp
#define TPZMHMixedMeshChannelControl_hpp

#include <stdio.h>

#include "TPZMHMeshControl.h"
#include "TPZMHMixedMeshControl.h"

/// class for creating TPZMHMM with Mixed Meshes
class TPZMHMixedMeshChannelControl : public TPZMHMixedMeshControl
{
    

public:
    
    TPZMHMixedMeshChannelControl() : TPZMHMixedMeshControl()
    {
        
    }
    
    TPZMHMixedMeshChannelControl(int dimension):TPZMHMixedMeshControl(dimension){
        
    }

    
    TPZMHMixedMeshChannelControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices):TPZMHMixedMeshControl( gmesh, coarseindices){
        
    }
    
    TPZMHMixedMeshChannelControl(TPZAutoPointer<TPZGeoMesh> gmesh):TPZMHMixedMeshControl(gmesh)
    {
        
    }
    
    
    TPZMHMixedMeshChannelControl(const TPZMHMixedMeshChannelControl &copy) : TPZMHMixedMeshControl(copy)
    {
       
         fFluxMesh = copy.fFluxMesh;
    }
    
    TPZMHMixedMeshChannelControl &operator=(const TPZMHMixedMeshChannelControl &cp)
    {
        fFluxMesh = cp.fFluxMesh;
        TPZMHMixedMeshControl::operator=(cp);
        return *this;
    }
    
    void BuildComputationalMesh(bool usersubstructure,bool OpenChannel,std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>>);
    
    void HideTheElements();
    
    int64_t WhichSubdomain(TPZCompEl *cel);
    
};

#endif /* TPZMHMixedMeshChannelControl_hpp */
