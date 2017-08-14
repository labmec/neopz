//
//  TPZMHMixedHybridMeshControl.hpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#ifndef TPZMHMixedHybridMeshControl_hpp
#define TPZMHMixedHybridMeshControl_hpp

#include <stdio.h>

#include "TPZMHMixedMeshControl.h"
class TPZMultiphysicsElement;

/// class for creating TPZMHMM with Mixed Meshes
class TPZMHMixedHybridMeshControl : public TPZMHMixedMeshControl
{
    
    

public:
    
    TPZMHMixedHybridMeshControl() : TPZMHMixedMeshControl(), fHDivWrapperMatId(105)
    {
        
    }
    
    TPZMHMixedHybridMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<long> &coarseindices);
    
    TPZMHMixedHybridMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<long> &coarseindices);
    
    
    TPZMHMixedHybridMeshControl(const TPZMHMixedHybridMeshControl &copy) : TPZMHMixedMeshControl(copy), fHDivWrapperMatId(copy.fHDivWrapperMatId)
    {
    }
    
    virtual ~TPZMHMixedHybridMeshControl()
    {
        
    }
    
    TPZMHMixedHybridMeshControl &operator=(const TPZMHMixedHybridMeshControl &cp)
    {
        TPZMHMixedMeshControl::operator=(cp);
        fHDivWrapperMatId = cp.fHDivWrapperMatId;
        return *this;
    }
    
    /// Insert Boundary condition objects that do not perform any actual computation
    virtual void InsertPeriferalMaterialObjects();
    

protected:
    
    int fHDivWrapperMatId;
    
    // create the elements domain per domain with approximation spaces disconnected from each other
    virtual void CreateInternalElements();

    // create the approximation space associated with the skeleton and restrain the connects
    virtual void CreateSkeleton();
        
    /// Create the interfaces between the pressure elements of dimension dim
    virtual void CreateMultiPhysicsInterfaceElements(int dim);


protected:
    
    // create the dim-1 pressure elements between the hdiv elements
    void CreatePressureInterfaces();
    
    // create a boundary element around each flux element to optimize calculations
    void CreateHDivWrappers();
    
    // Find the two connected flux elements to a pressure element
    void FindConnectedElements(TPZGeoElSide pressureindex, TPZVec<TPZCompElSide> &fluxconnected);
    
    /// group and condense the elements
    virtual void GroupandCondenseElements();

    /// group element H(div) elements with surrounding interface elements
    void GroupElements(TPZCompMesh *cmesh);
    
};

#endif /* TPZMHMixedHybridMeshControl_hpp */
