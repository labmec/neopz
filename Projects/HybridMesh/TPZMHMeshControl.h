//
//  TPZMHMeshControl.h
//  PZ
//
//  Created by Philippe Devloo on 5/3/14.
//
//

#ifndef __PZ__TPZMHMeshControl__
#define __PZ__TPZMHMeshControl__

#include <iostream>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"

/// class oriented towards the creation of multiscale hybrid meshes
class TPZMHMeshControl
{
    /// geometric mesh used to create the computational mesh
    TPZAutoPointer<TPZGeoMesh> fGMesh;
    
    /// computational MHM mesh being built by this class
    TPZAutoPointer<TPZCompMesh> fCMesh;
    
    /// material id associated with the skeleton elements
    int fSkeletonMatId;
    
    /// material id associated with the lagrange multiplier elements
    int fLagrangeMatIdLeft, fLagrangeMatIdRight;
    
    /// interpolation order of the internal elements
    int fpOrderInternal;
    
    /// indices of the geometric elements which define the skeleton mesh
    std::set<long> fCoarseIndices;
    
    /// indices of the skeleton elements and their left/right elements of the skeleton mesh
    std::map<long, std::pair<long,long> > fInterfaces;
    
public:
    TPZMHMeshControl()
    {
        
    }
    
    TPZMHMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<long> &coarseindices);
    
    TPZAutoPointer<TPZCompMesh> CMesh()
    {
        return fCMesh;
    }
    
    /// Create all data structures for the computational mesh
    void BuildComputationalMesh();
    
    /// will create 1D elements on the interfaces between the coarse element indices
    void CreateCoarseInterfaces(int matid);
    
    /// Print diagnostics
    void PrintDiagnostics(std::ostream &out);
    
private:
    /// will create a computational mesh using the coarse element indexes and its interface elements
    TPZCompMesh *CriaMalhaTemporaria();
    
    /// will create the internal elements, one coarse element at a time
    void CreateInternalElements();
    
    /// Add the boundary elements to the computational mesh
    void AddBoundaryElements();
    
    /// Add the boundary interface elements to the computational mesh
    void AddBoundaryInterfaceElements();
    
    /// will create the elements on the skeleton
    void CreateSkeleton();
    
    /// will create the interface elements between the internal elements and the skeleton
    void CreateInterfaceElements();
    
    /// verify if the element is a sibling of
    bool IsSibling(long son, long father);
    
    /// put the element side which face the boundary on the stack
    void AddElementBoundaries(long elseed, long compelindex, TPZStack<TPZCompElSide> &result);
    
    /// print the diagnostics for a subdomain
    void PrintSubdomain(long elindex, std::ostream &out);
    
    /// print the indices of the boundary elements and interfaces
    void PrintBoundaryInfo(std::ostream &out);
};

#endif /* defined(__PZ__TPZMHMeshControl__) */
