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
    
    /// computational mesh to represent the constant states
    TPZAutoPointer<TPZCompMesh> fCMeshLagrange;
    
    /// computational mesh to represent the constant states
    TPZAutoPointer<TPZCompMesh> fCMeshConstantPressure;
    
    /// computational mesh to contain the pressure elements
    // this mesh is the same as fCMesh if there are no lagrange multipliers assocated with the average pressure
    TPZAutoPointer<TPZCompMesh> fPressureFineMesh;
    
    /// material id associated with the skeleton elements
    int fSkeletonMatId;
    
    /// material id associated with the lagrange multiplier elements
    int fLagrangeMatIdLeft, fLagrangeMatIdRight;
    
    /// interpolation order of the internal elements
    int fpOrderInternal;
    
    /// interpolation order of the skeleton elements
    int fpOrderSkeleton;
    
    /// indices of the geometric elements which define the skeleton mesh
    std::set<long> fCoarseIndices;
    
    /// indices of the skeleton elements and their left/right elements of the skeleton mesh
    std::map<long, std::pair<long,long> > fInterfaces;
    
    /// flag to determine whether a lagrange multiplier is included to force zero average pressures in the subdomains
    /**
     * when imposing average pressure to be zero, a multiphysics mesh is created
     */
    bool fLagrangeAveragePressure;
    
public:
    TPZMHMeshControl() : fSkeletonMatId(-1), fLagrangeMatIdLeft(-1), fLagrangeMatIdRight(-1), fpOrderInternal(-1), fpOrderSkeleton(-1), fLagrangeAveragePressure(false)
    {
        
    }
    
    TPZMHMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<long> &coarseindices);
    
    TPZMHMeshControl(const TPZMHMeshControl &copy);
    
    TPZMHMeshControl &operator=(const TPZMHMeshControl &cp);
    
    TPZAutoPointer<TPZCompMesh> CMesh() const
    {
        return fCMesh;
    }
    
    TPZAutoPointer<TPZGeoMesh> GMesh() const
    {
        return fGMesh;
    }
    
    /// Set the porder for the internal elements
    void SetInternalPOrder(int order)
    {
        fpOrderInternal = order;
        fCMesh->SetDefaultOrder(order);
    }
    
    void SetSkeletonPOrder(int order)
    {
        fpOrderSkeleton = order;
    }
    
    /// Set the flag for creating Lagrange Dofs for the average pressure
    void SetLagrangeAveragePressure(bool flag)
    {
        fLagrangeAveragePressure = flag;
    }
    
    /// Create all data structures for the computational mesh
    void BuildComputationalMesh(bool usersubstructure);
    
    /// will create 1D elements on the interfaces between the coarse element indices
    void CreateCoarseInterfaces(int matid);
    
    /// Print diagnostics
    void PrintDiagnostics(std::ostream &out);
    
    /// Put the pointers to the meshes in a vector
    void GetMeshVec(TPZVec<TPZCompMesh *> &meshvec)
    {
        meshvec.Resize(3);
        meshvec[0] = fPressureFineMesh.operator->();
        meshvec[1] = fCMeshLagrange.operator->();
        meshvec[2] = fCMeshConstantPressure.operator->();
    }
    

    
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
    
    /// create the lagrange multiplier mesh, one element for each subdomain
    void CreateLagrangeMultiplierMesh();
    
    /// transform the computational mesh into a multiphysics mesh
    void TransferToMultiphysics();
    
    /// substructure the mesh
    void SubStructure();
    
    
    /// print the diagnostics for a subdomain
    void PrintSubdomain(long elindex, std::ostream &out);
    
    /// print the indices of the boundary elements and interfaces
    void PrintBoundaryInfo(std::ostream &out);
};

#endif /* defined(__PZ__TPZMHMeshControl__) */
