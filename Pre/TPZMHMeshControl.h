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

/// class oriented towards the creation of multiscale hybrid meshes - YES
class TPZMHMeshControl
{
    
protected:
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
    
    /// material id associated with the skeleton elements
    int fSecondSkeletonMatId;
    
    /// material id associated with the skeleton elements in a hybrid context
    int fPressureSkeletonMatId;
    
    /// material id associated with the lagrange multiplier elements
    int fLagrangeMatIdLeft, fLagrangeMatIdRight;
    
    /// interpolation order of the internal elements
    int fpOrderInternal;
    
    /// interpolation order of the skeleton elements
    int fpOrderSkeleton;
    
    /// indices of the geometric elements which define the skeleton mesh and their corresponding subcmesh indices
    std::map<long,long> fCoarseIndices;
    
    /// indices of the skeleton elements and their left/right geometric elements of the skeleton mesh
    std::map<long, std::pair<long,long> > fInterfaces;
    
    /// geometric index of the connects - subdomain where the connect will be internal
    TPZManVector<long> fConnectToSubDomainIdentifier;
    
    /// flag to determine whether a lagrange multiplier is included to force zero average pressures in the subdomains
    /**
     * when imposing average pressure to be zero, a multiphysics mesh is created
     */
    bool fLagrangeAveragePressure;
    
    /// flag to indicate whether we create a hybridized mesh
    bool fHybridize;
    
public:
    TPZMHMeshControl() : fSkeletonMatId(-1), fLagrangeMatIdLeft(-1), fLagrangeMatIdRight(-1), fpOrderInternal(-1), fpOrderSkeleton(-1), fLagrangeAveragePressure(false),
    fHybridize(false)
    {
        
    }
    
    TPZMHMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<long> &coarseindices);
    
    TPZMHMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<long> &coarseindices);
    
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
        if(fCMesh) fCMesh->SetDefaultOrder(order);
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
    
    /// Set the hybridization to true
    void Hybridize(int SecondSkeletonMatid, int PressureMatid)
    {
        fHybridize = true;
        // the three material ids must be different
        if (SecondSkeletonMatid == fSkeletonMatId || PressureMatid == fSkeletonMatId || SecondSkeletonMatid == PressureMatid) {
            DebugStop();
        }
        fSecondSkeletonMatId = SecondSkeletonMatid;
        fPressureSkeletonMatId = PressureMatid;
    }
    /// Create all data structures for the computational mesh
    void BuildComputationalMesh(bool usersubstructure);
    
    /// will create dim-1 geometric elements on the interfaces between the coarse element indices
    void CreateSkeletonElements(int skeletonmatid);
    
    /// divide the skeleton elements
    void DivideSkeletonElements(int ndivide);
    
    /// divide one skeleton element
    void DivideOneSkeletonElement(long index);
    
    /// print the data structure
    void Print(std::ostream &out);
    
    /// Print diagnostics
    void PrintDiagnostics(std::ostream &out);
    
    /// Put the pointers to the meshes in a vector
    void GetMeshVec(TPZVec<TPZCompMesh *> &meshvec)
    {
        if (fCMeshLagrange)
        {
            meshvec.Resize(3);
            meshvec[0] = fPressureFineMesh.operator->();
            meshvec[1] = fCMeshLagrange.operator->();
            meshvec[2] = fCMeshConstantPressure.operator->();
        }
        else
        {
            meshvec.resize(0);
        }
    }
    
    TPZVec<TPZAutoPointer<TPZCompMesh> > GetMeshes()
    {
        TPZManVector<TPZAutoPointer<TPZCompMesh>,3> result;
        if (fCMeshLagrange)
        {
            result.Resize(3);
            result[0] = fPressureFineMesh;
            result[1] = fCMeshLagrange;
            result[2] = fCMeshConstantPressure;
        }
        else
        {
            result.resize(0);
        }
        return result;
    }
    
    /// return the coarseindex to submesh index data structure
    std::map<long,long> &Coarse_to_Submesh()
    {
        return fCoarseIndices;
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
    
    /// hybridize the flux elements - each flux element becomes 5 elements
    void Hybridize();
    
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
    void SubStructure2();
    
protected:
    /// associates the connects of an element with a subdomain
    void SetSubdomain(TPZCompEl *cel, long subdomain, long offset = 0);
    
    /// associates the connects index with a subdomain
    void SetSubdomain(long connectindex, long subdomain, long offset = 0);
    
    /// returns to which subdomain a given element belongs
    // this method calls debugstop if the element belongs to two subdomains
    long WhichSubdomain(TPZCompEl *cel, long offset = 0);
    
    /// print the diagnostics for a subdomain
    void PrintSubdomain(long elindex, std::ostream &out);
    
    /// print the indices of the boundary elements and interfaces
    void PrintBoundaryInfo(std::ostream &out);
    
    /// identify connected elements to the skeleton elements
    // the computational mesh is determined by the element pointed to by the geometric element
    void ConnectedElements(long skeleton, std::pair<long,long> &leftright, std::map<long, std::list<TPZCompElSide> > &ellist);
};

#endif /* defined(__PZ__TPZMHMeshControl__) */
