//
//  TPZMHMixedMeshControl.hpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#ifndef TPZMHMixedMeshControl_hpp
#define TPZMHMixedMeshControl_hpp

#include <stdio.h>

#include "TPZMHMeshControl.h"

/// class for creating TPZMHMM with Mixed Meshes
class TPZMHMixedMeshControl : public TPZMHMeshControl
{
    
protected:
    /// computational mesh to contain the pressure elements
    // this mesh is the same as fCMesh if there are no lagrange multipliers assocated with the average pressure
    TPZAutoPointer<TPZCompMesh> fFluxMesh;
    
    /// computational mesh to contain the rotation elements
    TPZAutoPointer<TPZCompMesh> fRotationMesh;
    

public:
    
    TPZMHMixedMeshControl() : TPZMHMeshControl()
    {
        
    }
    
    TPZMHMixedMeshControl(int dimension);
    //    TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<int64_t> &coarseindices);
    
    TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices);
    
    TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh);
    
    
    TPZMHMixedMeshControl(const TPZMHMixedMeshControl &copy) : TPZMHMeshControl(copy)
    {
        fFluxMesh = copy.fFluxMesh;
    }
    
    TPZMHMixedMeshControl &operator=(const TPZMHMixedMeshControl &cp)
    {
        fFluxMesh = cp.fFluxMesh;
        TPZMHMeshControl::operator=(cp);
        return *this;
    }
    
    virtual ~TPZMHMixedMeshControl();
    /// Insert Boundary condition objects that do not perform any actual computation
    virtual void InsertPeriferalMaterialObjects() override;
    
    /// Insert the necessary H(div) material objects to create the flux mesh
    virtual void InsertPeriferalHdivMaterialObjects();
    
    /// Insert the necessary Pressure material objects to create the flux mesh
    virtual void InsertPeriferalPressureMaterialObjects();
    
    /// Insert the necessary Rotation material objects to create the flux mesh
    virtual void InsertPeriferalRotationMaterialObjects();
    
    /// Create all data structures for the computational mesh
    virtual void BuildComputationalMesh(bool usersubstructure) override;
    
    TPZAutoPointer<TPZCompMesh> FluxMesh()
    {
        return fFluxMesh;
    }
    
    /// Set the flag for creating Lagrange Dofs for the average pressure
    void SetLagrangeAveragePressure(bool flag)
    {
        if (flag == true) {
            DebugStop();
        }
    }
    
    /// Set the hybridization to true
    virtual void Hybridize(bool flag)
    {
        DebugStop();
    }

    /// Put the pointers to the meshes in a vector
    void GetMeshVec(TPZVec<TPZCompMesh *> &meshvec)
    {
        meshvec.Resize(2);
        meshvec[0] = fFluxMesh.operator->();
        meshvec[1] = fPressureFineMesh.operator->();
    }

    TPZVec<TPZAutoPointer<TPZCompMesh> > GetMeshes()
    {
        TPZManVector<TPZAutoPointer<TPZCompMesh>,3> result(2);
        result[0] = fFluxMesh;
        result[1] = fPressureFineMesh;
        return result;
    }
    

    /// print the data structure
    void Print(std::ostream &out);

    /// print in a user friendly manner
    virtual void PrintFriendly(std::ostream &out)
    {
        
    }
    
protected:
    
   

    /// Create the mesh of the flux approximation space
    void CreateHDivMHMMesh();
    
    /// Create the pressure mesh which is dual to the flux mesh
    virtual void CreatePressureMHMMesh();
    
    /// Create the rotation mesh to elasticity problem
    virtual void CreateRotationMesh();
    
    // create the elements domain per domain with approximation spaces disconnected from each other
    virtual void CreateInternalFluxElements();
    
    // create the approximation space associated with the skeleton and restrain the connects
    virtual void CreateSkeleton();
    

    /// Create the multiphysics mesh
    void CreateHDivPressureMHMMesh();
    
    /// build the multi physics mesh (not at the finest geometric mesh level
    void BuildMultiPhysicsMesh();

    /// Create the interfaces between the pressure elements of dimension dim
    virtual void CreateMultiPhysicsInterfaceElements(int dim);
    
    /// Create the multiphysics interface elements between elements of specified material id
    virtual void CreateMultiPhysicsInterfaceElements(int dim, int pressmatid, std::pair<int,int> skelmatid);
    
    /// put the elements in TPZSubCompMesh, group the elements and condense locally
    void HideTheElements();

    /// hybridize the flux elements with the given material id - each flux element creates
    /// a pressure element
    virtual void HybridizeSkeleton(int skeletonmatid, int pressurematid) override;
    
    /// switch the elements pointed to by the interface by lower dimensional elements
    void OptimizeInterfaceElements();
    
    /// group and condense the elements
    virtual void GroupandCondenseElements();
    
    /// delete the pressure elements leaving the geometric mesh without pointing to the computational mesh
    void DeletePressureElements();

    void AdjustBoundaryElements();
};

#endif /* TPZMHMixedMeshControl_hpp */
