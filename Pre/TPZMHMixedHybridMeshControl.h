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
    
    TPZMHMixedHybridMeshControl() : TPZMHMixedMeshControl()
    {
        
    }
    

    //    TPZMHMixedHybridMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<int64_t> &coarseindices);
    
    TPZMHMixedHybridMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices);
    
    TPZMHMixedHybridMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh);
    
    
    TPZMHMixedHybridMeshControl(const TPZMHMixedHybridMeshControl &copy) : TPZMHMixedMeshControl(copy), fHDivWrapperMatId(copy.fHDivWrapperMatId)
    {
    }
    
    virtual ~TPZMHMixedHybridMeshControl();
    
    TPZMHMixedHybridMeshControl &operator=(const TPZMHMixedHybridMeshControl &cp)
    {
        TPZMHMixedMeshControl::operator=(cp);
        fHDivWrapperMatId = cp.fHDivWrapperMatId;
        return *this;
    }
    
    /// Insert Boundary condition objects that do not perform any actual computation
    virtual void InsertPeriferalMaterialObjects() override;
    
    /// Insert the necessary H(div) material objects to create the flux mesh
    virtual void InsertPeriferalHdivMaterialObjects() override;
    
    /// Insert the necessary pressure material objects to create the pressure mesh
    virtual void InsertPeriferalPressureMaterialObjects() override;
    
    /// Create all data structures for the computational mesh
    virtual void BuildComputationalMesh(bool usersubstructure) override;
    
public:
    
    /// material id of the wrapper elements at the end of H(div) elements
    int fHDivWrapperMatId = 505;
    
    /// material id of the pressure Lagrange multipliers of dimension fDim-1
    int fPressureDim1MatId = 506;
    

protected:
    
    // create the elements domain per domain with approximation spaces disconnected from each other
    virtual void CreateInternalFluxElements();

    // create the approximation space associated with the skeleton and restrain the connects
    virtual void CreateSkeleton() override;
        
    /// Create the interfaces between the pressure elements of dimension dim
    virtual void CreateMultiPhysicsInterfaceElements(int dim);

    /// Create geometric Skeleton elements for flux and pressure in case of hybridization of skeleton fluxes
    void CreateSecondSkeletonGeometricElements();
    
    /// Create geometric elements for the pressure interface elements with mat id fPressureDim1MatId
    // this method relies upon the existence of the HDivWrapMatId elements
    void CreatePressureInterfaceGeometricElements();

    /// Control the numbering of the equations based on the lagrange multiplier level of the connects
    virtual void SetLagrangeMultiplierLevels() override;

public:
    
    
    /// verify the consistency of the datastructure
    //  only implemented in the hybrid Hdiv version
    virtual void CheckMeshConsistency();
    
    /// print the elements in a readable format
    virtual void PrintFriendly(std::ostream &out);


protected:

    /// Create the pressure mesh which is dual to the flux mesh
    virtual void CreatePressureMHMMesh();
    
    /// create the multiphysics mesh
    void CreateMultiphysicsMHMMesh();

    // create the dim-1 pressure elements between the hdiv elements
    // this method now relies on the existence of the geometric elements
    void CreatePressureInterfaces();
    
    /// hybridize the flux elements with the given material id - each flux element creates
    /// a pressure element
    virtual void HybridizeSkeleton(int skeletonmatid, int pressurematid);
    
    // create a boundary element around each flux element to optimize calculations
    void CreateHDivWrappers();
    
    // Find the two connected flux elements to a pressure element
    void FindConnectedElements(TPZGeoElSide &pressureindex, int domain, TPZVec<TPZCompElSide> &fluxconnected);
    
    /// group and condense the elements
    virtual void GroupandCondenseElements();

    /// group element H(div) elements with surrounding interface elements
    void GroupElements(TPZCompMesh *cmesh);
    
};

#endif /* TPZMHMixedHybridMeshControl_hpp */
