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
    
    TPZMHMixedHybridMeshControl(int dimension) : TPZMHMixedMeshControl(dimension)
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
    
    /// material id of the pressure Lagrange multipliers of dimension fDim-2
    int fPressureDim2MatId = 507;
    
    /// material id of the flow elements of dimension fDim-1
    std::set<int> fFractureFlowDim1MatId;
    
    /// Return true if the material id is related to a skeleton
    virtual bool IsSkeletonMatid(int matid) override
    {
        return matid == fSkeletonMatId || fSkeletonWithFlowMatId.find(matid) != fSkeletonWithFlowMatId.end();
    }

    /// material id of zero flux boundary condition, used for fracture elements
    int fHomogeneousNeumannBcMatId = 509;
    
    /// material id of skeleton elements which have axial flow
    std::set<int> fSkeletonWithFlowMatId;
    
    /// material id of the pressure corresponding to the skeleton with axial flow
    int fSkeletonWithFlowPressureMatId = 12;
    
protected:
    // create the elements domain per domain with approximation spaces disconnected from each other
    virtual void CreateInternalFluxElements() override;

    // create the approximation space associated with the skeleton and restrain the connects
    virtual void CreateSkeleton() override;
        
    /// Create the interfaces between the pressure elements of dimension dim
    virtual void CreateMultiPhysicsInterfaceElements(int dim) override;


public:
    
    void InsertFractureFlowMaterial(int matid)
    {
        if (fFractureFlowDim1MatId.find(matid) != fFractureFlowDim1MatId.end()) {
            DebugStop();
        }
        fFractureFlowDim1MatId.insert(matid);
    }
    
    // create HDiv approximation spaces associated with the pressure skeleton elements
    // initially we will create only HDiv elements associated with the coarse skeleton elements
    void CreateSkeletonAxialFluxes();
    
    // create HDiv approximation spaces associated with the pressure elements on internal side
    // initially we will create only HDiv elements associated with the coarse skeleton elements
    void CreateInternalAxialFluxes();
    
    /// Add axial flux to a pressure element
    /// the index of the associated geometric element is given and determines the material id of the element
    // the material id of the pressure element can be either fPressureSkeletonMatId or fPressureDim1MatId
    // An H(div) and a pressure element will be created with material Id of the geometric element
    // HDivWrapper boundary elements will also be created
    void CreateAxialFluxElement(TPZInterpolatedElement *PressureElement, int gelfluxmatid);
    
    /// split the fluxes between the flux elements adjacent to a fracture
    void SplitFluxElementsAroundFractures();
    

    /// Change the material id of the boundary elements associated with fracture flow
    // if the element is neighbour of a pressure boundary condition, maintain the pressure bc
    // if the element is neighbour of a flux boundary condition, apply a homogeneous flux condition
    void AdjustBoundaryConditionsOfFractures();
    
    /// verify the consistency of the datastructure
    //  only implemented in the hybrid Hdiv version
    virtual void CheckMeshConsistency() override;
    
    /// print the elements in a readable format
    virtual void PrintFriendly(std::ostream &out) override;

    /// Set the hybridization to true
    virtual void SetHybridize(bool flag) override
    {
        fHybridize = flag;
    }

protected:

    /// Create the pressure mesh which is dual to the flux mesh
    virtual void CreatePressureMHMMesh() override;
    

    // create the dim-1 pressure elements between the hdiv elements
    // this method now relies on the existence of the geometric elements
    void CreatePressureInterfaces();
    
    /// hybridize the flux elements with the given material id - each flux element creates
    /// a pressure element
    virtual void HybridizeSkeleton(int skeletonmatid, int pressurematid) override;
    
    /// Create lower dimension pressure elements (dim-2)
    /// The (dim-2) geometric elements have already been created.
    // when dim==2 create a point pressure element for each skeleton element
    // only create internal elements - we do not hybridize the contours
    // associate an L2 do nothing material
    void CreateLowerDimensionPressureElements();
    
    /// verifies if a HDiv wrapper needs to be created for a given element/side
    // the method will look if there is a neighbouring geometric element with material id fPressureDim1MatId
    bool NeedsHDivWrapper(TPZGeoElSide gelside);
    
    // create a boundary element around each flux element to optimize calculations
    void CreateHDivWrappers();
    
    // Find the two connected flux elements to a pressure element
    void FindConnectedElements(TPZGeoElSide &pressureindex, int domain, TPZVec<TPZCompElSide> &fluxconnected);
    
    /// group and condense the elements
    virtual void GroupandCondenseElements() override;

    /// group element H(div) elements with surrounding interface elements
    void GroupElements(TPZCompMesh *cmesh);
    
    /// for a given boundary wrap element
    /// determine if it has a neighbouring fracture element
    /// determine the appropriate boundary condition
    /// change the material id of the geometric element to apply the boundary condition
    void ApplyNeighbourBoundaryCondition(TPZGeoEl *gel);
    
};

#endif /* TPZMHMixedHybridMeshControl_hpp */
