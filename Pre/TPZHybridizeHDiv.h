//
//  TPZHybridizeHDiv.hpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 16/05/18.
//

#ifndef TPZHybridizeHDiv_hpp
#define TPZHybridizeHDiv_hpp

#include <stdio.h>
#include <iostream>
#include <map>

template<class T>
class TPZVec;

class TPZCompMesh;
class TPZMaterial;
class TPZCompElSide;
class TPZInterpolatedElement;
class TPZMultiphysicsCompMesh;

struct TPZHybridizeHDiv {
    
    // boundary condition for flux with zero contribution
    int fHDivWrapMatid = -10;
    // material id of the lagrange interface (either pressure or displacement)
    int fLagrangeInterface = -9;
    // material id of the interface elements
    int fInterfaceMatid = -8;
    // number of state variables
    int fNState = 1;
    
    TPZHybridizeHDiv() = default;
    
    TPZHybridizeHDiv(TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
    void ComputeNState(TPZVec<TPZCompMesh*>& meshvec_Hybrid);

    /// compute material ids for the periferal material objects
    void ComputePeriferalMaterialIds(TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
    /// set the periferal material ids
    void SetPeriferalMaterialIds(int HDivWrapMatid, int LagrangeInterface, int InterfaceMatid)
    {
        fHDivWrapMatid = HDivWrapMatid;
        fLagrangeInterface = LagrangeInterface;
        fInterfaceMatid = InterfaceMatid;
        
    }
    
    /// return true if a material id is a peripheral material
    bool IsPeriferalMaterialId(int matid)
    {
        return matid == fHDivWrapMatid || matid == fLagrangeInterface || matid == fInterfaceMatid;
    }
    /// split the connects between flux elements and create a dim-1 pressure element
    void HybridizeInternalSides(TPZVec<TPZCompMesh *> &meshvec_Hybrid);

    /// Create interface elements with material id InterfaceMatid
    void CreateInterfaceElements(TPZCompMesh *cmesh_Hybrid, TPZVec<TPZCompMesh *> &meshvec_Hybrid);

    /// Create interface elements with material id InterfaceMatid using the first and second
    /// meshes
    void CreateInterfaceElements(TPZMultiphysicsCompMesh *cmesh_Hybrid);
    
    /// create a multiphysics mesh for the hybrid formulation using the materials of another mesh and the given atomic meshes
    TPZCompMesh * CreateMultiphysicsMesh(TPZCompMesh *cmesh_HDiv, TPZVec<TPZCompMesh *> &meshvec_Hybrid, double Lagrange_term_multiplier = 1.);
    
    /// create a multiphysics hybridized mesh based on and input mesh
    TPZCompMesh * CreateMultiphysicsMesh(TPZMultiphysicsCompMesh *cmesh_HDiv, double Lagrange_term_multiplier = 1.);

    /// create a multiphysics hybridized mesh based on and input mesh
    void ReCreateMultiphysicsMesh(TPZMultiphysicsCompMesh *cmesh_HDiv, double Lagrange_term_multiplier = 1.);
    
    /// Associate elements with a volumetric element
    // elementgroup[el] = index of the element with which the element should be grouped
    // this method only gives effective result for hybridized hdiv meshes
    static void AssociateElements(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup);
    
    /// group and condense the elements
    static void GroupandCondenseElements(TPZCompMesh *cmesh_Hybrid);
    
    /// insert the material objects for HDivWrap and LagrangeInterface in the atomic meshes
    void InsertPeriferalMaterialObjects(TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
    /// insert the material objects for HDivWrap, LagrangeInterface and InterfaceMatid in the multiphysics mesh
    void InsertPeriferalMaterialObjects(TPZCompMesh *cmesh_Hybrid, double Lagrange_term_multiplier = 1.);
    
    /// clones the atomic meshes in meshvec_HDiv and creates a multiphysics hybrid mesh
    std::tuple<TPZCompMesh *, TPZVec<TPZCompMesh *> > Hybridize(TPZCompMesh *cmesh_Multiphysics, TPZVec<TPZCompMesh *> &meshvec_HDiv, bool group_elements=true, double Lagrange_term_multiplier = 1.);
    
    /// make a hybrid mesh from a H(div) multiphysics mesh
    TPZMultiphysicsCompMesh *Hybridize(TPZMultiphysicsCompMesh *multiphysics, bool group_elements=true, double Lagrange_term_multiplier = 1.);
    
    /// make a hybrid mesh from a H(div) multiphysics mesh
    void HybridizeGivenMesh(TPZMultiphysicsCompMesh &multiphysics, bool group_elements=true, double Lagrange_term_multiplier = 1.);
    
    /// verify the consistency of the solution of the flux mesh
    static void VerifySolutionConsistency(TPZCompMesh *fluxmesh, std::ostream &out);
private:
    
    std::tuple<int64_t,int> SplitConnects(const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> &meshvec_Hybrid);

public:
    
    static TPZCompElSide RightElement(TPZInterpolatedElement *intel, int side);

};

#endif /* TPZHybridizeHDiv_hpp */
