//
//  TPZHybridizeHDiv.hpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 16/05/18.
//

#ifndef TPZHybridizeHDiv_hpp
#define TPZHybridizeHDiv_hpp

#include <stdio.h>
#include <map>

template<class T>
class TPZVec;

class TPZCompMesh;
class TPZMaterial;
class TPZCompElSide;
class TPZInterpolatedElement;

struct TPZHybridizeHDiv {
    int HDivWrapMatid = -10;
    int LagrangeInterface = -9;
    int InterfaceMatid = -8;
    int NState = 1;
    
    TPZHybridizeHDiv() = default;
    
    TPZHybridizeHDiv(TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
    void ComputeNState(TPZVec<TPZCompMesh*>& meshvec_Hybrid);

    /// compute material ids for the periferal material objects
    void ComputePeriferalMaterialIds(TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
    /// return true if a material id is a peripheral material
    bool IsPeriferalMaterialId(int matid)
    {
        return matid == HDivWrapMatid || matid == LagrangeInterface || matid == InterfaceMatid;
    }
    /// split the connects between flux elements and create a dim-1 pressure element
    void HybridizeInternalSides(TPZVec<TPZCompMesh *> &meshvec_Hybrid);

    /// Create interface elements with material id InterfaceMatid
    void CreateInterfaceElements(TPZCompMesh *cmesh_Hybrid, TPZVec<TPZCompMesh *> &meshvec_Hybrid);

    /// create a multiphysics mesh for the hybrid formulation using the materials of another mesh and the given atomic meshes
    TPZCompMesh * CreateMultiphysicsMesh(TPZCompMesh *cmesh_HDiv, TPZVec<TPZCompMesh *> &meshvec_Hybrid, double Lagrange_term_multiplier = 1.);
    
    /// group and condense the elements
    static void GroupElements(TPZCompMesh *cmesh_Hybrid);
    
    /// insert the material objects for HDivWrap and LagrangeInterface in the atomic meshes
    void InsertPeriferalMaterialObjects(TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
    /// insert the material objects for HDivWrap, LagrangeInterface and InterfaceMatid in the multiphysics mesh
    void InsertPeriferalMaterialObjects(TPZCompMesh *cmesh_Hybrid, double Lagrange_term_multiplier = 1.);
    
    /// clones the atomic meshes in meshvec_HDiv and creates a multiphysics hybrid mesh
    std::tuple<TPZCompMesh *, TPZVec<TPZCompMesh *> > Hybridize(TPZCompMesh *cmesh_Multiphysics, TPZVec<TPZCompMesh *> &meshvec_HDiv, bool group_elements=true, double Lagrange_term_multiplier = 1.);
    
private:
    
    std::tuple<int64_t,int> SplitConnects(const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> &meshvec_Hybrid);

    static TPZCompElSide RightElement(TPZInterpolatedElement *intel, int side);

};

#endif /* TPZHybridizeHDiv_hpp */
