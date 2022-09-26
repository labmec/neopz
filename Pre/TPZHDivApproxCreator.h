#ifndef TPZHDIVAPPROXCREATOR_H
#define TPZHDIVAPPROXCREATOR_H

#include "TPZApproxCreator.h"

class TPZMaterial;
class TPZMultiphysicsCompMesh;

/// @class TPZHDivApproxCreator
/// @brief Manages the construction of approximation spaces with HDiv spaces
class TPZHDivApproxCreator : public TPZApproxCreator {

protected:
    
    /// Type of HDiv family to be used. Some of the options are HDivStandard, HDivKernel, and HDivConstant
    HDivFamily fHDivFam = HDivFamily::EHDivStandard;
    
public:
    
    /// Default constructor
    TPZHDivApproxCreator() = default;
    
    /// Constructor with gmesh
    TPZHDivApproxCreator(TPZGeoMesh * gmesh);

    /// Default destructor
    ~TPZHDivApproxCreator();

    /// Get/Set HDiv family
    HDivFamily &HdivFamily(){return fHDivFam;}
    /// Get HDiv family
    const HDivFamily &HdivFamily() const {return fHDivFam;}
    
    /// Driver function. Will create the atomic meshes (HDiv, L2, etc.) and an associate multiphysics mesh
    TPZMultiphysicsCompMesh *CreateApproximationSpace() override;

protected:

    /// Checks if the parameters provided do not violate any condition for mesh generation
    void CheckSetupConsistency() override;
        
    /// Groups the elements in data structure to be condensed
    /// @param mcmesh multiphysics compmesh with elements to be condensed
    void GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh) override;
    
private:

    /// Creates and HDiv approximation space/cmesh
    TPZCompMesh * CreateHDivSpace();
    
    /// Creates an L2 approximation space/cmesh
    /// @param pOrder polynomial order of the space
    /// @param lagLevel lagrange multiplier level for factorization
    TPZCompMesh * CreateL2Space(const int pOrder,const int lagLevel);
    
    /// Creates the multiphysics compmesh based on the vector of atomic meshes
    /// @param meshvector vector of atomic meshes. For instance HDiv and L2 for Darcy problem
    TPZMultiphysicsCompMesh * CreateMultiphysicsSpace(TPZManVector<TPZCompMesh*> meshvector);
    
    /// Creates the rotation space for elasticity problems
    TPZCompMesh * CreateRotationSpace(const int pOrder, const int lagLevel);
};


#endif
