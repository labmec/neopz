

#ifndef TPZHDIVAPPROXCREATOR_H
#define TPZHDIVAPPROXCREATOR_H

#include "TPZEnumApproxFamily.h"
#include "pzstack.h"
#include "pzgmesh.h"

class TPZMaterial;
class TPZMultiphysicsCompMesh;

enum class HybridizationType {ENone, EStandard, ESemi};
enum class ProblemType {ENone, EElastic, EDarcy, EStokes};

/// @class TPZHDivApproxCreator
/// @brief Manages the construction of approximation spaces with HDiv spaces
class TPZHDivApproxCreator
{
protected:
    
    /// Type of HDiv family to be used. Some of the options are HDivStandard, HDivKernel, and HDivConstant
    HDivFamily fHDivFam = HDivFamily::EHDivStandard;
    
    /// Sets what type of hybridization should be used between elements, if any
    HybridizationType fHybridType = HybridizationType::ENone;
    
    /// Variable controling the problem (Darcy, elasticity, etc.). This may influence the number of compmeshes
    ProblemType fProbType = ProblemType::ENone;

    /// Vector of materials. Will be used to to create the compmeshes
    TPZStack<TPZMaterial * > fMaterialVec;

    /// Pointer to geomesh
    TPZGeoMesh *fGeoMesh = nullptr;

    /// Default porder to be used by all the elements in the cmeshes
    int fDefaultPOrder = 1;
    
    /// @brief This sets hdiv+/hdiv++ internal/bubble POrder: 0- hdiv, 1- hdiv+, 2- hdiv++
    int fExtraInternalPOrder = 0; 

    /// Number of compmeshes which is related to the number of state variables. Default is Darcy that has Flux (Hdiv) and Pressure (L2) cmeshes
    int fNumMeshes = 2;

    /// Enhances space with constant fields that allow for the condensation of internal degrees of freedom
    bool fIsEnhancedSpaces = false;

public:
    
    /// Default constructor
    TPZHDivApproxCreator() = default;
    
    /// Constructor with gmesh
    TPZHDivApproxCreator(TPZGeoMesh * gmesh);

    /// Default destructor
    ~TPZHDivApproxCreator();
    
    /// Inserts a domain material to the vector of materials of this class. These are used to build the cmeshes
    /// @param mat pointer to materials to be added
    int InsertMaterialObject(TPZMaterial * mat);
    
    /// Inserts a bc material to the vector of materials of this class. These are used to build the cmeshes
    /// @param mat boundary condition material to be added
    int InsertMaterialObject(TPZBndCond * mat);

    /// Get/Set GeoMesh
    TPZGeoMesh *GeoMesh(){return fGeoMesh;}
    /// Get GeoMesh
    const TPZGeoMesh *GeoMesh() const {return fGeoMesh;}

    /// Get/Set HDiv family
    HDivFamily &HdivFamily(){return fHDivFam;}
    /// Get HDiv family
    const HDivFamily &HdivFamily() const {return fHDivFam;}

    /// Get/Set Hybridization type
    HybridizationType &HybridType(){return fHybridType;}
    /// Get Hybridization type
    const HybridizationType &HybridType() const {return fHybridType;}
    
    /// Get/Set Enhanced spaces
    bool &EnhancedSpaces(){return fIsEnhancedSpaces;}
    /// Get Enhanced spaces
    const bool &EnhancedSpaces() const {return fIsEnhancedSpaces;}

    /// Get/Set Problem type
    ProblemType &ProbType(){return fProbType;}
    /// Get Problem type
    const ProblemType &ProbType() const {return fProbType;}

    /// Set cmeshes default polynomial order
    void SetDefaultOrder(const int ord){
        fDefaultPOrder = ord;
    }
    
    /// Get cmeshes default polynomial order
    int &GetDefaultOrder(){return fDefaultPOrder;}

    /// Set extra internal porder
    void SetExtraInternalOrder(const int ord){
        fExtraInternalPOrder = ord;
    }
    /// Get extra internal porder
    int &GetExtraInternalOrder(){return fExtraInternalPOrder;}
    
    /// Driver function. Will create the atomic meshes (HDiv, L2, etc.) and an associate multiphysics mesh
    TPZMultiphysicsCompMesh *CreateApproximationSpace();

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
    TPZCompMesh * CreateRotationSpace();
   

};




#endif
