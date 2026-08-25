//
// Created by victor on 21/09/2022.
//

#ifndef PZ_TPZAPPROXCREATOR_H
#define PZ_TPZAPPROXCREATOR_H

#include "TPZEnumApproxFamily.h"
#include "pzstack.h"
#include "pzgmesh.h"

class TPZMaterial;
class TPZMultiphysicsCompMesh;
class TPZInterpolationSpace;
class TPZMultiphysicsElement;

enum class HybridizationType {ENone, EStandard, EStandardSquared, ESemi};
enum class ProblemType {ENone, EElastic, EDarcy, EStokes};

/// @struct TPZHybrid
/// @brief stores the geometric elements involved in a hybridization
struct TPZHybrid {
    /// geometric element associated with the lagrange multiplier
    int64_t fLagrange;
    /// left geometric element
    int64_t fLeft;
    /// right geometric element
    int64_t fRight;
    
    TPZHybrid() : fLagrange(-1), fLeft(-1), fRight(-1) {
        
    }
    TPZHybrid(int64_t lagrange, int64_t left, int64_t right) : fLagrange(lagrange), fLeft(left), fRight(right)
    {
        
    }
};
/// @class TPZApproxCreator
/// @brief Manages the construction of approximation spaces
class TPZApproxCreator {
protected:

    /// Sets what type of hybridization should be used between elements, if any
    HybridizationType fHybridType = HybridizationType::ENone;

    /// Variable controling the problem (Darcy, elasticity, etc.). This may influence the number of compmeshes
    ProblemType fProbType = ProblemType::ENone;

    /// Vector of materials. Will be used to to create the computational meshes
    std::map<int,TPZMaterial *> fMaterialVec;

    /// number of state variables
    int fNState = 1;
    /// Pointer to geomesh
    TPZGeoMesh *fGeoMesh = nullptr;

    /// Default porder to be used by all the elements in the cmeshes
    int fDefaultPOrder = 1;

    /// @brief This sets hdiv+/hdiv++ internal/bubble POrder: 0- hdiv, 1- hdiv+, 2- hdiv++
    int fExtraInternalPOrder = 0;

    /// Number of compmeshes which is related to the number of state variables. Default is Darcy that has Flux (Hdiv) and Pressure (L2) cmeshes
    int fNumMeshes = 2;

    /// Boolean telling internal equations should be condensed
    bool fShouldCondense = 1;

    /// Enhances space with constant fields that allow for the condensation of internal degrees of freedom
    bool fIsRBSpaces = false;

public:
    /// Struct with all the data regarding hybridization between elements
    struct HybridizationData{
        
        HybridizationData(const HybridizationData &copy) = default;
        
        HybridizationData() {
            
        }
        
        HybridizationData &operator=(const HybridizationData &copy) = default;
        
        /// compute the lagrange multiplier coeficients as a function of the problem type and hybridization
        void SetProblemHybridH1(ProblemType prob, HybridizationType hybrid);
        
        /// set the lagrange multiplier coefficients for hybrid Hdiv problems
        /// compute the lagrange multiplier coeficients as a function of the problem type and hybridization
        void SetProblemHybridHDiv(ProblemType prob, HybridizationType hybrid);

        /// Print the hybridization data for debugging
        void Print(std::ostream &out = std::cout) const;

        /// Matid of element at the border or higher domain element
        int fWrapMatId = -123456; // to check if it has been initialized

        /// Matid of element at the interface between elements
        int fLeftInterfaceMatId = -123456;
        int fRightInterfaceMatId = -123456;

        /// Matid of lagrange multiplier element
        int fLagrangeMatId = -123456;

        /// Material ids generated due second hybridization
        int fSecondLeftInterfaceMatId = -123456;

        /// Material ids generated due second hybridization
        int fSecondRightInterfaceMatId = -123456;

        /// Matid of second lagrange multiplier, useful in double hybridizations
        int fSecondLagrangeMatId = -123456;

        /// indicates whether the boundary conditions should be hybridized and how many times it should be
        int fHybridizeBCLevel = 0;
        
        /// Vector of multipliers associated with the interface elements.
        /// The multipliers are influenced by the type of problem and whether the hybridization is simple or squared
        /// The order is left, right, second left, second right
        TPZManVector<REAL,4> fMultipliers = {0.,0.,0.,0.};
        
        /// map linking the element index of an interface to left and right elements
        std::map<int64_t,TPZHybrid> fInterfaces;
    };

protected:
    /// Attribute of struct with all the data regarding hybridization between elements
    HybridizationData fHybridizationData;
    
    /// Add a hybrid geometric configuration
    virtual void AddHybridGeometry(int64_t interface_id, TPZHybrid &hybrid);

    /// Predominant GeoMesh dimension element
    MElementType fElementType = MElementType::ENoType;

public:
    
    /// Default constructor
    TPZApproxCreator() = default;

    /// Constructor with gmesh
    TPZApproxCreator(TPZGeoMesh * gmesh);

    /// Default destructor
    ~TPZApproxCreator() = default;

    /// Inserts a domain material to the vector of materials of this class. These are used to build the cmeshes
    /// @param mat pointer to materials to be added
    int InsertMaterialObject(TPZMaterial * mat);

    /// Inserts a bc material to the vector of materials of this class. These are used to build the cmeshes
    /// @param mat boundary condition material to be added
    int InsertMaterialObject(TPZBndCond * mat);
    
    void SetNState(int nstate) {
        fNState = nstate;
    }
    
    int NState() {
        return fNState;
    }

    /// Get/Set GeoMesh
    TPZGeoMesh *GeoMesh(){return fGeoMesh;}
    /// Get GeoMesh
    const TPZGeoMesh *GeoMesh() const {return fGeoMesh;}

    /// Get/Set Hybridization type
    HybridizationType HybridType() const {return fHybridType;}
    /// Set Hybridization type
    virtual void SetHybridType(HybridizationType hybrid) = 0;

    /// Get/Set Enhanced spaces
    bool &IsRigidBodySpaces(){return fIsRBSpaces;}
    /// Get Enhanced spaces
    const bool &IsRigidBodySpaces() const {return fIsRBSpaces;}
    
    const HybridizationData &HybridData() const {return fHybridizationData;}
    void SetHybridData(const HybridizationData &copy) {
        fHybridizationData = copy;
    }

    /// Get/Set Problem type
    virtual void SetProbType(ProblemType prob) = 0;
    
    /// Get Problem type
    const ProblemType ProbType() const {return fProbType;}

    /// return number of meshes
    const int NumMeshes() const {return fNumMeshes;}

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

    /// Set if should condense internal equations
    void SetShouldCondense( bool isCondensed){
        fShouldCondense = isCondensed;
    }

    /// Get if should condense internal equations
    bool GetShouldCondense( bool isCondensed){
        return fShouldCondense;
    }

    /// Driver function. Will create the atomic meshes (HDiv, L2, etc.) and an associate multiphysics mesh. Should be implemented in son classes
    virtual TPZMultiphysicsCompMesh *CreateApproximationSpace() = 0;

    /// Print attributes of class
    void Print(std::ostream &out = std::cout);
    
    /// Sets if should hybridize boundary
    void SetHybridizeBoundary() {
        if(fHybridType == HybridizationType::EStandardSquared) fHybridizationData.fHybridizeBCLevel = 1;
        else if(fHybridType == HybridizationType::EStandard || fHybridType == HybridizationType::ESemi) fHybridizationData.fHybridizeBCLevel = 1;
        else{
            std::cout << "ERROR! Please set a hybridization type before calling SetToHybridizeBoundary()" << std::endl;
            DebugStop();
        }
    }

protected:
    
    /// Checks if the hybridization neighborhood is correct
    void CheckNeighborhoodHybridization() const;

    /// Compute Periferal Material ids for additional elements generated during hybridization
    void ComputePeriferalMaterialIds(int base = 10);

    /// Adds the geo els related to the hybridization
    virtual void AddHybridizationGeoElements();
    
    /// Add geometric elements to represent squared hybridization
    virtual void AddHybridSquareGeoElements();

    ///This method checks if the current configuration is valid
    virtual void CheckSetupConsistency() = 0;

    ///Get BC material ids
    std::set<int> GetBCMatIds();

    ///Get material ids from which has the same dimension as the mesh
    std::set<int> GetVolumetricMatIds();

    /// Insert periferal material objects related to geometric. objects created during hybridization (for wrap and lagrange objects)
    void InsertWrapAndLagrangeMaterialObjects(TPZMultiphysicsCompMesh *mphys);

    /// Insert interface periferal material objects related to geometric objects created during hybridization
    virtual void InsertInterfaceMaterialObjects(TPZMultiphysicsCompMesh *mphys);

    /// Groups elements to be condensed
    virtual void GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh) = 0;

    /// Changes the internal order of the volumetric connect of all CompEls in a given CMesh
    void ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder) const;

public:
    typedef std::map<int64_t, int64_t> CtoMFCel;
    /// @brief Compute the constraints to orthogonalize the restraints
    // this function will compute restraints for all boundary flux connects
    static void ComputeOrthogonalizingRestraints(TPZMultiphysicsCompMesh &mfmesh, CtoMFCel &geltogel, const TPZApproxCreator::HybridizationData &hybridData);

        /// Compute projection directions for a geometric element
        //// Method related to the Semi hybridization.
    static void ProjectionDirections(TPZGeoEl *gel, TPZFMatrix<REAL> &projdir);

    /// Compute projection values for a set of relative coordinates
        //// Method related to the Semi hybridization.
    static void ProjectionValues(TPZVec<REAL> &xrelative, TPZFMatrix<REAL> &projdir, int ncorner, TPZFMatrix<REAL> &funcval);

    /// Compute the projection matrix for an interpolation space
        //// Method related to the Semi hybridization.
    static void ComputeProjectionMatrix(TPZInterpolationSpace *intel, TPZFMatrix<REAL> &projection);

    /// @brief compute the restraints of a connect for a geometric element
        //// Method related to the Semi hybridization.
    static void RestraintConnect(TPZInterpolationSpace *intel, TPZMultiphysicsElement *mfcel, int64_t newind1, int64_t newind2);

    /// @brief Copy the restraint of a connect to another connect
        //// Method related to the Semi hybridization.
    static void CopyRestraintConnect(TPZConnect &cfrom, TPZConnect &cto);

    /// @brief compute an equivalent B matrix
        //// Method related to the Semi hybridization.
    static void ComputeBMatrix(TPZInterpolationSpace *intel, TPZFMatrix<STATE> &B);

    /// @brief  compute orthogonalizing restraint
        //// Method related to the Semi hybridization.
    static void ComputeOrthogonalizingR(TPZFMatrix<REAL> &B, TPZFMatrix<REAL> &Restraint);


public:
    /// Determine the predominant element type in the GeoMesh
    void SetMeshElementType();
};


#endif //PZ_TPZAPPROXCREATOR_H
