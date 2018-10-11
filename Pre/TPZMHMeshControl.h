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
    
public:

    /// Specify the type of differential equation
    enum MProblemType {ENone, EScalar, EElasticity2D, EElasticity3D};

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
    
    /// Variable defining the type of problem
    MProblemType fProblemType = EScalar;
    
    /// number of state variables
    int fNState = 1;
    
public:
    /// material id associated with the skeleton elements
    int fSkeletonMatId = 4;
    
    /// material id associated with the skeleton elements
    int fSecondSkeletonMatId = 3;
    
    /// material id associated with the skeleton elements in a hybrid context
    int fPressureSkeletonMatId = 6;
    
    /// material id associated with the lagrange multiplier elements
    int fLagrangeMatIdLeft = 550, fLagrangeMatIdRight = 551;
    
    /// materials used for modeling the differential equation
    std::set<int> fMaterialIds;
    
    /// materials for boundary conditions
    std::set<int> fMaterialBCIds;
    
protected:
    
    /// interpolation order of the internal elements
    int fpOrderInternal = 1;
    
    /// interpolation order of the skeleton elements
    int fpOrderSkeleton = 1;
    
    /// material index of the skeleton wrap
    int fSkeletonWrapMatId = 500;
    
    /// material index of the boundary wrap
    int fBoundaryWrapMatId = 502;
    
    /// material index of the internal wrap
    int fInternalWrapMatId = 501;
    
    /// vector of coarse domain index associated with each geometric element
    TPZManVector<int64_t> fGeoToMHMDomain;
    
    /// indices of the geometric elements which define the skeleton mesh and their corresponding subcmesh indices
    std::map<int64_t,int64_t> fMHMtoSubCMesh;
    
    /// indices of the skeleton elements and their left/right geometric elements of the skeleton mesh
    std::map<int64_t, std::pair<int64_t,int64_t> > fInterfaces;
    
    /// geometric index of the connects - subdomain where the connect will be internal
    std::map<TPZCompMesh *,TPZManVector<int64_t> > fConnectToSubDomainIdentifier;
    
    /// flag to determine whether a lagrange multiplier is included to force zero average pressures in the subdomains
    /**
     * when imposing average pressure to be zero, a multiphysics mesh is created
     */
    bool fLagrangeAveragePressure = false;
    
    /// flag to indicate whether we create a hybridized mesh
    bool fHybridize = false;
    
    /// flag to indicate whether the lagrange multipliers should switch signal
    bool fSwitchLagrangeSign = false;
    
public:
    /// number of equations when not condensing anything
    int64_t fGlobalSystemSize;
    
    /// number of equations considering local condensation
    int64_t fGlobalSystemWithLocalCondensationSize;
    
    /// number of equations of the global system
    int64_t fNumeq;
public:
    TPZMHMeshControl() : fProblemType(EScalar), fNState(1), fGeoToMHMDomain(), fMHMtoSubCMesh(),
    fLagrangeAveragePressure(false),
    fHybridize(false), fSwitchLagrangeSign(false), fGlobalSystemSize(-1), fGlobalSystemWithLocalCondensationSize(-1), fNumeq(-1)
    {
        
    }

    TPZMHMeshControl(int dimension) : fProblemType(EScalar), fNState(1), fGeoToMHMDomain(), fMHMtoSubCMesh(),
    fLagrangeAveragePressure(false),
    fHybridize(false), fSwitchLagrangeSign(false), fGlobalSystemSize(-1), fGlobalSystemWithLocalCondensationSize(-1), fNumeq(-1)
    {
        fGMesh = new TPZGeoMesh;
        fGMesh->SetDimension(dimension);
        fCMesh = new TPZCompMesh(fGMesh);
        fPressureFineMesh = new TPZCompMesh(fGMesh);
    }

    /// constructor, indicating that the MHM approximation will use the elements indicated by coarseindices as the macro elements
    //    TPZMHMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<int64_t> &coarseindices);
    
    TPZMHMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &geotomhm);
    
    /// create the mhm object without defining the MHM partition
    TPZMHMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh);
    
    TPZMHMeshControl(const TPZMHMeshControl &copy);
    
    virtual ~TPZMHMeshControl();
    
    TPZMHMeshControl &operator=(const TPZMHMeshControl &cp);
    
    TPZAutoPointer<TPZCompMesh> CMesh() const
    {
        return fCMesh;
    }
    
    TPZAutoPointer<TPZGeoMesh> GMesh() const
    {
        return fGMesh;
    }
    
    TPZAutoPointer<TPZCompMesh> PressureMesh()
    {
        return fPressureFineMesh;
    }
    
    /// Define the MHM partition by the coarse element indices
    void DefinePartitionbyCoarseIndices(TPZVec<int64_t> &coarseindices);
    
    void DefineSkeleton(std::map<int64_t,std::pair<int64_t,int64_t> > &skeleton)
    {
        fInterfaces = skeleton;
#ifdef PZDEBUG
        // the skeleton elements link MHM domain indices (which need to exist)
        for (auto it = fInterfaces.begin(); it != fInterfaces.end(); it++) {
            if (fMHMtoSubCMesh.find(it->second.first) == fMHMtoSubCMesh.end()) {
                DebugStop();
            }
            if (it->second.second != it->first && fMHMtoSubCMesh.find(it->second.second) == fMHMtoSubCMesh.end()) {
                DebugStop();
            }
        }
#endif
    }
    
    /// Define the partitioning information of the MHM mesh
    void DefinePartition(TPZVec<int64_t> &partitionindex, std::map<int64_t,std::pair<int64_t,int64_t> > &skeleton);
    
    /// Set the problem type of the simulation
    void SetProblemType(MProblemType problem)
    {
        fProblemType = problem;
        switch(problem)
        {
            case EScalar:
                fNState = 1;
                break;
            case EElasticity2D:
                fNState = 2;
                break;
            case EElasticity3D:
                fNState = 3;
                break;
            default:
                DebugStop();
        }
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
    virtual void SetHybridize(bool flag = true)
    {
        fHybridize = flag;
    }
    
    virtual bool GetHybridize() {
        return fHybridize;
    }
    /// Create all data structures for the computational mesh
    virtual void BuildComputationalMesh(bool usersubstructure);
    
protected:
    /// will create dim-1 geometric elements on the interfaces between the coarse element indices
    virtual void CreateSkeletonElements();
    
    /// Insert material objects that do not perform any actual computation
    virtual void InsertPeriferalMaterialObjects();
public:
    /// divide the skeleton elements
    void DivideSkeletonElements(int ndivide);
    
    /// divide the boundary skeleton elements
    void DivideBoundarySkeletonElements();
    
    /// switch the sign of the lagrange multipliers
    void SwitchLagrangeMultiplierSign(bool sw)
    {
        fSwitchLagrangeSign = sw;
    }
    
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
    std::map<int64_t,int64_t> &Coarse_to_Submesh()
    {
        return fMHMtoSubCMesh;
    }
    
    /// verify the consistency of the datastructure
    //  only implemented in the hybrid Hdiv version
    virtual void CheckMeshConsistency()
    {
        
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
    void CreateInterfaceElements2();
    
    /// hybridize the flux elements with the given material id - each flux element creates
    /// a pressure element
    virtual void HybridizeSkeleton(int skeletonmatid, int pressurematid);
    
    /// verify if the element is a sibling of
    bool IsSibling(int64_t son, int64_t father);
    
    /// put the element side which face the boundary on the stack
    void AddElementBoundaries(int64_t elseed, int64_t compelindex, TPZStack<TPZCompElSide> &result);
    
    /// create the lagrange multiplier mesh, one element for each subdomain
    void CreateLagrangeMultiplierMesh();
    
    /// transform the computational mesh into a multiphysics mesh
    void TransferToMultiphysics();
    
    /// substructure the mesh
    void SubStructure();
    void SubStructure2();
    
    /// Create the wrap elements
    void BuildWrapMesh(int dim);
    
    /// Verify if the element side contains a wrap neighbour
    int HasWrapNeighbour(TPZGeoElSide gelside);
    
    /// Verify if the mesh datastructure is consistent
    /// if the current gelside has no father, then none of its neighbours should either
    void CheckDivisionConsistency(TPZGeoElSide gelside);
    
    /// Return the wrap material id (depends on being boundary, neighbour of skeleton or interior
    int WrapMaterialId(TPZGeoElSide gelside);
    
    /// Return true if the material id is related to a skeleton
    virtual bool IsSkeletonMatid(int matid)
    {
        return matid == fSkeletonMatId;
    }
    
    /// return true if the material id is related to a boundary
    virtual bool IsBoundaryMatid(int matid)
    {
        return fMaterialBCIds.find(matid) != fMaterialBCIds.end();
    }
    
    /// CreateWrapMesh of a given material id
    void CreateWrap(TPZGeoElSide gelside);
    
    /// CreateWrapMesh of a given material id
    void CreateWrap(TPZGeoElSide gelside, int wrapmaterial);
    
    /// Divide the wrap element while it has divided neighbours
    void DivideWrap(TPZGeoEl *wrapelement);
    
protected:
    /// associates the connects of an element with a subdomain
    void SetSubdomain(TPZCompEl *cel, int64_t subdomain);
    
    /// associates the connects index with a subdomain
    void SetSubdomain(TPZCompMesh *cmesh, int64_t connectindex, int64_t subdomain);
    
    /// returns to which subdomain a given element beint64_ts
    // this method calls debugstop if the element beint64_ts to two subdomains
    int64_t WhichSubdomain(TPZCompEl *cel);
    
    /// Subdomains are identified by computational mesh, this method will join
    // the subdomain indices of the connects to the multi physics mesh
    void JoinSubdomains(TPZVec<TPZCompMesh *> &meshvec, TPZCompMesh *multiphysicsmesh);
    
    /// print the diagnostics for a subdomain
    void PrintSubdomain(int64_t elindex, std::ostream &out);
    
    /// print the indices of the boundary elements and interfaces
    void PrintBoundaryInfo(std::ostream &out);
    
    /// identify connected elements to the skeleton elements
    // the computational mesh is determined by the element pointed to by the geometric element
    void ConnectedElements(int64_t skeleton, std::pair<int64_t,int64_t> &leftright, std::map<int64_t, std::list<TPZCompElSide> > &ellist);

    /// identify interface elements connected to the skeleton elements
    // the computational mesh is determined by the element pointed to by the geometric element
    void ConnectedInterfaceElements(int64_t skeleton, std::pair<int64_t,int64_t> &leftright, std::map<int64_t, std::list<TPZInterfaceElement *> > &ellist);
};

#endif /* defined(__PZ__TPZMHMeshControl__) */
