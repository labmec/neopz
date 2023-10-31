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
#include "TPZEnumApproxFamily.h"

/// class for creating TPZMHMM with Mixed Meshes
class TPZMHMixedMeshControl : public TPZMHMeshControl
{
    
protected:
    
    /// computational mesh to contain the rotation elements
    TPZAutoPointer<TPZCompMesh> fRotationMesh;
    
    /// distributed flux per element
    TPZAutoPointer<TPZCompMesh> fElementFluxMesh;
    /// average state value per element
    TPZAutoPointer<TPZCompMesh> fElementAverageMesh;
    
    HDivFamily fHDivFamily = HDivFamily::EHDivStandard;
public:
    
    enum LagrangeLevels {BoundHdiv = 0, CornerPressure = 1, ElementFLux = 2, Pressure = 3, Rotation = 4, InternalHdiv = 5, Elpressure = 6, DomainFlux = 7, OneElPressure = 8, SkeletonHdiv = 9, DomainPressure = 10 };
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
    
    
    /// Create all data structures for the computational mesh
    virtual void BuildComputationalMesh(bool usersubstructure) override;
    
    /// Set the flag for creating Lagrange Dofs for the average pressure
    void SetLagrangeAveragePressure(bool flag)
    {
        if (flag == true) {
            DebugStop();
        }
    }
    
    /// set the HDiv Family of the HDiv space
    void SetHDivFamily(HDivFamily fam)
    {
        if(fam == HDivFamily::EHDivKernel) DebugStop();
        fHDivFamily = fam;
    }
    
    HDivFamily HdivFamily() const
    {
        return fHDivFamily;
    }
    
    /// Set the hybridization to true
    virtual void Hybridize(bool flag)
    {
        DebugStop();
    }

    /// Put the pointers to the meshes in a vector
    void GetMeshVec(TPZVec<TPZCompMesh *> &meshvec)
    {
        meshvec.Resize(7);
        meshvec[0] = fFluxMesh.operator->();
        meshvec[1] = fPressureFineMesh.operator->();
        int count = 2;
        if(fRotationMesh)
        {
            meshvec[count++] =  fRotationMesh.operator->();
        }
        if(fElementFluxMesh)
        {
            meshvec[count++] = fElementFluxMesh.operator->();
            meshvec[count++] = fElementAverageMesh.operator->();
        }
        if(fCMeshDomainFlux)
        {
            meshvec[count++] = fCMeshDomainFlux.operator->();
            meshvec[count++] = fCMeshDomainPressure.operator->();
        }
        meshvec.Resize(count);
    }

    TPZVec<TPZAutoPointer<TPZCompMesh> > GetMeshes()
    {
        TPZManVector<TPZAutoPointer<TPZCompMesh>,7> result(7);
        result[0] = fFluxMesh;
        result[1] = fPressureFineMesh;
        int count = 2;
        if(fRotationMesh)
        {
            result[count++] =  fRotationMesh;
        }
        if(fElementFluxMesh)
        {
            result[count++] = fElementFluxMesh;
            result[count++] = fElementAverageMesh;
        }
        if(fCMeshDomainFlux)
        {
            result[count++] = fCMeshDomainFlux;
            result[count++] = fCMeshDomainPressure;
        }
        result.Resize(count);
        return result;
    }
    
    /// Create the meshes which represent element distributed flux and element average
    void CreateAverageElementMeshes();
    
    /// print the data structure
    void Print(std::ostream &out);

    /// print in a user friendly manner
    virtual void PrintFriendly(std::ostream &out)
    {
        
    }
    
    // set the Lagrange levels of the connects for the atomic meshes
    void InitializeLagrangeLevels();
    
    void ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder);
    
protected:
    
   

    /// Create the mesh of the flux approximation space
    void CreateFluxMesh();
    
    /// Create the pressure mesh which is dual to the flux mesh
    virtual void CreatePressureMHMMesh();
    
    /// Create the rotation mesh to elasticity problem
    virtual void CreateRotationMesh();
    
    // create the elements domain per domain with approximation spaces disconnected from each other
    virtual void CreateFluxElements();
    
    

    /// Create the multiphysics mesh
    void CreateHDivPressureMHMMesh();
    
    /// build the multi physics mesh (not at the finest geometric mesh level
    void BuildMultiPhysicsMesh();

    /// put the elements in TPZSubCompMesh, group the elements and condense locally
    void HideTheElements();

    /// hybridize the flux elements with the given material id - each flux element creates
    /// a pressure element
//    virtual void HybridizeSkeleton(int skeletonmatid, int pressurematid);
    
    /// group and condense the elements
    virtual void GroupandCondenseElements();
    
};

#endif /* TPZMHMixedMeshControl_hpp */
