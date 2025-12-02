//
//  TPZMultiphysicsCompMesh.h
//  pz
//
//  Created by Omar Durán on 3/21/19.
//

#ifndef TPZMultiphysicsCompMesh_h
#define TPZMultiphysicsCompMesh_h

#include <stdio.h>
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include "TPZMultiphysicsInterfaceEl.h"

class TPZMultiphysicsCompMesh : public TPZCompMesh {
    
    /// Vector of active physics: index vector
    /// Define wich space will be active in order to generate equations. Should be defined for each space that you want to use by 0: no active or 1: active
    ///The size have to be the same as the m_mesh_vector
    TPZManVector<int,5> m_active_approx_spaces;
    
    /// Vector of computational meshes
    TPZManVector<TPZCompMesh * , 7> m_mesh_vector;
    
public:
    
    /// Default constructor
    TPZMultiphysicsCompMesh();
    
    /// Constructor based on TPZGeoMesh pointer and vector of meshes
    TPZMultiphysicsCompMesh(TPZGeoMesh * gmesh, bool isComplex=false);
    
    /// Constructor based on TPZGeoMesh pointer and vector of meshes
    TPZMultiphysicsCompMesh(TPZAutoPointer<TPZGeoMesh>  gmesh, bool isComplex=false);
    
    /// Copy constructor
    TPZMultiphysicsCompMesh(const TPZMultiphysicsCompMesh &other);
    
    /// Assignement constructor
    TPZMultiphysicsCompMesh & operator=(const TPZMultiphysicsCompMesh &other);
    
    /// Destructor
    ~TPZMultiphysicsCompMesh();
    
    /// Automatic builder for the computational mesh structure
    /// Build the multiphysics space using the previously provided mesh vector and active spaces
    virtual void AutoBuild() override;

    /// Automatic builder for the computational mesh structure
    /// Build the multiphysics space using the previously provided mesh vector and active spaces only
    /// for the specified set of matids
    virtual void AutoBuild(const std::set<int> &matids) override;

    /// Set active approximation spaces
    // active_approx_spaces : vector of the size of mesh_vector containing value 0 or 1
    void BuildMultiphysicsSpace(TPZVec<int> & active_approx_spaces, const TPZVec<TPZCompMesh * > & mesh_vector);

    /// @brief Build the multiphysics space using the previously provided mesh vector and active spaces
    virtual void BuildMultiphysicsSpace();

    /// Build the multiphysics space assuming all spaces are active
    virtual void BuildMultiphysicsSpace(TPZVec<TPZCompMesh * > & mesh_vector);
    
    /// Build the multiphysics elements only for the given geometric element indexes
    virtual void BuildMultiphysicsSpace(TPZVec<TPZCompMesh * > & mesh_vector, const TPZVec<int64_t> &gelindexes);

    /// Build the multiphysics elements only for the given geometric element indexes
    /// This signature assumes the mesh vector and active spaces have been already set
    virtual void BuildMultiphysicsSpace(const TPZVec<int64_t> &gelindexes);

    /// Build the multiphysics space where the elements are associated with a material with memory
    void BuildMultiphysicsSpaceWithMemory(TPZVec<int> & active_approx_spaces, TPZVec<TPZCompMesh * > & mesh_vector);

    /// Build the multiphysics space where the elements are associated with a material with memory
    /// This signature assumes the mesh vector and active spaces have been already set
    virtual void BuildMultiphysicsSpaceWithMemory();

    /// Build the multiphysics space where the elements are associated with a material with memory
    void BuildMultiphysicsSpaceWithMemory(TPZVec<int> & active_approx_spaces, TPZVec<TPZCompMesh * > & mesh_vector, std::set<int> matsIdWithMem, std::set<int> matsIdNoMem);

    /// Build the multiphysics space where the elements are associated with a material with memory
    /// This signature assumes the mesh vector and active spaces have been already set
    virtual void BuildMultiphysicsSpaceWithMemory(std::set<int> matsIdWithMem, std::set<int> matsIdNoMem);

    /// @brief Set the mesh vector and active approximation spaces
    void SetMeshVectorAndActiveSpaces(const TPZVec<TPZCompMesh * > & mesh_vector, const TPZVec<int> & active_approx_spaces) {
        if(NElements() > 0 || NConnects() > 0) {
            std::cout << "TPZMultiphysicsCompMesh::SetMeshVectorAndActiveSpaces. The multiphysics mesh already has elements. Clear the mesh before setting a new mesh vector and active spaces." << std::endl;
            DebugStop();
        }
        m_active_approx_spaces = active_approx_spaces;
        m_mesh_vector          = mesh_vector;
        if (m_mesh_vector.size() != m_active_approx_spaces.size()) {
            std::cout<< "TPZMultiphysicsCompMesh:: The vector provided should have the same size." << std::endl;
            DebugStop();
        }
        SetNMeshes(m_mesh_vector.size());
    }

    /// @brief Load the solution from the meshes in the multiphysics mesh
    virtual void LoadSolutionFromMeshes();
    
    /// @brief Load the solution from the multiphysics mesh to the meshes in the multiphysics mesh
    virtual void LoadSolutionFromMultiPhysics();

    /// Get the vector of computational meshes
    const TPZVec<TPZCompMesh *> &MeshVector() const {
        return  m_mesh_vector;
    }

    /// Get the vector of active physics
    const TPZVec<int> &GetActiveApproximationSpaces() const {
        return m_active_approx_spaces;
    }

    /// delete the elements and connects
    virtual void CleanElementsConnects();


protected:

/// @brief  populate the Referred vector with the computational elements of cmesh
/// @param cmesh Computational mesh from which the references will be loaded
/// @param Referred Vector to store the referred computational elements by geometric element index
static void LoadReferred(TPZCompMesh *cmesh, TPZVec<TPZCompEl *> &Referred);

    /// add the elements from the atomic meshes to the multiphysics elements
    virtual void AddElements();
    /// add the connects from the atomic meshes
    virtual void AddConnects();
    
    template<class TVar>
    void LoadSolutionFromMeshesInternal();
    template<class TVar>
    void LoadSolutionFromMultiPhysicsInternal();
};

#endif /* TPZMultiphysicsCompMesh_h */
