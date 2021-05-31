/**
 * @file
 * @brief Contains the declaration of the TPZBuildmultiphysicsMesh class.
 */

#ifndef PZBUILDMULTIPHYSICSMESHH
#define PZBUILDMULTIPHYSICSMESHH

#include <iostream>

#include "pzcompel.h"
#include "pzmultiphysicselement.h"
class TPZLinearAnalysis;

/**
 * @brief This class has methods to build the mesh multiphysics
 * @author Agnaldo
 * @since 10/31/2011
 */

typedef std::pair<TPZCompMesh *, int64_t> atomic_index;

class TPZBuildMultiphysicsMesh {
	
public:
	TPZBuildMultiphysicsMesh();
	
	~TPZBuildMultiphysicsMesh();
	
	/**
	 * @brief Creating multiphysic elements into mphysics computational mesh. Method to add elements in the mesh multiphysics
	 * @param cmeshVec [in]:  pointer to an vector of meshes
	 * @param MFMesh [out]: my mesh multiphysics  
	 */
	static void AddElements(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);
	
	/*
	 * @brief Method to add connects in the mesh multiphisics
	 * @param cmeshVec [in]: pointer to an vector of meshes
	 * @param MFMesh [out]: my mesh multiphysics  
	 */
	static void AddConnects(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);
    
    /*
	 * @brief Methodo to append connects from mesh to a multiphysics mesh
	 * @param cmesh: pointer to an computational mesh
	 * @param MFMesh [out]: my mesh multiphysics
	 */
    static void AppendConnects(TPZCompMesh *cmesh, TPZCompMesh *MFMesh);
	
	/**
	 * @brief Transfer information from a specific set of meshes for the current mesh multiphysics
	 * @param cmeshVec [in] vector of meshes. Transfers the information
	 * @param MFMesh [out] mesh pointer who will receive the information
	 */	
	static void TransferFromMeshes(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);
	
    /**
     * @brief Transfer information from a specific set of meshes for the current mesh multiphysics
     * @param cmeshVec [in] vector of meshes. Transfers the information
     * @param MFMesh [out] mesh pointer who will receive the information
     */
    static void TransferFromMeshes(TPZVec<TPZAutoPointer<TPZCompMesh> > &cmeshVec, TPZAutoPointer<TPZCompMesh> MFMesh)
    {
        TPZManVector<TPZCompMesh *,4> meshVecPtr(cmeshVec.size());
        for (int i=0; i<cmeshVec.size(); i++) {
            meshVecPtr[i] = cmeshVec[i].operator->();
        }
        TransferFromMeshes(meshVecPtr, MFMesh.operator->());
    }
    
	/**
	 * @brief Transfer information from a specific mesh multiphysics for the current specific set of meshes 
	 * @param cmeshVec [out] vector of meshes that will receive the information.
	 * @param MFMesh [in] mesh pointer that Transfers the information 
	 */	
	static void TransferFromMultiPhysics(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);
		
    /**
     * @brief Transfer information from a specific mesh multiphysics for the current specific set of meshes
     * @param cmeshVec [out] vector of meshes that will receive the information.
     * @param MFMesh [in] mesh pointer that Transfers the information
     */
    static void TransferFromMultiPhysics(TPZVec<TPZAutoPointer<TPZCompMesh> > &cmeshVec, TPZAutoPointer<TPZCompMesh> MFMesh)
    {
        TPZManVector<TPZCompMesh *,4> meshVecPtr(cmeshVec.size());
        for (int i=0; i<cmeshVec.size(); i++) {
            meshVecPtr[i] = cmeshVec[i].operator->();
        }
        TransferFromMultiPhysics(meshVecPtr, MFMesh.operator->());
    }
    
    
    /**
     * @brief Show shape functions associated with connects of a multiphysics mesh
     */
    static void ShowShape(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh, TPZLinearAnalysis &analysis, const std::string &filename, TPZVec<int64_t> &equationindices);
	/**
	 * @brief Creating computational mesh with interface elements
	 * @param cmesh [in]: computational mesh
	 * @param MaterialIDs [in]: set of index materials
	 * @param LagrangeMat [in]: index Lagrange material
	 * @param InterfaceMat [in]: index interface material
	 */	
	static void BuildHybridMesh(TPZCompMesh *cmesh, std::set<int> &MaterialIDs, std::set<int> &BCMaterialIds, int LagrangeMat, int InterfaceMat);
	
	/**
	 * @brief Uniform refinement of the computational mesh
	 * @param cMesh [in]: computational mesh
	 * @param ndiv [in]: number of refinements
	 */
	static void UniformRefineCompMesh(TPZCompMesh  *cMesh, int ndiv, bool isLagrMult=true);
	
	/**
	 * @brief Uniform refinement of the computational element
	 * @param cMesh [in] : computational mesh
	 * @param indexEl [in]: index of the element
	 */
	static void UniformRefineCompEl(TPZCompMesh  *cMesh, int64_t indexEl, bool isLagrMult);
    
    /**
     *@brief Create skeleton elements of the wrap of me.
     *@param mfcel [in] : multifysics element
     *@param matskeleton [in]: Material Id for skeleton elements
     *@param ListGroupEl [out]: List of me and my elements of wrap
     */
    static void AddWrap(TPZMultiphysicsElement *mfcel, int matskeleton, TPZStack< TPZStack<TPZMultiphysicsElement *,7> > &ListGroupEl);
    
    /**
     * Compute the correspondence between the connect index in the multiphysics
     * mesh and the connect indexes in the atomic meshes
     */
    static void ComputeAtomicIndexes(TPZCompMesh *mesh, TPZVec<atomic_index> &indexes);
protected:
    template<class TVar>
    static void TransferFromMeshesT(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);
    template<class TVar>
    static void TransferFromMultiPhysicsT(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);
};

#endif
