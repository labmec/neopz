/**
 * @file
 * @brief Contains the declaration of the TPZBuildmultiphysicsMesh class.
 */

#ifndef PZBUILDMULTIPHYSICSMESHH
#define PZBUILDMULTIPHYSICSMESHH

#include <iostream>

#include "pzcompel.h"

/**
 * @brief This class has methods to build the mesh multiphysics
 * @author Agnaldo
 * @since 10/31/2011
 */

class TPZBuildMultiphysicsMesh{
	
	
public:
	TPZBuildMultiphysicsMesh();
	
	~TPZBuildMultiphysicsMesh();
	
	/*
	 *@brief Creating multiphysic elements into mphysics computational mesh.
	 Method to add elements in the mesh multiphysics
	 *@param cmeshVec [in]:  pointer to an vector of meshes
	 @param MFMesh [out]: my mesh multiphysics  
	 */
	static void AddElements(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh);
	
	/*
	 *@brief Method to add connects in the mesh multiphisics
	 *@param cmeshVec [in]: pointer to an vector of meshes
	 @param MFMesh [out]: my mesh multiphysics  
	 */
	static void AddConnects(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh);
	
	/**
	 * @brief Transfer information from a specific set of meshes for the current mesh multiphysics
	 * @param cmeshVec [in] vector of meshes. Transfers the information
	 * @param MFMesh [out] mesh pointer who will receive the information
	 */	
	static void TransferFromMeshes(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);
	
	/**
	 * @brief Transfer information from a specific mesh multiphysics for the current specific set of meshes 
	 * @param cmeshVec [out] vector of meshes that will receive the information.
	 * @param MFMesh [in] mesh pointer that Transfers the information 
	 */	
	static void TransferFromMultiPhysics(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);
		
	/**
	 * @brief Creating computational mesh with interface elements
	 * @param cmesh [in]: computational mesh
	 * @param MaterialIDs [in]: set of index materials
	 * @param LagrangeMat [in]: iindex Lagrange material
	 * @param InterfaceMat [in]: iindex interface material
	 */	
	static void BuildHybridMesh(TPZCompMesh *cmesh, std::set<int> &MaterialIDs, int LagrangeMat, int InterfaceMat);
	
	/**
	 * brief Uniform refinement of the computational mesh
	 * @param cMesh [in]: computational mesh
	 * @param ndiv [in]: number of refinements
	 */
	static void UniformRefineCompMesh(TPZCompMesh  *cMesh, int ndiv);
	
	/**
	 * brief Uniform refinement of the computational element
	 * @param cMesh [in] : computational mesh
	 * @param indexEl [in]: index of the element
	 */
	static void UniformRefineCompEl(TPZCompMesh  *cMesh, int indexEl);
		
};

#endif
