/**
 * @file
 * @brief Contains declaration of TPZFlowCompMesh class which is a computational mesh with additional data for CFD problems.
 */
//$Id: pzflowcmesh.h,v 1.15 2007-01-03 00:06:47 phil Exp $

#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzconslaw.h"
#include "pzerror.h"

#ifndef TPZFLOWCOMPMESH_H
#define TPZFLOWCOMPMESH_H

class TPZMaterial;
class TPZGeoMesh;
#include <iostream>
#include "pzcmesh.h"

/**
 * @brief Computational mesh with additional data for CFD problems. \ref CompMesh "Computational Mesh"
 * @ingroup CompMesh
 */
/**
 * This class implements the additions behavior and
 * data to handle CFD problems as computional meshes.
 */
class TPZFlowCompMesh : public TPZCompMesh
{
	
public:
	
	//TPZFlowCompMesh();/**empty constructor*/
	
	TPZFlowCompMesh(TPZGeoMesh* gr);
	
	TPZFlowCompMesh();
	
	~TPZFlowCompMesh(){};
	
	/**
	 * @brief According to each material settings, it returns
	 * the greatest velocity in the mesh.
	 */
	REAL MaxVelocityOfMesh();
	
	/**
	 * @brief Computes the current time step for the mesh.
	 * 
	 * User must call CollectFluidMaterials and
	 * SetCFL first.
	 */
	REAL ComputeTimeStep();
	
	/**
	 * @brief Informs all materials the CFL number.
	 * 
	 * User must call CollectFluidMaterials first.
	 */
	void SetCFL(REAL CFL);
	
	/**
	 * @brief Scales the CFL of all materials
	 *
	 */
	void ScaleCFL(REAL scale);
	
	/**
	 * @brief Should be called after all materials have been
	 * added to the mesh.
	 *
	 * This method collects all
	 * fluid materials and stores a pointer to them
	 * in the vector fFluidMaterial.
	 */
	void CollectFluidMaterials();
	
	/**
	 * @brief Informs the time at which
	 * the current solution in the computational
	 * mesh belongs, so that the materials can
	 * choose whether to contribute implicitly
	 * or explicitly.
	 */
	void SetContributionTime(TPZContributeTime time);
	
	/**
	 * @brief Sets the kind of residual to be computed
	 */
	void SetResidualType(TPZResidualType type);
	
	/**
	 * @brief Sets the forcing funtion for all fluid materials in the mesh.
	 */
	void SetFlowforcingFunction(TPZAutoPointer<TPZFunction> fp);
	
	/**
	 * @brief Creates the computational elements, and the degree of freedom nodes.
	 * In this reimplementation, also calls CollectFluidMaterials;
	 */
	virtual void AutoBuild();
	
	/**
	 * @brief Returns the first flow material in the mesh
	 *
	 */
	TPZAutoPointer<TPZMaterial> GetFlowMaterial();
	
	/**
	 * @brief Returns the number of Flow materials.
	 */
	int NFlowMaterials();
	
	
	/**
	 * @brief Returns the unique identifier for reading/writing objects to streams
	 */
	virtual int ClassId() const;
	/**
	 @brief Saves the element data to a stream
	 */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/**
	 @brief Read the element data from a stream
	 */
	virtual void Read(TPZStream &buf, void *context);
	
	/**
	 * @brief Adapt the solution vector to new block dimensions
	 */
	virtual void ExpandSolution2();
	
protected:
	
	/**
	 * @brief This vector of pointers represent the collection
	 * of all fluid materials in the mesh. 
	 *
	 * These are the
	 * materials that deserve special attention during
	 * the contribution processes.
	 */
	std::map<int, TPZAutoPointer< TPZMaterial> > fFluidMaterial;
	
	
};

#endif
