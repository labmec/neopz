//$Id: pzflowcmesh.h,v 1.10 2004-05-21 13:35:15 erick Exp $

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
using namespace std;
/*******************/
/**
 * This class implements the additions behavior and
 * data to handle CFD problems as computional meshes.
 *
 */
class TPZFlowCompMesh : public TPZCompMesh
{

public:

  //TPZFlowCompMesh();/**empty constructor*/

  TPZFlowCompMesh(TPZGeoMesh* gr);

  TPZFlowCompMesh();

  ~TPZFlowCompMesh(){};

  /**
   * According to each material settings, it returns
   * the greatest velocity in the mesh.
   */
  REAL MaxVelocityOfMesh();

  /**
   * Computes the current time step for the mesh.
   * User must call CollectFluidMaterials and
   * SetCFL first.
   */
  void ComputeTimeStep();

  /**
   * Informs all materials the CFL number.
   * User must call CollectFluidMaterials first.
   */
  void SetCFL(REAL CFL);

  /**
   * Scales the CFL of all materials
   *
   */
  void ScaleCFL(REAL scale);

  /**
   * Should be called after all materials have been
   * added to the mesh. This method collects all
   * fluid materials and stores a pointer to them
   * in the vector fFluidMaterial.
   */
  void CollectFluidMaterials();

   /**
    * Informs the time at which
    * the current solution in the computational
    * mesh belongs, so that the materials can
    * choose whether to contribute implicitly
    * or explicitly.
    */
  void SetContributionTime(TPZContributeTime time);

  /**
   * Sets the kind of residual to be computed
   */
  void SetResidualType(TPZResidualType type);

   /**
    * Sets the forcing funtion for all fluid materials in the mesh.
    */
  void SetFlowforcingFunction(void (*fp)(TPZVec<REAL> &loc,
					 TPZVec<REAL> &result));

  /**
   * Creates the computational elements, and the degree of freedom nodes.
   * In this reimplementation, also calls CollectFluidMaterials;
   */
  virtual void AutoBuild();

  /**
   * Returns the ith flow material in the mesh
   *
   */
  TPZMaterial * GetFlowMaterial(int i);

  /**
   * Returns the number of Flow materials.
   */
  int NFlowMaterials();


  /**
  * returns the unique identifier for reading/writing objects to streams
  */
  virtual int ClassId() const;
  /**
  Save the element data to a stream
  */
  virtual void Write(TPZStream &buf, int withclassid);

  /**
  Read the element data from a stream
  */
  virtual void Read(TPZStream &buf, void *context);



protected:

   /**
    * This vector of pointers represent the collection
    * of all fluid materials in the mesh. These are the
    * materials that deserve special attention during
    * the contribution processes.
    */
   TPZAdmChunkVector< TPZConservationLaw2 * > fFluidMaterial;
   

};

#endif

