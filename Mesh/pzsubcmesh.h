// -*- c++ -*-
// subcmesh.h: interface for the TPZSubCompMesh class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SUBCMESH_H
#define SUBCMESH_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <stdlib.h>
#include "pzcompel.h"
#include "pzcmesh.h"
#include "pzsmanal.h"
#include "pzvec.h"
#include "pzreal.h"

class TPZSubMeshFrontalAnalysis;
/**
 * Class TPZSubCompMesh derived from Computacional mesh and
 * computacional element classes.
 * @ingroup CompMesh
 * @ingroup CompElement
 */

class TPZSubCompMesh : 
	public TPZCompMesh,
	public TPZCompEl
{
protected:
  /**
   * Pointer to submesh analysis object. Defines the resolution type.
   */
  TPZSubMeshFrontalAnalysis *fAnalysis;

  /**
   * Pointer to external location index of the connection. 
   * If the connection hasn't external location return the local id. 
   */
  TPZVec<int> fConnectIndex;

  /**
   * Indexes of the external connections.
   * If the connection isn't external id is -1! 
   */
  TPZVec<int> fExternalLocIndex;


private:
  /**
   * Transfer one element from a submesh to another mesh. 
   */
  int TransferElementTo(TPZCompMesh * mesh, int elindex);
  
  /**
   * Transfer one element from a specified mesh to the current submesh. 
   */
  int TransferElementFrom(TPZCompMesh *mesh, int elindex);
  

public:
	TPZAnalysis * GetAnalysis();
  /**
   * Constructor. 
   * @param mesh reference mesh
   * @param index reference mesh element index to transfer to submesh
   */
  TPZSubCompMesh(TPZCompMesh &mesh, int &index);

  /**
   * Destructor. 
   */
  virtual ~TPZSubCompMesh();

  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
    cout << "TPZSubCompMesh::Clone should be implemented\n";
    return 0;
  }

  /**
   * Static function for validation tests.
   */
  static int main();

  /**
   * Set the analysis type.
   */
  void SetAnalysis();

  /**
   * This method will load the elements of the grid in their corresponding geometric
   * elements
   **/
  virtual void LoadElementReference();

  /**
   * This method will initiate the comparison between the current computational
   * mesh and the mesh which is referenced by the geometric mesh
   * @param var state variable number
   * @param matname pointer to material name
   **/
  virtual REAL CompareElement(int var, char *matname);	

  /**
   * Verifies the transfer possibility of the connection elindex from
   * the mesh to the current submesh.
   * @param mesh  pointer to given mesh
   * @param elindex given mesh element index
   */
  int IsAllowedElement(TPZCompMesh *mesh, int elindex);

 /**
   * @name Mesh
   * Methods derived from TPZCompMesh
   */

  //@{	
  // ==== Virtual methods - Computacional Mesh derived ====

  /**
   * Transfer one element form a submesh to another mesh. 
   */
  virtual int TransferElement(TPZCompMesh *mesh, int elindex);

  /**
   * Make all mesh connections internal mesh connections. 
   */
  virtual void MakeAllInternal();

  /** 
   * compute the number of internal equations
   */
  int NumInternalEquations();
	
  /**
   * This method computes the skyline of the system of equations
   * @param skyline vector where the skyline will be computed
   */
  virtual void Skyline(TPZVec<int> &skyline);

  /**
   * Returns the rootmesh who have the specified connection. 
   * @param local connection local index
   */
  virtual TPZCompMesh * RootMesh(int local);

  /**
   * Makes a specified connection a internal mesh connection. 
   * @param local connection local number to be processed
   */
  virtual void MakeInternal(int local);

   /**
    * Put an local connection in the supermesh - Supermesh is one
    * mesh who contains the analised submesh. 
    * @param local local index of the element to be trasnfered
    * @param super pointer to the destination mesh
    */
  virtual int PutinSuperMesh(int local, TPZCompMesh *super);

  /**
   * Get an external connection from the supermesh - Supermesh is one
   * mesh who contains the analised submesh. 
   * @param superind index of the element to be trasnfered
   * @param super pointer to the destination mesh
   */
  virtual int GetFromSuperMesh(int superind, TPZCompMesh *super);

  /**
   * Changes an local internal connection to a external connection in the father mesh. 
   * @param local makes the connect with index local an external node
   */
  void MakeExternal(int local);

  /**
   * Gives the commom father mesh of the specified mesh and the current submesh. 
   * @param mesh pointer to other mesh whosw want to know the commom father mesh
   */
  virtual TPZCompMesh * CommonMesh (TPZCompMesh *mesh);

  /**
   * Return the current submesh father mesh . 
   */
  virtual TPZCompMesh * FatherMesh();
  //@}

  /**
   * Optimize the connections positions on block.
   * void TPZSubCompMesh::PermuteExternalNodes(){
   */
  void PermuteExternalConnects();

  /**
   * Print the submesh information on the specified device/file out.
   * This method use the virtual method from Computacional Mesh class.
   * @param out indicates the device where the data will be printed
   */
  void Prints(ostream &out = cout);

  /**
   * @name Element
   * Virtual methods derived from TPZCompEl
   */
  //@{
  //  /**
  //     * Changes the current node index -inode- to the specified node- index. 
  //     * 
  //     */
	virtual void SetConnectIndex(int inode, int index);

  //	/**
  //     * Calculates the submesh stiffness matrix
  //     */
	virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);
    
  //	/**
  //     * Virtual Method! 
  //     */
	virtual int AllocateNewConnect(int blocksize);

  //    /**
  //     * Gives the id node  of one local node in a neighbour mesh. 
  //     */
	int NodeIndex(int nolocal, TPZCompMesh *neighbour);

  //    /**
  //     * Virtual Method! See declaration in the TPZCompEl class. 
  //     */
	virtual void SetMaterial(TPZMaterial *mat);

  //    /**
  //     * Virtual Method! See declaration in TPZCompEl class. 
  //     */
	virtual TPZMaterial *Material() const;

  //    /**
  //     * Virtual Method! See declaration in TPZCompEl class. 
  //     * The use of this method in submesh class return -1 == Error! 
  //     */
	virtual int Dimension() const;

  //    /**
  //     * Returns the connection index i. 
  //     */
	virtual int ConnectIndex(int i) const;

  //    /**
  //     * Returns the number of connections. 
  //     */
	virtual int NConnects() const;
	
  //  /**
  //     * Load the father mesh solution to all submesh connects -
  //	 * (internal and external).
  //     */
	virtual void LoadSolution();
	
	virtual void GetExternalConnectIndex (TPZVec<int> &extconn);
  //@}
};

#endif
