/**File : pzgclonemesh.h

Header file for class TPZGeoCloneMesh.
A TPZGeoCloneMesh defines a clone for geometrical mesh and contains a corresponding list of
	geometrical clone elements
	geometrical clone nodes
	coordinate systems
	elementwise defined boundary conditions ????
	nodal boundary conditions. ????
Various methods are defined to add, delete or loop over the items which are
contained within the TPZGeoMesh.
*/

#ifndef PZGEOCLONEMESHH
#define PZGEOCLONEMESHH

#include "pzgmesh.h"
#include "pzavlmap.h"
#include "pzadaptmesh.h"
#include "pzstack.h"

class TPZGeoNode;
class TPZGeoEl;



class  TPZGeoCloneMesh : public TPZGeoMesh
{

protected:
	/** Geometrical mesh associated*/
  	TPZGeoMesh 	*fGeoReference;

	/** Maps node id from cloned mesh to the mesh*/
	TPZAVLMap <int,int> fMapNodes;

	/** Maps element pointers from mesh to the cloned mesh*/
	TPZAVLMap <TPZGeoEl *,TPZGeoEl *> fMapElements;

	/** Maps element id from cloned mesh to mesh */
	TPZStack<TPZGeoEl *> fReferenceElement;

	/** Elements corresponding to the patch */
	TPZStack<TPZGeoEl *> fPatchReferenceElements;
	TPZStack<TPZGeoEl *> fPatchElements;

	/**
	 * Pointer to adapt mesh
	 * Its necessary to get element error vec
	 */
	TPZAdaptMesh * fAdaptMesh;

	/**
	 * Geometric Element Reference to get the clone mesh patch
	 */
	TPZGeoEl* fGeoElRef;    

	int IsPatchReferenceElement(TPZGeoEl *refpatch);
	int IsPatchElement(TPZGeoEl *refpatch);
	void AddBoundaryConditionElements(TPZGeoEl *eltoadd);

public:
  	/**Constructors and destructor*/
  	TPZGeoCloneMesh(TPZGeoMesh *mesh);

	/*Destructor*/
	virtual ~TPZGeoCloneMesh();

	/**
	 * Defines the mesh elements
	 * @param patch: elements to be cloned
	 */
	void SetElements(TPZStack<TPZGeoEl *> &patch, TPZGeoEl *ref);

	/**
	 * Add Elements to structure
	 */
	void AddElement(TPZGeoEl *eltoadd);

	/**
	 * Returns the reference element
	 */
	TPZGeoEl* ReferenceElement(int i);

	/**
	 * Only for debug purposes
	 */
	int NReference() {return fReferenceElement.NElements();}

	/**
	 * Print object data
	 */
	void Print (ostream & out);

	/**
	 * Return the index of an element
	 */
	int Index(TPZGeoEl *gel); 

	/**
	 * Return reference element index
	 */
	TPZGeoEl * GetMeshReferenceElement(){return fGeoElRef;}

	int IsPatchSon(TPZGeoEl *gel);

	/**
	 * Test and validation routines
	 */
	static int main();

protected:
	/**
	 * Creates a clone of a given node pointer
	 * return the id of the cloned node
	 * @param nod: node to be cloned
	 **/
	int CloneNode(int nodindex);

	/** 
	 * Creates an element clone and insert it into the mesh
	 * return the value of the cloned element
	 * @param org: geometric element to be cloned
	 **/
	int CloneElement(TPZGeoEl *orggel);

	/**
	 * Verifies if the specified node was created
	 * @param nodeindex: Node index to be verified
	 */
	int HasNode(int nodeindex);
	/**
	 * Verifies if a given element was created
	 * @param el: element to be verified
	 */
	int HasElement(TPZGeoEl *el);

private:
	/**
	 * Create a copy of the type of elementorgel
	 */
	TPZGeoEl * InitializeClone(TPZGeoEl* orgel);	
	
};


#endif



