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
#include "pzadaptmesh.h"
#include "pzstack.h"

#include <map>

class TPZGeoNode;
class TPZGeoEl;



class  TPZGeoCloneMesh : public TPZGeoMesh
{

protected:
	/** Geometrical mesh associated*/
  	TPZGeoMesh 	*fGeoReference;

	/** Maps node id from cloned mesh to the mesh*/
    std::map<int,int> fMapNodes;

	/** Maps element pointers from the original (reference) mesh to the cloned mesh
     */
    std::map<TPZGeoEl *,TPZGeoEl *> fMapElements;

	/** Maps element id from cloned mesh to mesh */
	TPZStack<TPZGeoEl *> fReferenceElement;

    /** Elements of the original (reference) mesh which defined the cloned mesh
     * fPatchReferenceElements is initialized in SetElements */
	TPZStack<TPZGeoEl *> fPatchReferenceElements;
    
    /** Elements of the cloned mesh corresponding to the patch defined in the original mesh */
	TPZStack<TPZGeoEl *> fPatchElements;

	/**
	 * Pointer to adapt mesh
	 * Its necessary to get element error vec
	 */
	TPZAdaptMesh * fAdaptMesh;

	/**
	 * Geometric Element int the cloned mesh corresponding to the root element in the original mesh
	 */
	TPZGeoEl* fGeoRoot;    

    /**
     * @brief returns true if the element belongs to the patch as defined in the original mesh (gel belongs to the original mesh)
     */
	int IsPatchReferenceElement(TPZGeoEl *gel);
    
    /**
     * @brief returns true if the element belongs to the patch elements in the cloned mesh */
	int IsPatchElement(TPZGeoEl *refpatch);
    
    /**
     * Given an element in the original mesh add all neighbours which correspond to boundary conditions to the clone
     */
	void AddBoundaryConditionElements(TPZGeoEl *eltoadd);

public:
  	/**Constructors and destructor*/
  	TPZGeoCloneMesh(TPZGeoMesh *mesh);

	/*Destructor*/
	virtual ~TPZGeoCloneMesh();

	/**
	 * Defines the mesh elements
	 * @param patch elements to be cloned
	 * @param ref element reference
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
	void Print (std::ostream & out);

	/**
	 * Return the index of an element
	 */
	int Index(TPZGeoEl *gel); 

	/**
	 * Return reference element index
	 */
	TPZGeoEl * GetMeshRootElement(){return fGeoRoot;}

	int IsPatchSon(TPZGeoEl *gel);

	/**
	 * Test and validation routines
	 */
	static int main();

protected:
	/**
	 * Creates a clone of a given node pointer
	 * return the id of the cloned node
	 * @param nodindex node to be cloned
	 **/
	int CloneNode(int nodindex);

	/** 
	 * Creates an element clone and insert it into the mesh
	 * return the value of the cloned element
	 * @param orggel geometric element to be cloned
	 **/
	int CloneElement(TPZGeoEl *orggel);

	/**
	 * Verifies if the specified node was created
	 * @param nodeindex Node index to be verified
	 */
	int HasNode(int nodeindex);
	/**
	 * Verifies if a given element was created
	 * @param el element to be verified
	 */
	int HasElement(TPZGeoEl *el);

private:
	/**
	 * Create a copy of the type of elementorgel
	 */
	TPZGeoEl * InitializeClone(TPZGeoEl* orgel);	
	
};


#endif



