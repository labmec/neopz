/**File : pzgclonemesh.h

Header file for class TPZGeoCloneMesh.
A TPZGeoCloneMesh defines a clone for geometrical mesh and contains a corresponding list of
	geometrical clone elements
	geometrical clone nodes
Various methods are defined to add, delete or loop over the items which are
contained within the TPZGeoMesh.
*/

#ifndef PZGEOCLONEMESHH
#define PZGEOCLONEMESHH

#include "pzgmesh.h"
#include "pzadaptmesh.h"
#include "pzstack.h"

#include <map>
#include <set>

class TPZGeoNode;
class TPZGeoEl;



class  TPZGeoCloneMesh : public TPZGeoMesh
{

protected:
	/** Geometrical mesh associated*/
  	TPZGeoMesh 	*fGeoReference;

	/** Maps node id from cloned mesh to the mesh*/
    std::map<int,int> fMapNodes;

	/** 
     * Maps element pointers from the original (reference) mesh to the cloned mesh
     */
    std::map<TPZGeoEl *,TPZGeoEl *> fMapElements;

	/** 
     * A vector of the size of the elements of the cloned mesh
     * Contains for each element in the cloned mesh a corresponding element in the reference mesh
     */
	TPZStack<TPZGeoEl *> fReferenceElement;
    
    /** 
     * Elements of the cloned mesh corresponding defining the patch
     */
    std::set<TPZGeoEl *> fPatchElements;

	/**
	 * Pointer to adapt mesh
	 * Its necessary to get element error vec
	 */
	TPZAdaptMesh * fAdaptMesh;

	/**
	 * Geometric Element int the cloned mesh corresponding to the root element in the original mesh
     * The root element is the element which generates the patch
     * The root elements form a partition of the original mesh
	 */
	TPZGeoEl* fGeoRoot;    

    /**
     * @brief returns true if the element belongs to the patch elements in the cloned mesh */
//	int IsPatchElement(TPZGeoEl *refpatch);
    
    /**
     * Given an element in the original mesh add all neighbours which are boundary condition elements
     */
	void AddBoundaryConditionElements(TPZGeoEl *eltoadd);

public:
  	/**Constructors and destructor*/
  	TPZGeoCloneMesh(TPZGeoMesh *mesh);

	/*Destructor*/
	virtual ~TPZGeoCloneMesh();

	/**
	 * Defines the mesh elements
	 * @param patch elements to be cloned - these are pointers in the original mesh
	 * @param ref element reference - element which generated the patch
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
	
	std::set<TPZGeoEl *> &PatchElements()
	{
		return fPatchElements;
	}
    /**
     * returns true if the element is a sibling of the set of patch elements
     */
	int IsPatchSon(TPZGeoEl *gel) const;
	int IsNeighBCForPatchSon(TPZGeoEl *gel) const;

	/**
	 * Test and validation routines
	 */
	static int main();

protected:
	/**
	 * Creates a clone of a given node defined by its index
	 * return the index of the cloned node
     * this method update fMapNodes
	 * @param nodindex node to be cloned
	 **/
	int CloneNode(int nodindex);

	/** 
	 * Creates an element clone and insert it into the mesh
     * updates fMapElements
     * updates fReferenceElements
	 * return the index of the cloned element
	 * @param orggel geometric element to be cloned
	 **/
	int CloneElement(TPZGeoEl *orggel);

	/**
	 * Verifies if the specified nodeindex was cloned
     * this method uses the datastructure fMapNodes
	 * @param nodeindex Node index in the original mesh
	 */
	int HasNode(int nodeindex);
	/**
	 * Verifies if a given element was cloned
     * this method uses the datastructure fMapElements
	 * @param el element pointer in the original mesh
	 */
	int HasElement(TPZGeoEl *el);

private:
	/**
	 * Create a copy of the type of elementorgel
	 */
	TPZGeoEl * InitializeClone(TPZGeoEl* orgel);	
	
};


#endif



