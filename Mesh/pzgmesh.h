/**File : pzgmes.h

Header file for class TPZGeoMesh.
A TPZGeoMesh defines a geometrical mesh and contains a corresponding list of
	geometrical elements
	geometrical nodes
	coordinate systems
	elementwise defined boundary conditions
	nodal boundary conditions.
Various methods are defined to add, delete or loop over the items which are
contained within the TPZGeoMesh.
*/

#ifndef PZGEOMESHH
#define PZGEOMESHH

#ifndef ISZERO
#define ISZERO
#define IsZero( x )    ( (x) < 1.e-10 && (x) > -1.e-10 )
#endif


#include <iostream>
#include <string.h>

#include "pzreal.h"
#include "pzgnode.h"
#include "pzbndcond.h"
#include "pzcmesh.h"
#include "pzadmchunk.h"

class TPZMaterial;
class TPZGeoNode;
struct TPZGeoNodeBC;
class TPZGeoEl;
struct TPZGeoElBC;
class TPZCosys;
class TPZMatrix;
class TPZCompMesh;

template<class T>
class TPZVec;
template<class T,class TT>
class TPZAVLMap;


class  TPZGeoMesh {

  /** TPZGeoMesh name for model identification*/
  char			fName[63];
  /** Computational mesh associated*/
  TPZCompMesh 	*fReference;

  /** List of pointers to finite elements*/
  TPZAdmChunkVector<TPZGeoEl *>		fElementVec;
  /** List of pointers to nodes*/
  TPZAdmChunkVector<TPZGeoNode>		fNodeVec;
  /** List of pointers to reference systems*/
  TPZAdmChunkVector<TPZCosys *>		fCosysVec;
  /** List of pointers to elements, boundary side and type*/
  TPZAdmChunkVector<TPZGeoElBC>		fBCElementVec;
  /** List of pointers to nodes and bc-type*/
  TPZAdmChunkVector<TPZGeoNodeBC>	fBCNodeVec;

  /**Maximum id used by all nodes of this mesh*/
  int fNodeMaxId;
  /**Maximum id used by all elements of this mesh*/
  int fElementMaxId;

 public:
  /**Constructors and destructor*/
  TPZGeoMesh();
  /*Destructor*/
  virtual ~TPZGeoMesh();

  /**Deletes all items in the TPZGeoMesh*/
  void CleanUp();

  /**Indicates that a node with id was created*/
  void SetNodeIdUsed(int id) { fNodeMaxId = (id > fNodeMaxId) ? id : fNodeMaxId; }
  /**Indicates that an element with id was created*/
  void SetElementIdUsed(int id) { fElementMaxId = (id > fElementMaxId) ? id : fElementMaxId; }

  /**Return ++fNodeMaxId*/
  int CreateUniqueNodeId() { return ++fNodeMaxId; }
  /**Return ++fElementMaxId*/
  int CreateUniqueElementId() { return ++fElementMaxId; }

  /**Number of nodes of the mesh*/
  int NNodes() {return fNodeVec.NElements();}
  /**Number of elements of the mesh*/
  int NElements() {return fElementVec.NElements();}

  void SetName(char *name);
  char* Name() { return fName; }

  /**Methods for handling pzlists*/
  TPZAdmChunkVector<TPZGeoEl *> &ElementVec() { return fElementVec; }
  TPZAdmChunkVector<TPZGeoNode> &NodeVec() { return fNodeVec; }
  TPZAdmChunkVector<TPZGeoElBC> &BCElementVec() { return fBCElementVec; }
  TPZAdmChunkVector<TPZGeoNodeBC> &BCNodeVec() { return fBCNodeVec; }
  TPZAdmChunkVector<TPZCosys *> &CosysVec() { return fCosysVec; }

  /**Resets all load references in elements and nodes*/
  void ResetReference();
  /**Restore all reference in elements from computational mesh criated from current
     geometrical mesh previously*/
  void RestoreReference(TPZCompMesh *cmesh);
  /**Sets the reference of the geometric grid to ref*/
  void SetReference(TPZCompMesh *ref) { fReference = ref;}
  /**Returns the currently loaded computational grid*/
  TPZCompMesh *Reference() {return fReference;}

  /** Print the information of the grid to an ostream*/
virtual  void Print(ostream & out = cout);

  /**Returns the nearest node to the coordinate
     this method is VERY INEFICIENT*/
  TPZGeoNode* FindNode(TPZVec<REAL> &co);

  /**Build the connectivity of the grid*/
  void BuildConnectivity();

  /**Fills the nodep vector with pointers to the nodes identified by their indexes*/
  void GetNodePtr(TPZVec<int> &nos,TPZVec<TPZGeoNode *> &nodep);

  /**GetBoundaryElements returns all elements beweeen NodFrom and NodTo counterclock wise this
     method uses the connectivity of the elements BuildConnectivity should be called to initialize
     the connectivity information this method will only work for grid with 2-D topology the current
     version will only work for a grid with only one level*/
  void GetBoundaryElements(int IndexNodeFrom,int IndexNodeTo,TPZStack<TPZGeoEl *> &ElementVec,TPZStack<int> &Sides);

    /**Find the element with ids elid*/
  TPZGeoEl *FindElement(int elid);//Cedric

  /**Returns the index of the given element into the fElementVec*/
  int ElementIndex(TPZGeoEl * gel); //Cesar 2002-05-02
  
   /**Returns the index of the given node into the fNodeVec*/
  int NodeIndex(TPZGeoNode * nod);//Cesar 2002-05-02



 private:

  /**Find all elements in elmap or neighbour of elements in elmap which contain a node*/
  void BuildElementsAroundNode(int currentnode,TPZAVLMap<int,TPZGeoEl *> &elmap);

  /**Method which works only for two dimensional topologies!
     Find, within elmap, the element which has currentnode as its first boundary side node*/
  void FindElement(TPZAVLMap<int,TPZGeoEl *> &elmap,int currentnode,TPZGeoEl* &candidate,int &candidateside);
};


#endif



