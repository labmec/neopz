//$Id: pzgmesh.h,v 1.12 2004-04-05 20:23:13 longhin Exp $

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
#include <string>
#include <map>

#include "pzreal.h"
#include "pzeltype.h"
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
class TPZRefPattern;

template<class T>
class TPZVec;
//template<class T,class TT>
//class TPZAVLMap;

template <class TShape, class TGeo> class TPZGeoElRefPattern;

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

  typedef map<pair<int,int>, int> InterfaceMaterialsMap;
  
  InterfaceMaterialsMap fInterfaceMaterials;

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
  TPZCompMesh *Reference() const {return fReference;}

  /** Print the information of the grid to an ostream*/
virtual  void Print(ostream & out = cout);

  /**Returns the nearest node to the coordinate
     this method is VERY INEFICIENT*/
  TPZGeoNode* FindNode(TPZVec<REAL> &co);


  /**Alternative method for computing the connectivity*/
  void BuildConnectivity2();

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

  /**
   * Generic method for creating a geometric element. Putting this method centrally facilitates
   * the modification of the element type all through the code
   * @param type element topology
   * @param cornerindexes indexes of the corner nodes of the element
   * @param index index of the element in the vector of element pointers
   * @param reftype defines the type of refinement : 0 -> uniform 1-> refinement pattern
   */
  virtual  TPZGeoEl *CreateGeoElement(MElementType type,TPZVec<int> &cornerindexes,int matid,int &index, int reftype = 0);

//  virtual void DeleteElement(int gelindex);
  /**
   * Centralized method to delete elements
   * @param gel pointer to the element to be deleted
   * @param index index of the element
   */
  void DeleteElement(TPZGeoEl *gel, int index = -1);

  /**
   * Add an interface material associated to left and right element materials.
   * If pair<left, right> already exist, nothing is made and method returns 0.
   * If material is inserted in geomesh method returns 1.
   * @since Feb 05, 2004
   */
  int AddInterfaceMaterial(int leftmaterial, int rightmaterial, int interfacematerial);

  /**
   * Returns the interface material associated to left and right element materials.
   * @since Feb 05, 2004
   */
  int InterfaceMaterial(int leftmaterial, int rightmaterial);

  /**
   * Delete all interface materials in map.
   * @since Feb 05, 2004
   */
  void ClearInterfaceMaterialsMap();

 private:

  /**Find all elements in elmap or neighbour of elements in elmap which contain a node*/
  //void BuildElementsAroundNode(int currentnode,TPZAVLMap<int,TPZGeoEl *> &elmap);
 void BuildElementsAroundNode(int currentnode,map<int,TPZGeoEl *> &elmap);

  /**Method which works only for two dimensional topologies!
     Find, within elmap, the element which has currentnode as its first boundary side node*/
  //void FindElement(TPZAVLMap<int,TPZGeoEl *> &elmap,int currentnode,TPZGeoEl* &candidate,int &candidateside);
 	void FindElement(map<int,TPZGeoEl *> &elmap,int currentnode,TPZGeoEl* &candidate,int &candidateside);
protected: // Protected attributes
  /** Maps all refinement pattern objects in the mesh */
  map<int,map<string,TPZRefPattern *> > fRefPatterns;
public:
  void InsertRefPattern(TPZRefPattern *refpat);
  TPZRefPattern * GetRefPattern(int eltype, string name);
  /** Verifies if the side based refinement pattern exists. If the refinement pattern doesn't exists return a Null refinement Pattern. */
  TPZRefPattern * GetRefPattern (TPZGeoEl *gel, int side);
};

inline int TPZGeoMesh::AddInterfaceMaterial(int leftmaterial, int rightmaterial, int interfacematerial){
  pair<int, int> leftright(leftmaterial, rightmaterial);
  pair<int, int> rightleft(rightmaterial, leftmaterial);
  InterfaceMaterialsMap::iterator w, e;
  e = fInterfaceMaterials.end();

  w = fInterfaceMaterials.find(leftright);
  if (w == e) { //pair leftright does not exist yet
    w = fInterfaceMaterials.find(rightleft);
    if (w == e){ //pair rightleft does not exist too
      fInterfaceMaterials[leftright] = interfacematerial;
      return 1;
    }
  }
  return 0;
}

inline int TPZGeoMesh::InterfaceMaterial(int leftmaterial, int rightmaterial){
  pair<int, int> leftright(leftmaterial, rightmaterial);
  pair<int, int> rightleft(rightmaterial, leftmaterial);
  InterfaceMaterialsMap::iterator w, e;
  e = fInterfaceMaterials.end();
  w = fInterfaceMaterials.find(leftright);
  if (w != e)
    return w->second;

  w = fInterfaceMaterials.find(rightleft);
  if (w != e)
    return w->second;

  PZError << "\nTPZGeoMesh::InterfaceMaterial - Interface material not found " << endl;
  return -9999;
}

inline void TPZGeoMesh::ClearInterfaceMaterialsMap(){
  InterfaceMaterialsMap::iterator b, e;
  b = fInterfaceMaterials.begin();
  e = fInterfaceMaterials.end();
  fInterfaceMaterials.erase(b, e);    
}

#endif
