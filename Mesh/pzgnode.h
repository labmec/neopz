/**File : pzgnode.h

Header file for class TPZGeoNode. TPZGeoNode defines a geometrical node.

*/

#ifndef  TPZGEONODEH
#define  TPZGEONODEH

#include <iostream>
#include "pzreal.h"

using namespace std;

template<class T>
class TPZVec;

class TPZGeoMesh;
class TPZGeoEl;

class TPZGeoNode {

  /**Identity of node*/
  int		fId;
  /**Node coordinates*/
  REAL	fCoord[3];

 public:
  /**Constructor to new node with id predefined*/
  TPZGeoNode(int id,TPZVec<REAL> &xp,TPZGeoMesh &mesh);
  /**Constructor to new node */
  TPZGeoNode();
  /**Constructor copy*/
  TPZGeoNode(TPZGeoNode &node);
  /**Destructor*/
  ~TPZGeoNode() { }

  /**Return the identity of the current node*/
  int Id() { return fId; }

  void SetNodeId(int id) { fId = id;}

  /**Initialize the data structure of the node. Creates a unique id for the node automatically*/
  void Initialize(TPZVec<REAL> &coord,TPZGeoMesh &mesh);
  /**Initialize the data structure of the node. Assumes that the id provided by the user
     is unique for the mesh*/
  void Initialize(int id,TPZVec<REAL> &coord,TPZGeoMesh &mesh);

  /**Return i-th coordinate of the current node*/
  REAL Coord(int i);

  /**Set all coordinates into the current node. It gets the dim values from x */
  void SetCoord(REAL *x,int dim = 3);
  /**Set the i-th coordinate for current node*/
  void SetCoord(int i,REAL coord);

  void Print(ostream & out = cout);
};


/**TPZGeoNodeBC

TPZGeoNodeBc defines a boundary condition applied to a geometrical node
*******       *******/

struct TPZGeoNodeBC {
  TPZGeoNode	 	*fNode;
  int			fBCId;
  TPZGeoEl *fGeoEl;
  int fGeoElSide;

  TPZGeoNodeBC() {
    fNode = 0;
    fBCId = 0;
    fGeoEl = 0;
    fGeoElSide = -1;
  }

  TPZGeoNodeBC(TPZGeoNode *node,int bcid, TPZGeoEl *gel, int gelside) {
    fNode = node;
    fBCId = bcid;
    fGeoEl = gel;
    fGeoElSide = gelside;
  }
};

#endif

