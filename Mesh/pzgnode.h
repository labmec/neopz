//$Id: pzgnode.h,v 1.8 2005-02-28 22:08:52 phil Exp $

/**File : pzgnode.h

Header file for class TPZGeoNode. TPZGeoNode defines a geometrical node.

*/

#ifndef  TPZGEONODEH
#define  TPZGEONODEH

#include <iostream>
#include "pzreal.h"
#include "pzerror.h"
#include "pzsave.h"
#include "pzstream.h"
#include "pzmeshid.h"

using namespace std;

template<class T>
class TPZVec;

class TPZGeoMesh;
class TPZGeoEl;

/// Implements a geometric node in the pz environment
/**
A geometric node is a place holder for 3d coordinates and a global Id
Note that the global id will influence the orientation of the shape functions
It is very important that the global Ids will not be duplicated within a single mesh
@ingroup geometry
*/
class TPZGeoNode : public TPZSaveable {

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
virtual  ~TPZGeoNode() { }
  
  /**
  return the id of the class (used for writing reading the object)
  */
  virtual int ClassId() const {
    return TPZGEONODEID;
  }
  
  /**
  Read the object from disk
  */
  virtual void Read(TPZStream &buf, void *context) {
    TPZSaveable::Read(buf,context);
    buf.Read(&fId,1);
    buf.Read(fCoord,3);
  }
  
  /**
  Write the object to disk
  */
  virtual void Write(TPZStream &buf, int withclassid) {
    TPZSaveable::Write(buf,withclassid);
    buf.Write(&fId,1);
    buf.Write(fCoord,3);
  }

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

template class TPZRestoreClass<TPZGeoNode,TPZGEONODEID>;

inline REAL TPZGeoNode::Coord(int i) {
#ifndef NODEBUG
  if(i > 2 || i < 0) {
    PZError << "Not exist (TPZGeoNode) coordinate " << i << endl;
    return 1.e12;
  }
#endif
  return fCoord[i];
}

///TPZGeoNodeBc defines a boundary condition applied to a geometrical node
/**
\deprecated
@ingroup geometry
*/

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

