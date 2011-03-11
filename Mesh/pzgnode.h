//$Id: pzgnode.h,v 1.15 2011-03-11 13:27:49 fortiago Exp $

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
  TPZGeoNode(const TPZGeoNode &node);

  TPZGeoNode & operator=(const TPZGeoNode &node);

  /**Destructor*/
virtual  ~TPZGeoNode() { }
  
  /**
  return the id of the class (used for writing reading the object)
  */
  virtual int ClassId() const;
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
  int Id() const
 {
	return fId; 
 }

  void SetNodeId(int id) { fId = id;}

  /**Initialize the data structure of the node. Creates a unique id for the node automatically*/
  void Initialize(TPZVec<REAL> &coord,TPZGeoMesh &mesh);
  /**Initialize the data structure of the node. Assumes that the id provided by the user
     is unique for the mesh*/
  void Initialize(int id,TPZVec<REAL> &coord,TPZGeoMesh &mesh);

  /**Initialize the node with data from a node from a different mesh */
  void Initialize(const TPZGeoNode &node, TPZGeoMesh &mesh);  
  /**Return i-th coordinate of the current node*/
  REAL Coord(int i) const;

  /**Set all coordinates into the current node. It gets the dim values from x */
  void SetCoord(const TPZVec<REAL> &x);


  /**Set the i-th coordinate for current node*/
  void SetCoord(int i,REAL coord);
	
	/**
	 * fill the coordinates of the node
	 */
	void GetCoordinates(TPZVec<REAL> &co);

  void Print(std::ostream & out = std::cout);
};


inline REAL TPZGeoNode::Coord(int i) const {
#ifndef NODEBUG
  if(i > 2 || i < 0) {
    PZError << "Not exist (TPZGeoNode) coordinate " << i << std::endl;
    return 1.e12;
  }
#endif
  return fCoord[i];
}

#endif

