/**
 * @file
 * @brief Contains declaration of TPZGeoNode class which defines a geometrical node.
 */
//$Id: pzgnode.h,v 1.15 2011-03-11 13:27:49 fortiago Exp $

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

/**
 * @brief Implements a geometric node in the pz environment. \ref geometry "Geometry"
 * @ingroup geometry
 */
/**
 * A geometric node is a place holder for 3d coordinates and a global Id \n
 * Note that the global id will influence the orientation of the shape functions \n
 * It is very important that the global Ids will not be duplicated within a single mesh
 */
class TPZGeoNode : public TPZSaveable {
	
	/** @brief Identity of node*/
	int		fId;
	/** @brief Node coordinates*/
	REAL	fCoord[3];
	
public:
	/** @brief Constructor to new node with id predefined*/
	TPZGeoNode(int id,TPZVec<REAL> &xp,TPZGeoMesh &mesh);
	/** @brief Constructor to new node */
	TPZGeoNode();
	/** @brief Constructor copy*/
	TPZGeoNode(const TPZGeoNode &node);
	
	TPZGeoNode & operator=(const TPZGeoNode &node);
	
	/** @brief Destructor*/
	virtual  ~TPZGeoNode() { }
	
	/**
	 @brief Returns the id of the class (used for writing reading the object)
	 */
	virtual int ClassId() const;
	/**
	 @brief Reads the object from disk
	 */
	virtual void Read(TPZStream &buf, void *context) {
		TPZSaveable::Read(buf,context);
		buf.Read(&fId,1);
		buf.Read(fCoord,3);
	}
	
	/**
	 @brief Writes the object to disk
	 */
	virtual void Write(TPZStream &buf, int withclassid) {
		TPZSaveable::Write(buf,withclassid);
		buf.Write(&fId,1);
		buf.Write(fCoord,3);
	}
	
	/** @brief Returns the identity of the current node*/
	int Id() const
	{
		return fId; 
	}
	
	void SetNodeId(int id) { fId = id;}
	
	/** @brief Initialize the data structure of the node. Creates a unique id for the node automatically*/
	void Initialize(TPZVec<REAL> &coord,TPZGeoMesh &mesh);
	/** @brief Initializes the data structure of the node. Assumes that the id provided by the user
     is unique for the mesh*/
	void Initialize(int id,TPZVec<REAL> &coord,TPZGeoMesh &mesh);
	
	/** @brief Initializes the node with data from a node from a different mesh */
	void Initialize(const TPZGeoNode &node, TPZGeoMesh &mesh);  
	/** @brief Returns i-th coordinate of the current node*/
	REAL Coord(int i) const;
	
	/** @brief Sets all coordinates into the current node. It gets the dim values from x */
	void SetCoord(const TPZVec<REAL> &x);
	
	
	/** @brief Set the i-th coordinate for current node*/
	void SetCoord(int i,REAL coord);
	
	/**
	 * @brief Fill the coordinates of the node
	 */
	void GetCoordinates(TPZVec<REAL> &co);
	/**
	 * @brief Print the node data into out.
	 */	
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

