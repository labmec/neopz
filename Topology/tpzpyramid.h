//
// C++ Interface: tpzpyramid
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PZTOPOLOGYTPZPYRAMID_H
#define PZTOPOLOGYTPZPYRAMID_H

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzquad.h"
#include "pzeltype.h"

class TPZIntPoints;
class TPZIntPyram3D;

namespace pztopology {

/**
@author Philippe R. B. Devloo
*/
class TPZPyramid{
public:

    enum {NSides = 19, NCornerNodes = 5, Dimension = 3};

    TPZPyramid();

    ~TPZPyramid();

 /**
  * returns all sides whose closure contains side
  * @param side smaller dimension side
  * @param high vector which will contain all sides whose closure contain sidefrom
  */
static void HigherDimensionSides(int side, TPZStack<int> &high);
/**
  * return the number of nodes (not connectivities) associated with a side
  */
static int NSideNodes(int side);
 /**
  * returns the local node number of the node "node" along side "side"
  */
static int SideNodeLocId(int side, int node);
/**
 * returns the barycentric coordinates in the master element space of the original element
 */
 static void CenterPoint(int side, TPZVec<REAL> &center);
/**volume of the master element*/
static REAL RefElVolume(){return (4./3.);}
 /**
  * returns the dimension of the side
  */
static int SideDimension(int side);
/**
  * returns the transformation which takes a point from the side sidefrom ot
  * the side sideto
  * @param sidefrom side where the point resides
  * @param sideto side whose closure contains sidefrom
  */
static TPZTransform SideToSideTransform(int sidefrom, int sideto);
/**
 * Returns the transformation which projects a point from the interior of the element to the side
 * @param side side to which the point will be tranformed (0<=side<=20)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
static TPZTransform TransformElementToSide(int side);
/**
 * Returns the transformation which transform a point from the side to the interior of the element
 * @param side side from which the point will be tranformed (0<=side<=20)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
static TPZTransform TransformSideToElement(int side);


/**
  * Create an integration rule 
  * @param order order of the integration rule to be created
  * @param side side to create integration rule
  */
static TPZIntPoints * CreateSideIntegrationRule(int side, int order);

typedef TPZIntPyram3D IntruleType;


  /**
   * return the type of the element as specified in file pzeltype.h
   */
static MElementType Type();// { return EPiramide;}

/**
  * return the type of the element as specified in file pzeltype.h
  */
static MElementType Type(int side);

/**
 * Number of connects of the element (21)
 * @return number of connects of the element
 */
static int NConnects();
 /**
  * return the number of nodes (not connectivities) associated with a side
  */
static int NSideConnects(int side);
 /**
  * returns the local connect number of the connect "c" along side "side"
  */
static int SideConnectLocId(int side, int c);


protected:
/** 
 * Data structure which defines the pyramid transformations and topology
 */
static int FaceNodes[5][4];
/** 
 * Data structure which defines the pyramid transformations and topology
 */
static int SideNodes[8][2];
/** 
 * Data structure which defines the pyramid transformations and topology
 */
static int ShapeFaceId[5][4];

};

}

#endif
