//
// C++ Interface: tpzcube
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PZTOPOLOGYTPZCUBE_H
#define PZTOPOLOGYTPZCUBE_H


#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzeltype.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

class TPZIntPoints;
class TPZIntCube3D;
class TPZGraphElQ3dd;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

namespace pztopology {

/**
@author Philippe R. B. Devloo
*/
class TPZCube{
public:

    enum {NSides = 27, NCornerNodes = 8, Dimension = 3};

    TPZCube();

    virtual ~TPZCube();


static void LowerDimensionSides(int side,TPZStack<int> &smallsides);
static void LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget);

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
static REAL RefElVolume(){return 8.0;}

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
 * @param side side to which the point will be tranformed (0<=side<=26)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
static TPZTransform TransformElementToSide(int side);
/**
 * Returns the transformation which transform a point from the side to the interior of the element
 * @param side side from which the point will be tranformed (0<=side<=26)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
static TPZTransform TransformSideToElement(int side);

/** Verifies if the parametric point pt is in the element parametric domain
*/
static bool IsInParametricDomain(TPZVec<REAL> &pt, REAL tol = 1e-6);

static TPZIntPoints *CreateSideIntegrationRule(int side, int order);

  typedef TPZIntCube3D IntruleType;
  typedef TPZGraphElQ3dd GraphElType;

  /**
   * return the type of the element as specified in file pzeltype.h
   */
static MElementType Type();// { return ECube;}

/**
 * return the type of the element as specified in file pzeltype.h
 */
static MElementType Type(int side);

/**
 * Number of connects of the element (27)
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

/// function pointer which determines the type of computational element
/**
 * function pointer which determines what type of computational element will be created
 */
static TPZCompEl *(*fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index);


protected:
/** 
 * Data structure which defines the hexahedral transformations
 */
static int FaceNodes[6][4];

/** 
 * Data structure which defines the hexahedral transformations
 */
static int SideNodes[12][2];

/** 
 * Data structure which defines the hexahedral transformations
 */
static int ShapeFaceId[6][2];

};

}

#endif
