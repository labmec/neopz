//
// C++ Interface: tpztetrahedron
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PZTOPOLOGYTPZTETRAHEDRON_H
#define PZTOPOLOGYTPZTETRAHEDRON_H


#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzeltype.h"

class TPZIntPoints;
class TPZIntTetra3D;
class TPZGraphElT3d;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

namespace pztopology {

/**
@author Philippe R. B. Devloo
*/
class TPZTetrahedron{
public:

    enum {NSides = 15, NCornerNodes = 4, Dimension = 3};

    TPZTetrahedron();

    virtual ~TPZTetrahedron();

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
    static REAL RefElVolume(){return (1./6.);}
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
    
    /**
    * Create an integration rule 
    * @param order order of the integration rule to be created
    * @param side side to create integration rule
    */
    static TPZIntPoints * CreateSideIntegrationRule(int side, int order);

    typedef TPZIntTetra3D IntruleType;
    typedef TPZGraphElT3d GraphElType;


  /**
   * return the type of the element as specified in file pzeltype.h
   */
  static MElementType Type();// { return ETetraedro;}


/**
 * return the type of the element as specified in file pzeltype.h
 */
static MElementType Type(int side) ;


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
    * Data structure which defines the tetrahedral transformations and topology
    */
    static int FaceNodes[4][3];
    /** 
    * Data structure which defines the tetrahedral transformations and topology
    */
    static int SideNodes[6][2];
    /** 
    * Data structure which defines the tetrahedral transformations and topology
    */
    static int ShapeFaceId[4][3];

};

}

#endif
