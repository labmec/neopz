// -*- c++ -*-
// $Id: pzgeopyramid.h,v 1.5 2005-04-25 02:07:50 phil Exp $

// TPZGeoPiramid.h: interface for the TPZGeoQuad class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TPZGEOTETRAPIRAMIDH
#define TPZGEOTETRAPIRAMIDH

#include "pzvec.h"
#include "pzeltype.h"

class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;
class TPZIntPyram3D;

namespace pzgeom {

/// implements the geometry of pyramid element
class TPZGeoPyramid  
{
public:

	enum {NNodes = 5, NSides = 19};

  /**
   * return the type of the element as specified in file pzeltype.h
   */
static MElementType Type() { return EPiramide;}

/**
  * return the type of the element as specified in file pzeltype.h
  */
static MElementType Type(int side) {
  switch(side) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
      return EPoint;
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
      return EOned;
    case 13:
      return EQuadrilateral;
    case 14:
    case 15:
    case 16:
    case 17:
      return ETriangle;
    case 18:
      return EPiramide;
    default:
      return ENoType;
  }
}

/**
 * returns the type name of the element
 */
static std::string TypeName() { return "Pyramid";} 

	/** implementation of two-dimensional bilinear interpolation*/
	static  void Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi);

	/**Computes the jacobian*/
	static  void Jacobian(TPZFMatrix & coord, TPZVec<REAL>& par, TPZFMatrix &jacobian, TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);

	/**Computes the geometric location*/
	static  void X(TPZFMatrix & coord, TPZVec<REAL>& par, TPZVec<REAL> &result);

	/**
	* Method which creates a geometric boundary condition 
	* element based on the current geometric element, 
	* a side and a boundary condition number
	*/
	static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);

	/**
	* Create an integration rule 
	* @param order order of the integration rule to be created
	* @param side side to create integration rule
	*/
	static TPZIntPoints * CreateSideIntegrationRule(int side, int order);

	typedef TPZIntPyram3D IntruleType;
};

};

#endif 
