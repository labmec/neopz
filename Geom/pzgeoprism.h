// -*- c++ -*-
// $Id: pzgeoprism.h,v 1.7 2005-04-25 02:07:50 phil Exp $

// TPZGeoPrism.h: interface for the TPZGeoQuad class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TPZGEOPRISMH
#define TPZGEOPRISMH

#include "pzvec.h"
#include "pzeltype.h"

class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;
class TPZIntPrism3D;

#include <string>


namespace pzgeom {

/// implements the geometry of a prism element
class TPZGeoPrism  
{
public:

	enum {NNodes = 6, NSides = 21};

  /**
   * return the type of the element as specified in file pzeltype.h
   */
  static MElementType Type() { return EPrisma;}

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
    case 5:
      return EPoint;
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
      return EOned;
    case 15:
      return ETriangle;
    case 16:
    case 17:
    case 18:
      return EQuadrilateral;
    case 19:
      return ETriangle;
    case 20:
      return EPrisma;
    default:
      return ENoType;
  }
}

/**
 * returns the type name of the element
 */
static std::string TypeName() { return "Prism";} 

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

  typedef TPZIntPrism3D IntruleType;

};

};
#endif 
