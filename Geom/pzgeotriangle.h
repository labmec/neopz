// -*- c++ -*-
// $ Id: $
// TPZGeoQuad.h: interface for the TPZGeoQuad class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TPZGEOTRIANGLEH
#define TPZGEOTRIANGLEH

#include "pzvec.h"
#include "pzeltype.h"

#include <string>

class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;
class TPZGraphElTd;
class TPZIntTriang;

namespace pzgeom {

/// implements the geometry of a triangle element
class TPZGeoTriangle  
{
public:

	enum {NNodes = 3, NSides = 7};

  /**
   * return the type of the element as specified in file pzeltype.h
   */
  static MElementType Type();// { return ETriangle;}

/**
  * return the type of the element as specified in file pzeltype.h
  */
static MElementType Type(int side) ;/*{
  switch(side) {
    case 0:
    case 1:
    case 2:
      return EPoint;
    case 3:
    case 4:
    case 5:
      return EOned;
    case 6:
      return ETriangle;
    default:
      return ENoType;
  }
}*/

/**
 * returns the type name of the element
 */
static std::string TypeName() { return "Triangle";} 

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

  typedef TPZIntTriang IntruleType;
  typedef TPZGraphElTd GraphElType;


};

};
#endif 
