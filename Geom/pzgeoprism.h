// -*- c++ -*-
// $Id: pzgeoprism.h,v 1.6 2005-02-28 22:08:03 phil Exp $

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
