// -*- c++ -*-
// $Id: pzgeopyramid.h,v 1.3 2003-10-20 02:13:25 phil Exp $

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

class TPZGeoPyramid  
{
public:

	enum {NNodes = 5, NSides = 19};

  /**
   * return the type of the element as specified in file pzeltype.h
   */
  static MElementType Type() { return EPiramide;}

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

#endif 
