//HEADER FILE FOR CLASS TPZGeoCube

#ifndef TPZGEOPOINTH
#define TPZGEOPOINTH

#include "pzvec.h"
#include "pzeltype.h"

class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;

class TPZGeoPoint {

public:
	enum {NNodes = 1, NSides = 1};

  /**
   * return the type of the element as specified in file pzeltype.h
   */
  static MElementType Type() { return EPoint;}

	static void X(TPZFMatrix &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);

	static void Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi);

	static void Jacobian(TPZFMatrix nodes,TPZVec<REAL> &param,TPZFMatrix &jacobian,
				TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);

	static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);

	static TPZIntPoints *CreateSideIntegrationRule(int side, int order);

	static int NSubElements();
};
#endif

