//HEADER FILE FOR CLASS TPZGeoCube

#ifndef TPZGEOCUBEH
#define TPZGEOCUBEH


#include "pzvec.h"
#include "pzeltype.h"

class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;

class TPZGeoCube {

public:
	enum {NNodes = 8, NSides = 27};

  /**
   * return the type of the element as specified in file pzeltype.h
   */
  static MElementType Type() { return ECube;}

static void X(TPZFMatrix &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);

static void Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi);

static void Jacobian(TPZFMatrix &nodes,TPZVec<REAL> &param,TPZFMatrix &jacobian,
				TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);

static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);

static TPZIntPoints *CreateSideIntegrationRule(int side, int order);

};
#endif

