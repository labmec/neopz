// -*- c++ -*-
// $Id: TPZGeoCube.h,v 1.4 2005-02-28 22:08:04 phil Exp $

//HEADER FILE FOR CLASS TPZGeoCube

#ifndef TPZGEOCUBEH
#define TPZGEOCUBEH


#include "pzvec.h"
#include "pzeltype.h"

class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;
class TPZIntCube3D;
class TPZGraphElQ3dd;

namespace pzgeom {

/// implements the geometry of hexahedra element
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

  typedef TPZIntCube3D IntruleType;
  typedef TPZGraphElQ3dd GraphElType;
};

};
#endif

