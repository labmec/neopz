// -*- c++ -*-
// $ Id: $
//HEADER FILE FOR CLASS TPZGeoCube

#ifndef TPZGEOLINEARH
#define TPZGEOLINEARH


#include "pzvec.h"
#include "pzeltype.h"
#include "pzfmatrix.h"

#include <string>


class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;
class TPZInt1d;
class TPZGraphEl1dd;
namespace pzgeom {

/// implements the geometry of a one dimensional linear element
class TPZGeoLinear {

public:
	enum {NNodes = 2, NSides = 3};

  /**
   * return the type of the element as specified in file pzeltype.h
   */
static MElementType Type() { return EOned;}

/**
  * return the type of the element as specified in file pzeltype.h
  */
static MElementType Type(int side) {
  switch(side) {
    case 0:
    case 1:
      return EPoint;
    case 2:
      return EOned;
    default:
      return ENoType;
  }
}

/**
 * returns the type name of the element
 */
static std::string TypeName() { return "Linear";} 

static void X(TPZFMatrix &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);

static void Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi);

static void Jacobian(TPZFMatrix &nodes,TPZVec<REAL> &param,TPZFMatrix &jacobian,
				TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);

static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);

static TPZIntPoints *CreateSideIntegrationRule(int side, int order);

static int NSubElements();

 typedef TPZInt1d IntruleType;
 typedef TPZGraphEl1dd GraphElType;


};

inline void TPZGeoLinear::Jacobian(TPZFMatrix &coord,TPZVec<REAL> &param,TPZFMatrix &jacobian,
							TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) {


//   REAL mod1 = (coord(0,1)-coord(0,0))*0.5;
//   jacobian(0,0) = mod1;
//   detjac = mod1;
//   jacinv(0,0) = 1./detjac;
//   axes(0,0) = 1.;

  //VERSAO FUNCIONAL
  int ic;
  REAL v1[3];
  REAL mod1 = 0.;
  for(ic=0; ic<3; ic++) {
    v1[ic] = (coord(ic,1)-coord(ic,0))*0.5;
    mod1 += v1[ic]*v1[ic];
  }
  mod1 = sqrt(mod1);
  jacobian(0,0) = mod1;
  detjac = mod1;
  jacinv(0,0) = 1./mod1;

  for(ic=0; ic<3; ic++) {
    axes(0,ic) = v1[ic]/mod1;
  }

  // VERSAO ORIGINAL
//   TPZFNMatrix<9> phi(NNodes,1);
//   TPZFNMatrix<18> dphi(1,NNodes);
//   Shape(param,phi,dphi);

//   int ic;
//   TPZManVector<REAL,3> v1(3,0.);
//   REAL mod1 = 0.;

//   for(int i=0; i < NNodes; i++) {
//     for(ic = 0; ic < 3; ic++) {
//       v1[ic] += coord(ic,i)*dphi(0,i);
//     }
//   }

//   for(ic=0; ic<3; ic++) {
//     mod1 += v1[ic]*v1[ic];
//   }
//   mod1 = sqrt(mod1);
//   jacobian(0,0) = mod1;
//   detjac = mod1;
//   jacinv(0,0) = 1./detjac;


  //  axes.Zero();
//   for(ic=0; ic<3; ic++) {
//     axes(0,ic) = v1[ic]/mod1;
//   }
}

inline void TPZGeoLinear::X(TPZFMatrix &coord,TPZVec<REAL> &loc,TPZVec<REAL> &result){

  int ic;
  REAL xi = loc[0];
  for(ic=0; ic<3; ic++) result[ic] = coord(ic,0)*(1.-xi)*0.5+coord(ic,1)*(1.+xi)*0.5;

//   TPZFNMatrix<9> phi(NNodes,1);
//   TPZFNMatrix<18> dphi(1,NNodes);
//   Shape(loc,phi,dphi);
//   int in,ic;
//   for(in=0; in<3; in++) result[in] = 0.;
//   for(in = 0; in < NNodes; in++) {
//     for(ic=0; ic<3 ; ic++) {
//       result[ic] += coord(ic,in)*phi(in,0);
//     }
//   }

}

};
#endif

