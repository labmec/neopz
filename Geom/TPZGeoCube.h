// -*- c++ -*-
// $Id: TPZGeoCube.h,v 1.7 2007-04-13 18:59:40 phil Exp $

//HEADER FILE FOR CLASS TPZGeoCube

#ifndef TPZGEOCUBEH
#define TPZGEOCUBEH


#include "pzvec.h"
#include "pzeltype.h"
#include "pznoderep.h"

#include <string>


class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;
class TPZIntCube3D;
class TPZGraphElQ3dd;
class TPZGeoMesh;

namespace pzgeom {

/// implements the geometry of hexahedra element
  class TPZGeoCube : public TPZNodeRep<8> {

public:
	enum {NNodes = 8, NSides = 27};

  /**
  * Constructor with list of nodes
   */
 TPZGeoCube(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes>(nodeindexes)
 {
 }
  
  /**
  * Empty constructor
   */
 TPZGeoCube() : TPZNodeRep<NNodes>()
 {
 }
  
  /**
  * Constructor with node map
   */
 TPZGeoCube(const TPZGeoCube &cp,
                std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes>(cp,gl2lcNdMap)
 {
 }
  
  /**
  * Copy constructor
   */
 TPZGeoCube(const TPZGeoCube &cp) : TPZNodeRep<NNodes>(cp)
 {
 }


  /**
   * return the type of the element as specified in file pzeltype.h
   */
static MElementType Type();// { return ECube;}

/**
 * return the type of the element as specified in file pzeltype.h
 */
static MElementType Type(int side);/* {
  switch(side) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
      return EPoint;
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
    case 16:
    case 17:
    case 18:
    case 19:        
      return EOned;
    case 20:
    case 21:
    case 22:
    case 23:
    case 24:
    case 25:
      return EQuadrilateral;
    case 26:
      return ECube;
    default:
      return ENoType;
  }
}
*/

/**
 * returns the type name of the element
 */
static std::string TypeName() { return "Hexa";} 

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

