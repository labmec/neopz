// -*- c++ -*-
// $ Id: $
//HEADER FILE FOR CLASS TPZGeoCube

#ifndef TPZGEOPOINTH
#define TPZGEOPOINTH

#include "pzvec.h"
#include "pzeltype.h"

class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;
class TPZInt1Point;
class TPZGraphEl1dd;

#include <string>


/// groups all classes which model the geometry
/**
* Objects of this class implement the mapping between the master element
* and deformed element
* These classes are used as template arguments of @seealso TPZGeoElement and
* @seealso TPZIntelGen
*/
namespace pzgeom {

/// implements the geometry of a point element
class TPZGeoPoint {

public:
	enum {NNodes = 1, NSides = 1};

  /**
   * return the type of the element as specified in file pzeltype.h
   */
static MElementType Type() ;//{ return EPoint;}


/**
  * return the type of the element as specified in file pzeltype.h
  */
static MElementType Type(int side) ;/*{
  switch(side) {
    case 0:
      return EPoint;
    default:
      return ENoType;
  }
}
*/

/**
 * returns the type name of the element
 */
static std::string TypeName() { return "Point";} 

	static void X(TPZFMatrix &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);

	static void Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi);

	static void Jacobian(TPZFMatrix nodes,TPZVec<REAL> &param,TPZFMatrix &jacobian,
				TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);

	static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);

	static TPZIntPoints *CreateSideIntegrationRule(int side, int order);

	static int NSubElements();

	typedef TPZInt1Point IntruleType;
	//	typedef TPZGraphEl1dd GraphElType;
};

};
#endif

