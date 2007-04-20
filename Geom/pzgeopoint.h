// -*- c++ -*-
// $ Id: $
//HEADER FILE FOR CLASS TPZGeoCube

#ifndef TPZGEOPOINTH
#define TPZGEOPOINTH

#include "pznoderep.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "pzintel.h"
#include "tpzpoint.h"

class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;
class TPZInt1Point;
class TPZGraphEl1dd;
class TPZGeoMesh;

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
  class TPZGeoPoint : public TPZNodeRep<1, pztopology::TPZPoint> {

public:
	enum {NNodes = 1};

  /**
  * Constructor with list of nodes
   */
 TPZGeoPoint(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes, pztopology::TPZPoint>(nodeindexes)
 {
 }
  
  /**
  * Empty constructor
   */
 TPZGeoPoint() : TPZNodeRep<NNodes, pztopology::TPZPoint>()
 {
 }
  
  /**
  * Constructor with node map
   */
 TPZGeoPoint(const TPZGeoPoint &cp,
                std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes, pztopology::TPZPoint>(cp,gl2lcNdMap)
 {
 }
  
  /**
  * Copy constructor
   */
 TPZGeoPoint(const TPZGeoPoint &cp) : TPZNodeRep<NNodes, pztopology::TPZPoint>(cp)
 {
 }


/**
 * returns the type name of the element
 */
        static std::string TypeName() { return "Point";} 

	static void X(TPZFMatrix &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);

	static void Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi);

	static void Jacobian(TPZFMatrix nodes,TPZVec<REAL> &param,TPZFMatrix &jacobian,
				TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);

	static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);

//	static int NSubElements();

	//	typedef TPZGraphEl1dd GraphElType;
};

};
#endif

