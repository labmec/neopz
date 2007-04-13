// -*- c++ -*-
// $ Id: $

// TPZGeoTetrahedra.h: interface for the TPZGeoQuad class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TPZGEOTETRAHEDRAH
#define TPZGEOTETRAHEDRAH

#include "pznoderep.h"
#include "pzvec.h"
#include "pzeltype.h"


#include <string>


class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;
class TPZIntTetra3D;
class TPZGraphElT3d;
class TPZGeoMesh;

namespace pzgeom {

/// implements the geometry of a tetrahedral element
class TPZGeoTetrahedra  : public TPZNodeRep<4>
{
public:

	enum {NNodes = 4, NSides = 15};
  /**
  * Constructor with list of nodes
   */
 TPZGeoTetrahedra(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes>(nodeindexes)
 {
 }
  
  /**
  * Empty constructor
   */
 TPZGeoTetrahedra() : TPZNodeRep<NNodes>()
 {
 }
  
  /**
  * Constructor with node map
   */
 TPZGeoTetrahedra(const TPZGeoTetrahedra &cp,
                std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes>(cp,gl2lcNdMap)
 {
 }
  
  /**
  * Copy constructor
   */
 TPZGeoTetrahedra(const TPZGeoTetrahedra &cp) : TPZNodeRep<NNodes>(cp)
 {
 }


  /**
   * return the type of the element as specified in file pzeltype.h
   */
  static MElementType Type();// { return ETetraedro;}


/**
 * return the type of the element as specified in file pzeltype.h
 */
static MElementType Type(int side) ;
/*{
  switch(side) {
    case 0:
    case 1:
    case 2:
    case 3:
      return EPoint;
    case 4:
    case 5:
    case 6:
    case 7:
      return ETriangle;
    case 8:
      return ETetraedro;
    default:
      return ENoType;
  }
}*/

/**
 * returns the type name of the element
 */
static std::string TypeName() { return "Tetra";} 

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

	typedef TPZIntTetra3D IntruleType;
        typedef TPZGraphElT3d GraphElType;
};

};
#endif 
