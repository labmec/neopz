// -*- c++ -*-

// $ Id: $

// TPZGeoQuad.h: interface for the TPZGeoQuad class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TPZGEOQUADH
#define TPZGEOQUADH

#include "pznoderep.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "tpzquadrilateral.h"

#include <string>


class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {

/// implements the geometry of a quadrilateral element
class TPZGeoQuad  : public TPZNodeRep<4, pztopology::TPZQuadrilateral>
{
public:

	enum {NNodes = 4};
	enum {NVectors = 18};
  /**
  * Constructor with list of nodes
   */
 TPZGeoQuad(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(nodeindexes)
 {
 }
  
  /**
  * Empty constructor
   */
 TPZGeoQuad() : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>()
 {
 }
  
  /**
  * Constructor with node map
   */
 TPZGeoQuad(const TPZGeoQuad &cp,
                std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp,gl2lcNdMap)
 {
 }
  
  /**
  * Copy constructor
   */
 TPZGeoQuad(const TPZGeoQuad &cp) : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
 {
 }

  /**
  * Copy constructor
   */
 TPZGeoQuad(const TPZGeoQuad &cp, TPZGeoMesh &) : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
 {
 }

/**
 * returns the type name of the element
 */
static std::string TypeName() { return "Quad";} 

  /** implementation of two-dimensional bilinear interpolation*/
static  void Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi);
	
	/** implementation of normal vector to Hdiv space*/
	/** construct the normal vector for element hdiv*/
static void VecHdiv(TPZFMatrix & coord,TPZFMatrix &fNormalVec,TPZVec<int> & fVectorSide);
	/** computes the vecorial product of the two vectors*/ 
static void VectorialProduct(TPZVec<REAL> &v1, TPZVec<REAL> &v2,TPZVec<REAL> &result);
	/** computes */
static void ComputeNormal(TPZVec<REAL> &p1, TPZVec<REAL> &p2,TPZVec<REAL> &p3,TPZVec<REAL> &result);
	
	

  /**Computes the jacobian*/
static  void Jacobian(TPZFMatrix & coord, TPZVec<REAL>& par, TPZFMatrix &jacobian, TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);

  /**Computes the geometric location*/
static  void X(TPZFMatrix & coord, TPZVec<REAL>& par, TPZVec<REAL> &result);

/**
 * returns the projection of a given point from "NSide - 1" side to "side".
 */
static bool MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix &JacToSide);

  /**
   * Method which creates a geometric boundary condition 
   * element based on the current geometric element, 
   * a side and a boundary condition number
   */
static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);


};

};
#endif 
