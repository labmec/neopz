// -*- c++ -*-
// $Id: pzgeoprism.h,v 1.14 2009-11-23 19:30:53 phil Exp $

// TPZGeoPrism.h: interface for the TPZGeoQuad class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TPZGEOPRISMH
#define TPZGEOPRISMH

#include "pznoderep.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "tpzprism.h"

class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

#include <string>


namespace pzgeom {

/// implements the geometry of a prism element
class TPZGeoPrism : public TPZNodeRep<6, pztopology::TPZPrism>  
{
public:

	enum {NNodes = 6};
  /**
  * Constructor with list of nodes
   */
 TPZGeoPrism(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes, pztopology::TPZPrism>(nodeindexes)
 {
 }
  
  /**
  * Empty constructor
   */
 TPZGeoPrism() : TPZNodeRep<NNodes, pztopology::TPZPrism>()
 {
 }
  
  /**
  * Constructor with node map
   */
 TPZGeoPrism(const TPZGeoPrism &cp,
                std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes, pztopology::TPZPrism>(cp,gl2lcNdMap)
 {
 }
  
  /**
  * Copy constructor
   */
 TPZGeoPrism(const TPZGeoPrism &cp) : TPZNodeRep<NNodes, pztopology::TPZPrism>(cp)
 {
 }

  /**
  * Copy constructor
   */
 TPZGeoPrism(const TPZGeoPrism &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZPrism>(cp)
 {
 }

/**
 * returns the type name of the element
 */
static std::string TypeName() { return "Prism";} 

  /** implementation of two-dimensional bilinear interpolation*/
static  void Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi);

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

protected:
   /**
    * This method apply an infinitesimal displacement in some points
    * to fix singularity problems when using MapToSide() method!
    * This points are CornerNodes, when projected in the opposing side
    */
    static void FixSingularity(int side, TPZVec<REAL>& OriginalPoint, TPZVec<REAL>& ChangedPoint);

	
public:
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
									  TPZVec<int>& nodeindexes,
									  int matid,
									  int& index);
	
	
};

};
#endif 
