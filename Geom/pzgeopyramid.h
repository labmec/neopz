// -*- c++ -*-
// $Id: pzgeopyramid.h,v 1.12 2009-11-23 19:30:53 phil Exp $

// TPZGeoPiramid.h: interface for the TPZGeoQuad class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TPZGEOTETRAPIRAMIDH
#define TPZGEOTETRAPIRAMIDH

#include "pznoderep.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "tpzpyramid.h"

class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {

/// implements the geometry of pyramid element
class TPZGeoPyramid  : public TPZNodeRep<5, pztopology::TPZPyramid>
{
public:

	enum {NNodes = 5};

  /**
  * Constructor with list of nodes
   */
 TPZGeoPyramid(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes, pztopology::TPZPyramid>(nodeindexes)
 {
 }
  
  /**
  * Empty constructor
   */
 TPZGeoPyramid() : TPZNodeRep<NNodes, pztopology::TPZPyramid>()
 {
 }
  
  /**
  * Constructor with node map
   */
 TPZGeoPyramid(const TPZGeoPyramid &cp,
                std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes, pztopology::TPZPyramid>(cp,gl2lcNdMap)
 {
 }
  
  /**
  * Copy constructor
   */
 TPZGeoPyramid(const TPZGeoPyramid &cp) : TPZNodeRep<NNodes, pztopology::TPZPyramid>(cp)
 {
 }

  /**
  * Copy constructor
   */
 TPZGeoPyramid(const TPZGeoPyramid &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZPyramid>(cp)
 {
 }

/**
 * returns the type name of the element
 */
static std::string TypeName() { return "Pyramid";} 

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
