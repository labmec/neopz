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
#include "tpztetrahedron.h"

#include <string>


class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {

/// implements the geometry of a tetrahedral element
class TPZGeoTetrahedra  : public TPZNodeRep<4,pztopology::TPZTetrahedron>
{
public:

	enum {NNodes = 4};
  /**
  * Constructor with list of nodes
   */
 TPZGeoTetrahedra(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(nodeindexes)
 {
 }
  
  /**
  * Empty constructor
   */
 TPZGeoTetrahedra() : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>()
 {
 }
  
  /**
  * Constructor with node map
   */
 TPZGeoTetrahedra(const TPZGeoTetrahedra &cp,
                std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp,gl2lcNdMap)
 {
 }

  /**
  * Copy constructor
   */
 TPZGeoTetrahedra(const TPZGeoTetrahedra &cp) : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp)
 {
 }

  /**
  * Copy constructor
   */
 TPZGeoTetrahedra(const TPZGeoTetrahedra &cp, TPZGeoMesh &) : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp)
 {
 }

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
 * returns the projection of a given point from "NSide - 1" side to "side".
 */
static void MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix &JacToSide);

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

};

};
#endif 
