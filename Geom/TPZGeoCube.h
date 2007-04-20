// -*- c++ -*-
// $Id: TPZGeoCube.h,v 1.8 2007-04-20 18:31:10 caju Exp $

//HEADER FILE FOR CLASS TPZGeoCube

#ifndef TPZGEOCUBEH
#define TPZGEOCUBEH


#include "pzvec.h"
#include "pzeltype.h"
#include "pznoderep.h"
#include "tpzcube.h"

#include <string>


class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {

/// implements the geometry of hexahedra element
  class TPZGeoCube : public TPZNodeRep<8, pztopology::TPZCube> {

public:
	enum {NNodes = 8};

  /**
  * Constructor with list of nodes
   */
 TPZGeoCube(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes, pztopology::TPZCube>(nodeindexes)
 {
 }
  
  /**
  * Empty constructor
   */
 TPZGeoCube() : TPZNodeRep<NNodes, pztopology::TPZCube>()
 {
 }
  
  /**
  * Constructor with node map
   */
 TPZGeoCube(const TPZGeoCube &cp,
                std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes, pztopology::TPZCube>(cp,gl2lcNdMap)
 {
 }
  
  /**
  * Copy constructor
   */
 TPZGeoCube(const TPZGeoCube &cp) : TPZNodeRep<NNodes, pztopology::TPZCube>(cp)
 {
 }



/**
 * returns the type name of the element
 */
static std::string TypeName() { return "Hexa";} 

static void X(TPZFMatrix &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);

static void Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi);

static void Jacobian(TPZFMatrix &nodes,TPZVec<REAL> &param,TPZFMatrix &jacobian,
				TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);

static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);

};

};
#endif

