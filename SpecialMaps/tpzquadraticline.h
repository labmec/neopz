#ifndef TPZQUADRATICLINE_H
#define TPZQUADRATICLINE_H

// #include "pzfmatrix.h"
// #include "pzvec.h"
#include "TPZGeoLinear.h"
// #include "pzgmesh.h"
#include "pzgeoel.h"
// #include "tpzquadrilateral.h"
#include "pznoderep.h"

#include <iostream>

using namespace std;
using namespace pzgeom;
using namespace pztopology;

  /**
  / Class made by Paulo Cesar de Alvarenga Lucci (Caju)
  / LabMeC - FEC - UNICAMP
  / 2007
 */

class TPZQuadraticLine : public TPZNodeRep<3,TPZLine> {

public:

enum {NNodes = 3};

  bool IsLinearMapping() const {
      return false;
  }

  TPZQuadraticLine(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes)
  {
  }

  TPZQuadraticLine() : TPZNodeRep<NNodes,pztopology::TPZLine>()
  {
  }

  TPZQuadraticLine(const TPZQuadraticLine &cp,std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap)
  {
  }

  TPZQuadraticLine(const TPZQuadraticLine &cp) : TPZNodeRep<NNodes,pztopology::TPZLine>(cp)
  {
  }

  TPZQuadraticLine(const TPZQuadraticLine &cp, TPZGeoMesh &) : TPZNodeRep<NNodes,pztopology::TPZLine>(cp)
  {
  }

  /**
  * returns the type name of the element
  */
  static std::string TypeName() { return "Line";} 

  /**
  * Method which creates a geometric boundary condition 
  * element based on the current geometric element, 
  * a side and a boundary condition number
  */
  static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
	
public:
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
									  TPZVec<int>& nodeindexes,
									  int matid,
									  int& index);
	
	

  static void Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi);

  static void X(TPZFMatrix &coord, TPZVec<REAL> &par, TPZVec<REAL> &result);

  static void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);
};

#endif
