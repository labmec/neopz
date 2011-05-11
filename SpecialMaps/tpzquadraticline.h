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

  /**
  / Class made by Paulo Cesar de Alvarenga Lucci (Caju)
  / LabMeC - FEC - UNICAMP
  / 2007
 */

class TPZQuadraticLine : public pzgeom::TPZNodeRep<3,pztopology::TPZLine> {

public:

enum {NNodes = 3};

  bool IsLinearMapping() const {
      return false;
  }

  TPZQuadraticLine(TPZVec<int> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes)
  {
  }

  TPZQuadraticLine() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>()
  {
  }

  TPZQuadraticLine(const TPZQuadraticLine &cp,std::map<int,int> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap)
  {
  }

  TPZQuadraticLine(const TPZQuadraticLine &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp)
  {
  }

  TPZQuadraticLine(const TPZQuadraticLine &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp)
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
