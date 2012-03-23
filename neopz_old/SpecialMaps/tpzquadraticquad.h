#ifndef TPZQUADRATICQUAD_H
#define TPZQUADRATICQUAD_H

// #include "pzfmatrix.h"
// #include "pzvec.h"
#include "pzgeoquad.h"
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

class TPZQuadraticQuad : public TPZNodeRep<8,TPZQuadrilateral> {

public:

    enum {NNodes = 8};

     bool IsLinearMapping() const {
          return false;
     }

     TPZQuadraticQuad(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(nodeindexes)
     {
     }

     TPZQuadraticQuad() : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>()
     {
     }

     TPZQuadraticQuad(const TPZQuadraticQuad &cp,std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp,gl2lcNdMap)
     {
     }

     TPZQuadraticQuad(const TPZQuadraticQuad &cp) : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
     {
     }

     TPZQuadraticQuad(const TPZQuadraticQuad &cp, TPZGeoMesh &) : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
     {
     }

     /**
      * returns the type name of the element
      */
     static std::string TypeName() { return "Quad";} 

     /**
      * Method which creates a geometric boundary condition 
      * element based on the current geometric element, 
      * a side and a boundary condition number
      */
     static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);

     static void Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi);

     static void X(TPZFMatrix &coord, TPZVec<REAL> &par, TPZVec<REAL> &result);

     static void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);

	
public:
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
									  TPZVec<int>& nodeindexes,
									  int matid,
									  int& index);
	
};

#endif
