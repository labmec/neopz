#ifndef TPZQUADRATICTETRA_H
#define TPZQUADRATICTETRA_H

// #include "pzfmatrix.h"
// #include "pzvec.h"
#include "pzgeotetrahedra.h"
// #include "tpztetrahedron.h"
// #include "pzgmesh.h"
#include "pzgeoel.h"
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

class TPZQuadraticTetra : public TPZNodeRep<10,TPZTetrahedron> {

public:

    enum {NNodes = 10};

     bool IsLinearMapping() const {
          return false;
     }

     TPZQuadraticTetra(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(nodeindexes) {
     }

     TPZQuadraticTetra() : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>() {
     }

     TPZQuadraticTetra(const TPZQuadraticTetra &cp,std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp,gl2lcNdMap) {
     }

     TPZQuadraticTetra(const TPZQuadraticTetra &cp) : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp) {
     }

     TPZQuadraticTetra(const TPZQuadraticTetra &cp, TPZGeoMesh &) : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp) {
     }

virtual	~TPZQuadraticTetra();
     /**
     * returns the type name of the element
     */
     static std::string TypeName() {
          return "Tetra";
     }

     /**
          * Method which creates a geometric boundary condition 
          * element based on the current geometric element, 
          * a side and a boundary condition number
          */
     static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);

     static void Shape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi);

     static void X(TPZFMatrix &coord, TPZVec<REAL> &loc,TPZVec<REAL> &result);

     static void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &param, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);
	
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
