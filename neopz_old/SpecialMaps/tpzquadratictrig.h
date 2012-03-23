#ifndef TPZQUADRATICTRIG_H
#define TPZQUADRATICTRIG_H

#include "pznoderep.h"
#include "tpztriangle.h"


class TPZGeoEl;

using namespace std;
using namespace pzgeom;
using namespace pztopology;

  /**
  / Class made by Paulo Cesar de Alvarenga Lucci (Caju)
  / LabMeC - FEC - UNICAMP
  / 2007
 */

class TPZQuadraticTrig : public TPZNodeRep<6,TPZTriangle> {

     public:

     enum {NNodes = 6};

     bool IsLinearMapping() const {
          return false;
     }

     TPZQuadraticTrig(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes,pztopology::TPZTriangle>(nodeindexes) {
     }

     TPZQuadraticTrig() : TPZNodeRep<NNodes,pztopology::TPZTriangle>() {
     }

     TPZQuadraticTrig(const TPZQuadraticTrig &cp,std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp,gl2lcNdMap) {
     }

     TPZQuadraticTrig(const TPZQuadraticTrig &cp) : TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp) {
     }

     TPZQuadraticTrig(const TPZQuadraticTrig &cp, TPZGeoMesh &) : TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp) {
     }
     /**
     * returns the type name of the element
     */
     static std::string TypeName() { return "Triangle";}

     /**
     * Method which creates a geometric boundary condition 
     * element based on the current geometric element, 
     * a side and a boundary condition number
     */
     static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);

     static void Shape(TPZVec<REAL> &param,TPZFMatrix &phi,TPZFMatrix &dphi);

     static void X(TPZFMatrix &coord, TPZVec<REAL> &par, TPZVec< REAL > &result);

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
