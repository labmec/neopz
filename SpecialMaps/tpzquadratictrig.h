#ifndef TPZQUADRATICTRIG_H
#define TPZQUADRATICTRIG_H

#include "pznoderep.h"
#include "tpztriangle.h"


class TPZGeoEl;

/**
 / LabMeC - FEC - UNICAMP
 / 2007
 */

namespace pzgeom
{
	
	/**
	 * @author Paulo Cesar de Alvarenga Lucci (Caju)
	 * @ingroup geometry
	 * @brief Defines a triangular geometric element with quadratic map
	 */
	class TPZQuadraticTrig : public pzgeom::TPZNodeRep<6,pztopology::TPZTriangle> {
		
	public:
		
		enum {NNodes = 6};
		
		bool IsLinearMapping() const {
			return false;
		}
		
		TPZQuadraticTrig(TPZVec<int> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTriangle>(nodeindexes) {
		}
		
		TPZQuadraticTrig() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTriangle>() {
		}
		
		TPZQuadraticTrig(const TPZQuadraticTrig &cp,std::map<int,int> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp,gl2lcNdMap) {
		}
		
		TPZQuadraticTrig(const TPZQuadraticTrig &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp) {
		}
		
		TPZQuadraticTrig(const TPZQuadraticTrig &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp) {
		}
		/**
		 * @brief Returns the type name of the element
		 */
		static std::string TypeName() { return "Triangle";}
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, 
		 * a side and a boundary condition number
		 */
		static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
		
		static void Shape(TPZVec<REAL> &param,TPZFMatrix &phi,TPZFMatrix &dphi);
		
		static void X(TPZFMatrix &coord, TPZVec<REAL> &par, TPZVec< REAL > &result);
		
		static void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);
		
	public:
		/**
		 * @brief Creates a geometric element according to the type of the father element
		 */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<int>& nodeindexes,
										  int matid,
										  int& index);
		
		
	};
    
};

#endif
