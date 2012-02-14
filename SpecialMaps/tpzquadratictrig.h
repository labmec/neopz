/**
 * @file
 * @brief Contains the TPZQuadraticTrig class which defines a triangular geometric element with quadratic map.
 */
#ifndef TPZQUADRATICTRIG_H
#define TPZQUADRATICTRIG_H

#include "pznoderep.h"
#include "tpztriangle.h"

//class TPZGeoEl;

	
	/**
	 * @author Paulo Cesar de Alvarenga Lucci (Caju)
	 * @since 2007
	 * @ingroup geometry
	 * @brief Defines a triangular geometric element with quadratic map. \ref geometry "Geometry"
	 */
namespace pzgeom {
    
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
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Triangle";}
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, 
		 * a side and a boundary condition number
		 */
		static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
		
		static void Shape(TPZVec<REAL> &param,TPZFMatrix &phi,TPZFMatrix &dphi);
		
		/* brief compute the coordinate of a point given in parameter space */
        void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
		
        /* @brief compute the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
        }
        
		static void X(TPZFMatrix &coord, TPZVec<REAL> &par, TPZVec< REAL > &result);
		
		static void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);
		
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<int>& nodeindexes,
										  int matid, int& index);

	};

};

#endif
