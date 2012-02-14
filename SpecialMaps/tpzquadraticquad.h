/**
 * @file
 * @brief Contains the TPZQuadraticQuad class which defines a quadrilateral geometric element with quadratic map.
 */
#ifndef TPZQUADRATICQUAD_H
#define TPZQUADRATICQUAD_H

#include "pzgeoquad.h"
// #include "pzgmesh.h"
#include "pzgeoel.h"
// #include "tpzquadrilateral.h"
#include "pznoderep.h"

#include <iostream>

namespace pzgeom {
    
	/**
	 * @ingroup geometry
	 * @brief Defines a quadrilateral geometric element with quadratic map. \ref geometry "Geometry"
	 * @author Paulo Cesar de Alvarenga Lucci (Caju)
	 * @since 2007
	 */
	class TPZQuadraticQuad : public pzgeom::TPZNodeRep<8,pztopology::TPZQuadrilateral> {
		
	public:
		
		enum {NNodes = 8};
		
		bool IsLinearMapping() const {
			return false;
		}
		
		TPZQuadraticQuad(TPZVec<int> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(nodeindexes)
		{
		}
		
		TPZQuadraticQuad() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>()
		{
		}
		
		TPZQuadraticQuad(const TPZQuadraticQuad &cp,std::map<int,int> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp,gl2lcNdMap)
		{
		}
		
		TPZQuadraticQuad(const TPZQuadraticQuad &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
		{
		}
		
		TPZQuadraticQuad(const TPZQuadraticQuad &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
		{
		}
		
		/**
		 * @brief Returns the type name of the element
		 */
		static std::string TypeName() { return "Quad";} 
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, 
		 * a side and a boundary condition number
		 */
		static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
		
		static void Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi);
		
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
        
		static void X(TPZFMatrix &coord, TPZVec<REAL> &par, TPZVec<REAL> &result);
		
		static void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);
		
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
