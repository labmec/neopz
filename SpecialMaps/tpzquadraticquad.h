/**
 * @file
 * @brief Contains the TPZQuadraticQuad class which defines a quadrilateral geometric element with quadratic map.
 */
 
#ifndef TPZQUADRATICQUAD_H
#define TPZQUADRATICQUAD_H

#include "pzgeoquad.h"
#include "pzgeoel.h"
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
		/** @brief Number of nodes */
		enum {NNodes = 8};
        
        //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
        
		/** @brief Constructor from node indexes */
		TPZQuadraticQuad(TPZVec<long> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(nodeindexes)
		{
		}
		/** @brief Default constructor */
		TPZQuadraticQuad() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>()
		{
		}
		/** @brief Constructor over node map */
		TPZQuadraticQuad(const TPZQuadraticQuad &cp,std::map<long,long> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp,gl2lcNdMap)
		{
		}
		/** @brief Copy constructor */
		TPZQuadraticQuad(const TPZQuadraticQuad &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
		{
		}
		/** @brief Copy constructor */
		TPZQuadraticQuad(const TPZQuadraticQuad &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
		{
		}
		
		/**
		 * @brief Returns the type name of the element
		 */
		static std::string TypeName() { return "Quad";}
        
        static bool IsLinearMapping(int side)
        {
            return false;
        }
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, 
		 * a side and a boundary condition number
		 */
		static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
		
        template<class T>
		static void Shape(TPZVec<T> &x,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
		
		/* brief compute the coordinate of a point given in parameter space */
        template<class T>
        void X(const TPZGeoEl &gel,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
        
        template<class T>
        void GradX(const TPZGeoEl &gel, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
        {
            DebugStop();
        }
		
        /* @brief compute the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
        }
        
        template<class T>
		static void X(TPZFMatrix<REAL> &coord, TPZVec<T> &par, TPZVec<T> &result);
		
		static void Jacobian(TPZFMatrix<REAL> &coord, TPZVec<REAL> &par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv);
		
		/**
		 * @brief Creates a geometric element according to the type of the father element
		 */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<long>& nodeindexes,
										  int matid,
										  long& index);
        
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

};

};

#endif
