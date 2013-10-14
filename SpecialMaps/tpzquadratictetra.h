/**
 * @file
 * @brief Contains the TPZQuadraticTetra class which defines a tetrahedral geometric element with quadratic map.
 */
 
#ifndef TPZQUADRATICTETRA_H
#define TPZQUADRATICTETRA_H

#include "pzgeotetrahedra.h"
#include "pzgeoel.h"
#include "pznoderep.h"

#include <iostream>

namespace pzgeom {
    
	/**
	 * @author Paulo Cesar de Alvarenga Lucci (Caju)
	 * @since 2007
	 * @ingroup geometry
	 * @brief Defines a tetrahedral geometric element with quadratic map. \ref geometry "Geometry"
	 */
	class TPZQuadraticTetra : public pzgeom::TPZNodeRep<10,pztopology::TPZTetrahedron> {
		
	public:
		/** @brief Number of nodes */
		enum {NNodes = 10};
        
        //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
        
		/** @brief It is quadratic mapping */
		bool IsLinearMapping() const {
			return false;
		}
		/** @brief Constructor for node indexes given */
		TPZQuadraticTetra(TPZVec<long> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(nodeindexes) {
		}
		/** @brief Default constructor */
		TPZQuadraticTetra() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>() {
		}
		/** @brief Constructor for node map given */
		TPZQuadraticTetra(const TPZQuadraticTetra &cp,std::map<long,long> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp,gl2lcNdMap) {
		}
		/** @brief Copy constructor */
		TPZQuadraticTetra(const TPZQuadraticTetra &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp) {
		}
		/** @brief Copy constructor */		
		TPZQuadraticTetra(const TPZQuadraticTetra &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp) {
		}
		/** @brief Destructor */
		virtual	~TPZQuadraticTetra();
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() {
			return "Tetra";
		}
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, 
		 * a side and a boundary condition number
		 */
		static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
		
		static void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
		
		/** @brief compute the coordinate of a point given in parameter space */
        void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
		
        /** @brief compute the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
        }
        
		static void X(TPZFMatrix<REAL> &coord, TPZVec<REAL> &loc,TPZVec<REAL> &result);
		
		static void Jacobian(TPZFMatrix<REAL> &coord, TPZVec<REAL> &param, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv);
		
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<long>& nodeindexes,
										  int matid,
										  long& index);	
};

};

#endif
