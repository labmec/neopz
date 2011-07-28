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

/**
 * @author Paulo Cesar de Alvarenga Lucci (Caju)
 * LabMeC - FEC - UNICAMP
 * 2007
 * @ingroup geometry
 * @brief Defines a tetrahedral geometric element with quadratic map
 */

namespace pzgeom
{
	
	class TPZQuadraticTetra : public pzgeom::TPZNodeRep<10,pztopology::TPZTetrahedron> {
		
	public:
		
		enum {NNodes = 10};
		
		bool IsLinearMapping() const {
			return false;
		}
		
		TPZQuadraticTetra(TPZVec<int> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(nodeindexes) {
		}
		
		TPZQuadraticTetra() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>() {
		}
		
		TPZQuadraticTetra(const TPZQuadraticTetra &cp,std::map<int,int> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp,gl2lcNdMap) {
		}
		
		TPZQuadraticTetra(const TPZQuadraticTetra &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp) {
		}
		
		TPZQuadraticTetra(const TPZQuadraticTetra &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp) {
		}
		
		virtual	~TPZQuadraticTetra();
		/**
		 * @brief Returns the type name of the element
		 */
		static std::string TypeName() {
			return "Tetra";
		}
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, 
		 * a side and a boundary condition number
		 */
		static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
		
		static void Shape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi);
		
		static void X(TPZFMatrix &coord, TPZVec<REAL> &loc,TPZVec<REAL> &result);
		
		static void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &param, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);
		
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
