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
        typedef pztopology::TPZTetrahedron Top;
		/** @brief Number of nodes */
		enum {NNodes = 10};
                
                public:
int ClassId() const override;

        
        //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
        
		/** @brief It is quadratic mapping */
        static bool IsLinearMapping(int side)
        {
            return false;
        }
		/** @brief Constructor for node indexes given */
		TPZQuadraticTetra(TPZVec<int64_t> &nodeindexes) : 
        TPZRegisterClassId(&TPZQuadraticTetra::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(nodeindexes) {
		}
		/** @brief Default constructor */
		TPZQuadraticTetra() : TPZRegisterClassId(&TPZQuadraticTetra::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>() {
		}
		/** @brief Constructor for node map given */
		TPZQuadraticTetra(const TPZQuadraticTetra &cp,std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZQuadraticTetra::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp,gl2lcNdMap) {
		}
		/** @brief Copy constructor */
		TPZQuadraticTetra(const TPZQuadraticTetra &cp) : 
        TPZRegisterClassId(&TPZQuadraticTetra::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp) {
		}
		/** @brief Copy constructor */		
		TPZQuadraticTetra(const TPZQuadraticTetra &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZQuadraticTetra::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp) {
		}
		/** @brief Destructor */
		virtual	~TPZQuadraticTetra();
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() {
			return "QuadraticTetra";
		}
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, 
		 * a side and a boundary condition number
		 */
		// static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
		
        /** @brief Compute the shape being used to construct the X mapping from local parametric coordinates  */
        static void Shape(TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
            TShape(loc, phi, dphi);
        }
		
        
        template<class T>
        static void TShape(const TPZVec<T> &param,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
        
        template<class T>
        static void X(const TPZFMatrix<REAL> &coord, TPZVec<T> &par, TPZVec< T > &result);
        
        /** @brief Compute gradient of X mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
		
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		// static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
		// 								  TPZVec<int64_t>& nodeindexes,
		// 								  int matid,
		// 								  int64_t& index);	
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
    };

};

#endif
