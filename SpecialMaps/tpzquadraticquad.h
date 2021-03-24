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
        typedef pztopology::TPZQuadrilateral Top;
		/** @brief Number of nodes */
		enum {NNodes = 8};
                
                public:
int ClassId() const override;

        
        //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
        
		/** @brief Constructor from node indexes */
		TPZQuadraticQuad(TPZVec<int64_t> &nodeindexes) : 
        TPZRegisterClassId(&TPZQuadraticQuad::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(nodeindexes)
		{
		}
		/** @brief Default constructor */
		TPZQuadraticQuad() : TPZRegisterClassId(&TPZQuadraticQuad::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>()
		{
		}
		/** @brief Constructor over node map */
		TPZQuadraticQuad(const TPZQuadraticQuad &cp,std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZQuadraticQuad::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp,gl2lcNdMap)
		{
		}
		/** @brief Copy constructor */
		TPZQuadraticQuad(const TPZQuadraticQuad &cp) : TPZRegisterClassId(&TPZQuadraticQuad::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
		{
		}
		/** @brief Copy constructor */
		TPZQuadraticQuad(const TPZQuadraticQuad &cp, TPZGeoMesh &) :
        TPZRegisterClassId(&TPZQuadraticQuad::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
		{
		}
		
		/**
		 * @brief Returns the type name of the element
		 */
		static std::string TypeName() { return "QuadraticQuad";}
        
        static bool IsLinearMapping(int side)
        {
            return false;
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
		
		/* brief compute the coordinate of a point given in parameter space */
//        template<class T>
//        void X(const TPZGeoEl &gel,TPZVec<T> &loc,TPZVec<T> &result) const
//        {
//            TPZFNMatrix<3*NNodes> coord(3,NNodes);
//            CornerCoordinates(gel, coord);
//            X(coord,loc,result);
//        }
        
        /** @brief Compute gradient of x mapping from local parametric coordinates */
//        template<class T>
//        void GradX(const TPZGeoEl &gel, TPZVec<T> &loc, TPZFMatrix<T> &gradx) const
//        {
//            TPZFNMatrix<3*NNodes> coord(3,NNodes);
//            CornerCoordinates(gel, coord);
//            GradX(coord,loc,gradx);
//        }
        
        template<class T>
        static void TShape(const TPZVec<T> &param,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
        
        template<class T>
		static void X(const TPZFMatrix<REAL> &coord, TPZVec<T> &par, TPZVec<T> &result);
        
        /** @brief Compute gradient of X mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
		
		/**
		 * @brief Creates a geometric element according to the type of the father element
		 */
		// static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
		// 								  TPZVec<int64_t>& nodeindexes,
		// 								  int matid,
		// 								  int64_t& index);
        
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

};

};

#endif
