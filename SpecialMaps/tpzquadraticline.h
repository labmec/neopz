/**
 * @file
 * @brief Contains the TPZQuadraticLine class which defines a linear geometric element with quadratic map.
 */
#ifndef TPZQUADRATICLINE_H
#define TPZQUADRATICLINE_H

#include "TPZGeoLinear.h"
#include "pzgeoel.h"
#include "pznoderep.h"

#include <iostream>

namespace pzgeom {

	/**
	 * @ingroup geometry
	 * @brief Defines a linear geometric element with quadratic map. \ref geometry "Geometry"
	 * @author Paulo Cesar de Alvarenga Lucci
	 * @since 2007
	 */
	class TPZQuadraticLine : public pzgeom::TPZNodeRep<3,pztopology::TPZLine> 
    {		
        public:
        typedef pztopology::TPZLine Top;
		enum {NNodes = 3};
                
                public:
int ClassId() const override;

        
        //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);

		
		TPZQuadraticLine(TPZVec<int64_t> &nodeindexes) :
        TPZRegisterClassId(&TPZQuadraticLine::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes)
		{
		}
		
		TPZQuadraticLine() : TPZRegisterClassId(&TPZQuadraticLine::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>()
		{
		}
		
		TPZQuadraticLine(const TPZQuadraticLine &cp,std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZQuadraticLine::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap)
		{
		}
		
		TPZQuadraticLine(const TPZQuadraticLine &cp) : TPZRegisterClassId(&TPZQuadraticLine::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp)
		{
		}
		
		TPZQuadraticLine(const TPZQuadraticLine &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZQuadraticLine::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp)
		{
		}
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "QuadraticLine";}
        
        static bool IsLinearMapping(int side)
        {
            return false;
        }
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, \n
		 * a side and a boundary condition number
		 */
		// static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
		
		/**
		 * @brief Creates a geometric element according to the type of the father element
		 */
		// static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
		// 								  TPZVec<int64_t>& nodeindexes,
		// 								  int matid,
		// 								  int64_t& index);
		
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
        
        
        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        template<class T>
        static void TShape(const TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
        
        template<class T>
		static void X(const TPZFMatrix<REAL> &coord, TPZVec<T> &par, TPZVec<T> &result);
        
        /** @brief Compute gradient of X mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
        
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

	};
    
    
};

#endif
