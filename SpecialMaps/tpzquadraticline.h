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
		
		enum {NNodes = 3};
        
        //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);

		
		TPZQuadraticLine(TPZVec<long> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes)
		{
		}
		
		TPZQuadraticLine() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>()
		{
		}
		
		TPZQuadraticLine(const TPZQuadraticLine &cp,std::map<long,long> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap)
		{
		}
		
		TPZQuadraticLine(const TPZQuadraticLine &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp)
		{
		}
		
		TPZQuadraticLine(const TPZQuadraticLine &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp)
		{
		}
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Line";}
        
        static bool IsLinearMapping(int side)
        {
            return false;
        }
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, \n
		 * a side and a boundary condition number
		 */
		static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
		
		/**
		 * @brief Creates a geometric element according to the type of the father element
		 */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<long>& nodeindexes,
										  int matid,
										  long& index);
		
        /** @brief Compute the shape being used to construct the X mapping from local parametric coordinates  */
        static void Shape(TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
            TShape(loc, phi, dphi);
        }

        
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
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            int nrow = coord.Rows();
            int ncol = coord.Cols();
            TPZFMatrix<T> nodes(nrow,ncol);
            for(int i = 0; i < nrow; i++)
            {
                for(int j = 0; j < ncol; j++)
                {
                    nodes(i,j) = coord(i,j);
                }
            }
            
            GradX(nodes,par,gradx);
        }
        
        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        template<class T>
        static void TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
		
        /* @brief compute the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
        }
        
        template<class T>
		static void X(TPZFMatrix<REAL> &coord, TPZVec<T> &par, TPZVec<T> &result);
        
        /** @brief Compute gradient of X mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<T> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
		
		static void Jacobian(TPZFMatrix<REAL> &coord, TPZVec<REAL> &par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv);
        
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

	};
    
    
};

#endif
