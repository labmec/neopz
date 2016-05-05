/**
 * @file
 * @brief Contains the TPZQuadraticPyramid class which defines a pyramid geometric element with quadratic map.
 */
#ifndef TPZQUADRATICPYRAMID_H
#define TPZQUADRATICPYRAMID_H


#include "pzgeopyramid.h"
#include "pzgeoel.h"
#include "pznoderep.h"

#include <iostream>

/**
 * @author Paulo Cesar de Alvarenga Lucci (Caju)
 * @since 2011
 * @brief Defines a pyramid geometric element with quadratic map. \ref geometry "Geometry"
 * @ingroup geometry
 */

namespace pzgeom {
    
    class TPZQuadraticPyramid : public pzgeom::TPZNodeRep<13,pztopology::TPZPyramid> {
        
    public:
        
        enum {NNodes = 13};
        
        //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
        
        static bool IsLinearMapping(int side)
        {
            return false;
        }
        
        TPZQuadraticPyramid(TPZVec<long> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPyramid>(nodeindexes)
        {
        }
        
        TPZQuadraticPyramid() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPyramid>()
        {
        }
        
        TPZQuadraticPyramid(const TPZQuadraticPyramid &cp,std::map<long,long> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPyramid>(cp,gl2lcNdMap)
        {
        }
        
        TPZQuadraticPyramid(const TPZQuadraticPyramid &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPyramid>(cp)
        {
        }
        
        TPZQuadraticPyramid(const TPZQuadraticPyramid &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPyramid>(cp)
        {
        }
        
        /** @brief Returns the type name of the element */
        static std::string TypeName() { return "Pyramid";} 
        
        /** @brief Compute the shape being used to construct the X mapping from local parametric coordinates  */
        static void Shape(TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
            TShape(loc, phi, dphi);
        }
        
        /** @brief compute the coordinate of a point given in parameter space */
        template<class T>
        void X(const TPZGeoEl &gel,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
        
        /** @brief Compute gradient of x mapping from local parametric coordinates */
        template<class T>
        void GradX(const TPZGeoEl &gel, TPZVec<T> &loc, TPZFMatrix<T> &gradx) const
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
            
            GradX(nodes,loc,gradx);
        }
        
        template<class T>
        static void TShape(TPZVec<T> &param,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
		
        /* @brief compute the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
        }
        
        template<class T>
        static void X(TPZFMatrix<REAL> &coord, TPZVec<T> &par, TPZVec< T > &result);
        
        /** @brief Compute gradient of X mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(TPZFMatrix<T> &nodes,TPZVec<T> &par, TPZFMatrix<T> &gradx);
        
        static void Jacobian(TPZFMatrix<REAL> &coord, TPZVec<REAL> &par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv);
        
        
        /** @brief Creates a geometric element according to the type of the father element */
        static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                          TPZVec<long>& nodeindexes,
                                          int matid, long& index);
        
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

        TPZGeoEl *CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);	
    };
    
};

#endif
