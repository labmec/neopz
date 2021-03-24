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
        typedef pztopology::TPZPyramid Top;
        enum {NNodes = 13};
        
        public:
int ClassId() const override;

        
        //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
        
        static bool IsLinearMapping(int side)
        {
            return false;
        }
        
        TPZQuadraticPyramid(TPZVec<int64_t> &nodeindexes) : 
        TPZRegisterClassId(&TPZQuadraticPyramid::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZPyramid>(nodeindexes)
        {
        }
        
        TPZQuadraticPyramid() : TPZRegisterClassId(&TPZQuadraticPyramid::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZPyramid>()
        {
        }
        
        TPZQuadraticPyramid(const TPZQuadraticPyramid &cp,std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZQuadraticPyramid::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZPyramid>(cp,gl2lcNdMap)
        {
        }
        
        TPZQuadraticPyramid(const TPZQuadraticPyramid &cp) : 
        TPZRegisterClassId(&TPZQuadraticPyramid::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZPyramid>(cp)
        {
        }
        
        TPZQuadraticPyramid(const TPZQuadraticPyramid &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZQuadraticPyramid::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZPyramid>(cp)
        {
        }
        
        /** @brief Returns the type name of the element */
        static std::string TypeName() { return "QuadraticPyramid";} 
        
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
        static void GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &par, TPZFMatrix<T> &gradx);

        /** @brief Creates a geometric element according to the type of the father element */
        // static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
        //                                   TPZVec<int64_t>& nodeindexes,
        //                                   int matid, int64_t& index);
        
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

        // TPZGeoEl *CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);	
    };
    
};

#endif
