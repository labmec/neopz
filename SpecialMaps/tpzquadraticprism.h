/**
 * @file
 * @brief Contains the TPZQuadraticPrism class which defines a prism geometric element with quadratic map.
 */
#ifndef TPZQUADRATICPRISM_H
#define TPZQUADRATICPRISM_H


#include "pzgeoprism.h"
#include "pzgeoel.h"
#include "pznoderep.h"

#include <iostream>

/**
 * @author Paulo Cesar de Alvarenga Lucci (Caju)
 * @since 2011
 * @brief Defines a prism geometric element with quadratic map. \ref geometry "Geometry"
 * @ingroup geometry
 */

namespace pzgeom {
    
class TPZQuadraticPrism : public pzgeom::TPZNodeRep<15,pztopology::TPZPrism> {
	
public:
    typedef pztopology::TPZPrism Top;
	enum {NNodes = 15};
        
int ClassId() const override;

    //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
	

	TPZQuadraticPrism(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZQuadraticPrism::ClassId),
    pzgeom::TPZNodeRep<NNodes,pztopology::TPZPrism>(nodeindexes)
	{
	}
	
	TPZQuadraticPrism() : TPZRegisterClassId(&TPZQuadraticPrism::ClassId),
    pzgeom::TPZNodeRep<NNodes,pztopology::TPZPrism>()
	{
	}
	
	TPZQuadraticPrism(const TPZQuadraticPrism &cp,std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZQuadraticPrism::ClassId), pzgeom::TPZNodeRep<NNodes,pztopology::TPZPrism>(cp,gl2lcNdMap)
	{
	}
	
	TPZQuadraticPrism(const TPZQuadraticPrism &cp) : TPZRegisterClassId(&TPZQuadraticPrism::ClassId),pzgeom::TPZNodeRep<NNodes,pztopology::TPZPrism>(cp)
	{
	}
	
	TPZQuadraticPrism(const TPZQuadraticPrism &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZQuadraticPrism::ClassId),pzgeom::TPZNodeRep<NNodes,pztopology::TPZPrism>(cp)
	{
	}
	
	/** @brief Returns the type name of the element */
	static std::string TypeName() { return "QuadraticPrism";}
    
    static bool IsLinearMapping(int side)
    {
        return false;
    }
	
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
	// 								  TPZVec<int64_t>& nodeindexes,
	// 								  int matid, int64_t& index);
	
    static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

	// TPZGeoEl *CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);	
};

};

#endif
