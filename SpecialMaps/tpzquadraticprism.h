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
	
	enum {NNodes = 15};
    
    //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
	

	TPZQuadraticPrism(TPZVec<long> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPrism>(nodeindexes)
	{
	}
	
	TPZQuadraticPrism() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPrism>()
	{
	}
	
	TPZQuadraticPrism(const TPZQuadraticPrism &cp,std::map<long,long> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPrism>(cp,gl2lcNdMap)
	{
	}
	
	TPZQuadraticPrism(const TPZQuadraticPrism &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPrism>(cp)
	{
	}
	
	TPZQuadraticPrism(const TPZQuadraticPrism &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPrism>(cp)
	{
	}
	
	/** @brief Returns the type name of the element */
	static std::string TypeName() { return "Prism";}
    
    static bool IsLinearMapping(int side)
    {
        return false;
    }
	
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
	
	
	/** @brief Creates a geometric element according to the type of the father element */
	static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
									  TPZVec<long>& nodeindexes,
									  int matid, long& index);
	
    static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

	TPZGeoEl *CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);	
};

};

#endif
