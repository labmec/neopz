/**
 * @file
 * @brief Contains the TPZQuadraticCube class which defines a cube geometric element with quadratic map.
 */
 
#ifndef TPZQUADRATICCUBE_H
#define TPZQUADRATICCUBE_H


#include "TPZGeoCube.h"
#include "pzgeoel.h"
#include "pznoderep.h"

#include <iostream>

/**
 * @author Paulo Cesar de Alvarenga Lucci (Caju)
 * @since 2011
 * @brief Defines a cube geometric element with quadratic map. \ref geometry "Geometry"
 * @ingroup geometry
 */

class TPZQuadraticCube : public pzgeom::TPZNodeRep<20,pztopology::TPZCube> {
	
public:
	/** @brief Number of nodes (3 by edge) */
	enum {NNodes = 20};
    
    //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
    
	/** @brief It is not linear mapping, is quadratic */
	bool IsLinearMapping() const {
		return false;
	}
	/** @brief Constructor from node indexes */
	TPZQuadraticCube(TPZVec<long> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZCube>(nodeindexes)
	{
	}
	/** @brief Default constructor */
	TPZQuadraticCube() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZCube>()
	{
	}
	/** @brief Copy constructor from node map */
	TPZQuadraticCube(const TPZQuadraticCube &cp,std::map<long,long> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZCube>(cp,gl2lcNdMap)
	{
	}
	/** @brief Copy constructor */
	TPZQuadraticCube(const TPZQuadraticCube &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZCube>(cp)
	{
	}
	/** @brief Copy constructor */
	TPZQuadraticCube(const TPZQuadraticCube &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZCube>(cp)
	{
	}
	
	/** @brief Returns the type name of the element */
	static std::string TypeName() { return "Hexa";} 
	
    /* brief compute the coordinate of a point given in parameter space */
    void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
    {
        TPZFNMatrix<3*NNodes> coord(3,NNodes);
        CornerCoordinates(gel, coord);
        X(coord,loc,result);
    }
    
    /* @brief compute the jacobian of the map between the master element and deformed element */
    void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
    {
        TPZFNMatrix<3*NNodes> coord(3,NNodes);
        CornerCoordinates(gel, coord);
        Jacobian(coord, param, jacobian, axes, detjac, jacinv);
    }
    
	static void Shape(TPZVec<REAL> &x,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
	
	static void X(TPZFMatrix<REAL> &coord, TPZVec<REAL> &par, TPZVec<REAL> &result);
	
	static void Jacobian(TPZFMatrix<REAL> &coord, TPZVec<REAL> &par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv);
	
	
	/** @brief Creates a geometric element according to the type of the father element */
	static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
									  TPZVec<long>& nodeindexes,
									  int matid, long& index);
	
	TPZGeoEl *CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);	
};

#endif
