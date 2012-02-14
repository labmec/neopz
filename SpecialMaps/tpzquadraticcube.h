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
	
	enum {NNodes = 20};
	
	bool IsLinearMapping() const {
		return false;
	}
	
	TPZQuadraticCube(TPZVec<int> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZCube>(nodeindexes)
	{
	}
	
	TPZQuadraticCube() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZCube>()
	{
	}
	
	TPZQuadraticCube(const TPZQuadraticCube &cp,std::map<int,int> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZCube>(cp,gl2lcNdMap)
	{
	}
	
	TPZQuadraticCube(const TPZQuadraticCube &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZCube>(cp)
	{
	}
	
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
    void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) const
    {
        TPZFNMatrix<3*NNodes> coord(3,NNodes);
        CornerCoordinates(gel, coord);
        Jacobian(coord, param, jacobian, axes, detjac, jacinv);
    }
    
	static void Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi);
	
	static void X(TPZFMatrix &coord, TPZVec<REAL> &par, TPZVec<REAL> &result);
	
	static void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);
	
	
	/** @brief Creates a geometric element according to the type of the father element */
	static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
									  TPZVec<int>& nodeindexes,
									  int matid, int& index);
	
	TPZGeoEl *CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);	
};

#endif
