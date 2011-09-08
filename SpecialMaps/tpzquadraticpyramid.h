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
	
	bool IsLinearMapping() const {
		return false;
	}
	
	TPZQuadraticPyramid(TPZVec<int> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPyramid>(nodeindexes)
	{
	}
	
	TPZQuadraticPyramid() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPyramid>()
	{
	}
	
	TPZQuadraticPyramid(const TPZQuadraticPyramid &cp,std::map<int,int> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPyramid>(cp,gl2lcNdMap)
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
	
	static void Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi);
	
	static void X(TPZFMatrix &coord, TPZVec<REAL> &par, TPZVec<REAL> &result);
	
	static void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);
	
	
	/** @brief Creates a geometric element according to the type of the father element */
	static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
									  TPZVec<int>& nodeindexes,
									  int matid, int& index);
	
	TPZGeoEl *CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);	
};

};

#endif
