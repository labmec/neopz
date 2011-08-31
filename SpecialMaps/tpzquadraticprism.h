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

class TPZQuadraticPrism : public pzgeom::TPZNodeRep<15,pztopology::TPZPrism> {
	
public:
	
	enum {NNodes = 15};
	
	bool IsLinearMapping() const {
		return false;
	}
	
	TPZQuadraticPrism(TPZVec<int> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPrism>(nodeindexes)
	{
	}
	
	TPZQuadraticPrism() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPrism>()
	{
	}
	
	TPZQuadraticPrism(const TPZQuadraticPrism &cp,std::map<int,int> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZPrism>(cp,gl2lcNdMap)
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
