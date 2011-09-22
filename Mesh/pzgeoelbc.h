/**
 * @file
 * @brief Contains declaration of TPZGeoElBC class, it is a structure to help the construction of geometric element along a side.
 */
//$Id: pzgeoelbc.h,v 1.6 2006-10-17 01:38:03 phil Exp $

#ifndef PZGEOELBCH
#define PZGEOELBCH

#include <iostream>
#include "pzsave.h"


class TPZGeoMesh;
class TPZGeoEl;
class TPZGeoElSide;

/** 
 * @brief Structure to help the construction of geometric elements along side of a given geometric element. \ref geometry "Geometry"
 * @ingroup geometry
 */
struct TPZGeoElBC {
	
private:
	
	/** @brief Pointer to the geometric element created in the class constructor */
	TPZGeoEl *fCreatedElement;
	
public:
	
	/** @brief Creates a geometric element along side of el. 
	 *
	 * The new geometric element is inserted in mesh and a pointer to it is stored here. 
	 */
	TPZGeoElBC(TPZGeoEl *el,int side,int matid);
	
	/** @brief Creates a geometric element along side of el. 
	 *  
	 * The new geometric element is inserted in mesh and a pointer to it is stored here. 
	 */
	TPZGeoElBC(TPZGeoElSide &elside,int matid);
	
	/** @brief Recovers pointer to the geometric element created */
	TPZGeoEl * CreatedElement(){ return this->fCreatedElement; }
	
};

#endif
