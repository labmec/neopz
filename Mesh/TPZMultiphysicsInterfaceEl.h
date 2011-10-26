/*
 *  TPZMultiphysicsInterfaceEl.h
 *  PZ
 *
 *  Created by Agnaldo on 10/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TPZMULTIPHYSICSINTERFACEELH
#define TPZMULTIPHYSICSINTERFACEELH 

#include <iostream>

#include "pzcompel.h"


/**
 @brief Computes the contribution over an interface between two discontinuous elements. \ref CompElement "Computational Element"
  */
class TPZMultiphysicsInterfaceElement : public TPZCompEl {

protected:

	/**Element vector the left of the normal a interface*/
	TPZManVector<TPZCompElSide*, 10> 	fLeftElSideVec;
		
	/**Element vector the right of the normal a interface*/
	TPZManVector<TPZCompElSide*, 10> 	fRightElSideVec;
	
public:
	TPZMultiphysicsInterfaceElement();	
		
	TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int &index);
		
	virtual ~TPZMultiphysicsInterfaceElement();
	
	
	/**
	 * @brief Compute the transform of a paramenter point in the multiphysic interface element to a parameter point in the neighbor super element
	 ** @param Neighbor [in] may be this->LeftElementSide() or this->RightElementSide()
	 * @param transf: vector of Transforms 
	 **/
	void ComputeSideTransform(TPZManVector<TPZCompElSide> &Neighbor, TPZManVector<TPZTransform> &transf);
	
	/**
	 * @brief Maps qsi coordinate at this master element to qsi coordinate at neighbor master element.
	 * @param Neighbor [in] may be this->LeftElementSide() or this->RightElementSide()
	 * @param qsi [in] is the point at this element master
	 * @param NeighIntPoint [out] is the point at neighbor element master. X[qsi] is equal to X[NeighIntPoint]
	 */
	void MapQsi(TPZManVector<TPZCompElSide> &Neighbor, TPZVec<REAL> &qsi, TPZVec<REAL> &NeighIntPoint);
			
};

#endif