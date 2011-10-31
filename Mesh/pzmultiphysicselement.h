/*
 *  tpzmultiphysicselement.h
 *  PZ
 *
 *  Created by Agnaldo on 9/30/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PZMULTIPHYSICSELEMENTH
#define PZMULTIPHYSICSELEMENTH 

#include <iostream>

#include "pzcompel.h"


class TPZMultiphysicsElement : public TPZCompEl {
	
	
public:
	TPZMultiphysicsElement() : TPZCompEl()
	{
	}
	
	TPZMultiphysicsElement(TPZCompMesh &mesh, TPZGeoEl *ref, int &index) : TPZCompEl(mesh, ref, index)
	{
	}
	
	virtual ~TPZMultiphysicsElement()
	{
	}
	
	virtual void AddElement(TPZCompEl *cel, int mesh) = 0;
	
	virtual TPZCompEl *ReferredElement(int mesh) = 0;
	
	virtual void SetConnectIndexes(TPZVec<int> &indexes) = 0;
	
	virtual void AffineTransform(TPZManVector<TPZTransform> &trVec) const = 0;
	
	virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef) = 0;
	
	
	virtual void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension)=0;
	
};

#endif