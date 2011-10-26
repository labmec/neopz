/*
 *  TPZMultiphysicsInterfaceEl.cpp
 *  PZ
 *
 *  Created by Agnaldo on 10/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZMultiphysicsInterfaceEl.h"


TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement() : TPZCompEl(),fLeftElSideVec(0), fRightElSideVec(0)
{
}

TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int &index) : TPZCompEl(mesh, ref, index),fLeftElSideVec(0), fRightElSideVec(0)
{
}

TPZMultiphysicsInterfaceElement::~TPZMultiphysicsInterfaceElement(){
}

void TPZMultiphysicsInterfaceElement::ComputeSideTransform(TPZManVector<TPZCompElSide> &Neighbor, TPZManVector<TPZTransform> &transf)
{
	std::cout<<"Falta implementar ....."<<std::endl;
}//ComputeSideTransform