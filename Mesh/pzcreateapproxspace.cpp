/**
 * @file
 * @brief Contains the implementation of the functions to creates computational elements.
 */
/*
 *  pzcreatecontinuous.cpp
 *  NeoPZ
 *  Created by Philippe Devloo on 20/11/09.
 *  Copyright 2009 UNICAMP. All rights reserved.
 */

#include "pzcreateapproxspace.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzelctemp.h"

#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapepiram.h"
#include "pzshapepoint.h"
#include "pzshapeprism.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapetriang.h"

#include "pzgeopoint.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"
#include "pzgeoquad.h"
#include "pzgeotetrahedra.h"
#include "pzgeotriangle.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"

using namespace pzshape;

TPZCompEl *CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapePoint>(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapeLinear>(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapeQuad>(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapeTriang>(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapeCube>(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapePrism>(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapePiram>(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapeTetra>(mesh,gel,index);
	return NULL;
}

void TPZCompMesh::SetAllCreateFunctionsDiscontinuous(){
	
	pzgeom::TPZGeoPoint::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoLinear::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoQuad::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoTriangle::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoPrism::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoTetrahedra::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoPyramid::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoCube::fp = TPZCompElDisc::CreateDisc;
}


void TPZCompMesh::SetAllCreateFunctionsContinuous(){
	pzgeom::TPZGeoPoint::fp =  CreatePointEl;
	pzgeom::TPZGeoLinear::fp =  CreateLinearEl;
	pzgeom::TPZGeoQuad::fp = CreateQuadEl;
	pzgeom::TPZGeoTriangle::fp =  CreateTriangleEl;
	pzgeom::TPZGeoPrism::fp = CreatePrismEl;
	pzgeom::TPZGeoTetrahedra::fp = CreateTetraEl;
	pzgeom::TPZGeoPyramid::fp = CreatePyramEl;
	pzgeom::TPZGeoCube::fp = CreateCubeEl;
	
	
}

#include "pzelchdiv.h"

void TPZCompMesh::SetAllCreateFunctionsHDiv(){
	
    pzgeom::TPZGeoPoint::fp =  CreateHDivPointEl;
	pzgeom::TPZGeoLinear::fp =  CreateHDivLinearEl;
	pzgeom::TPZGeoQuad::fp = CreateHDivQuadEl;
	pzgeom::TPZGeoTriangle::fp =  CreateHDivTriangleEl;
	pzgeom::TPZGeoPrism::fp = CreateHDivPrismEl;
	pzgeom::TPZGeoTetrahedra::fp = CreateHDivTetraEl;
	pzgeom::TPZGeoPyramid::fp = CreateHDivPyramEl;
	pzgeom::TPZGeoCube::fp = CreateHDivCubeEl;
}


#include "pzreferredcompel.h"
#include "pzelctemp.h"
void TPZCompMesh::SetAllCreateFunctionsDiscontinuousReferred(){
	
	pzgeom::TPZGeoPoint::fp =  CreateReferredDisc;
	pzgeom::TPZGeoLinear::fp =  CreateReferredDisc;
	pzgeom::TPZGeoQuad::fp = CreateReferredDisc;
	pzgeom::TPZGeoTriangle::fp =  CreateReferredDisc;
	pzgeom::TPZGeoPrism::fp = CreateReferredDisc;
	pzgeom::TPZGeoTetrahedra::fp = CreateReferredDisc;
	pzgeom::TPZGeoPyramid::fp = CreateReferredDisc;
	pzgeom::TPZGeoCube::fp = CreateReferredDisc;
}

void TPZCompMesh::SetAllCreateFunctionsContinuousReferred(){
	
	pzgeom::TPZGeoPoint::fp =  CreateReferredPointEl;
	pzgeom::TPZGeoLinear::fp =  CreateReferredLinearEl;
	pzgeom::TPZGeoQuad::fp = CreateReferredQuadEl;
	pzgeom::TPZGeoTriangle::fp =  CreateReferredTriangleEl;
	pzgeom::TPZGeoPrism::fp = CreateReferredPrismEl;
	pzgeom::TPZGeoTetrahedra::fp = CreateReferredTetraEl;
	pzgeom::TPZGeoPyramid::fp = CreateReferredPyramEl;
	pzgeom::TPZGeoCube::fp = CreateReferredCubeEl;
	
}

#include "pzmultiphysicscompel.h"
void TPZCompMesh::SetAllCreateFunctionsMultiphysicElem(){
	
	pzgeom::TPZGeoPoint::fp =  CreateMultiphysicsPointEl;
	pzgeom::TPZGeoLinear::fp =  CreateMultiphysicsLinearEl;
	pzgeom::TPZGeoTriangle::fp =  CreateMultiphysicsTriangleEl;
	pzgeom::TPZGeoQuad::fp = CreateMultiphysicsQuadEl;
	pzgeom::TPZGeoCube::fp = CreateMultiphysicsCubeEl;
	pzgeom::TPZGeoPrism::fp = CreateMultiphysicsPrismEl;
	pzgeom::TPZGeoTetrahedra::fp = CreateMultiphysicsTetraEl;
	pzgeom::TPZGeoPyramid::fp = CreateMultiphysicsPyramEl;
}


