/*
 *  Instantiate.cpp
 *  SubstructEigen
 *
 *  Created by Philippe Devloo on 30/11/08.
 *  Copyright 2008 UNICAMP. All rights reserved.
 *
 */

#include "Instantiate.h"

#include "tpzcurvedtriangle.h"
#include "pzgeoelrefless.h"
#include "tpzgeoelmapped.h"
#include "tpzgeoelrefpattern.h"

#include "pzgeoelrefless.h.h"
#include "tpzgeoelrefpattern.h.h"




///CreateGeoElement -> TPZCurvedTriangle
template< >
TPZGeoEl *TPZGeoElRefLess<TPZCurvedTriangle >::CreateGeoElement(MElementType type, TPZVec<int>& nodeindexes, int matid, int& index)
{
    TPZGeoMesh &mesh = *(this->Mesh());
    if(!&mesh) return 0;
    return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

#define TPZGEOELEMENTCURVEDTRIANGLEID 312
template<>
int TPZGeoElRefPattern<TPZCurvedTriangle>::ClassId() const {
	return TPZGEOELEMENTCURVEDTRIANGLEID;
}
template class 
TPZRestoreClass< TPZGeoElRefPattern<TPZCurvedTriangle>, TPZGEOELEMENTCURVEDTRIANGLEID>;

template<>
TPZCompEl *(*TPZGeoElRefLess<TPZCurvedTriangle>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateTriangleEl;

template class TPZGeoElRefLess<TPZCurvedTriangle>;


#include "pznoderep.h.h"
template class pzgeom::TPZNodeRep<3,TPZCurvedTriangle>;

