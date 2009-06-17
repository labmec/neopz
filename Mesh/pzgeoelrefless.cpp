/***************************************************************************
                          pzgeoelrefless.cc  -  description
                             -------------------
    begin                : Fri Dec 12 2003
    copyright            : (C) 2003 by phil
    email                : phil@localhost
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "tpzint1point.h"
#include "pzintel.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
//#include "TPZRefPattern.h"
#include "pzvec.h"
#include "pzmanvector.h"
//#include "pzstack.h"



#include "pzelctemp.h"

using namespace pzgeom;
using namespace pzshape;
using namespace pztopology;

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"
#include "tpzcube.h"
#include "tpztetrahedron.h"
#include "tpzprism.h"
#include "tpzpyramid.h"

#ifndef WIN32
#include "pzgeoelrefless.h.h"
#endif

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


template<>
TPZCompEl *(*TPZGeoElRefLess<TPZGeoPoint>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreatePointEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZGeoLinear>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateLinearEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZGeoQuad>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateQuadEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZGeoTriangle>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateTriangleEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZGeoCube>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateCubeEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZGeoPrism>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreatePrismEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZGeoTetrahedra>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateTetraEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZGeoPyramid>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreatePyramEl;
//HDiv

 template<>
 void TPZGeoElRefLess<pzgeom::TPZGeoQuad>::VecHdiv(TPZFMatrix &coordinate, TPZFMatrix &normalvec,TPZVec<int> &sidevector )
 {
 pzgeom::TPZGeoQuad::VecHdiv(coordinate,normalvec,sidevector);
 }
 
template<>
void TPZGeoElRefLess<pzgeom::TPZGeoTriangle>::VecHdiv(TPZFMatrix &coordinate, TPZFMatrix &normalvec,TPZVec<int> &sidevector )
{
	pzgeom::TPZGeoTriangle::VecHdiv(coordinate,normalvec,sidevector);
}


template class TPZGeoElRefLess<TPZGeoCube>;
template class TPZGeoElRefLess<TPZGeoLinear>;
template class TPZGeoElRefLess<TPZGeoQuad>;
template class TPZGeoElRefLess<TPZGeoTriangle>;
template class TPZGeoElRefLess<TPZGeoPrism>;
template class TPZGeoElRefLess<TPZGeoTetrahedra>;
template class TPZGeoElRefLess<TPZGeoPyramid>;
template class TPZGeoElRefLess<TPZGeoPoint>;

int main_refless()
{
return 0;
}







