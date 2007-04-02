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

#include "pzgeoelrefless.h.h"

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

TPZCompEl *CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoPoint,TPZShapePoint>(mesh,gel,index);
}
TPZCompEl *CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoLinear,TPZShapeLinear>(mesh,gel,index);
}
TPZCompEl *CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoQuad,TPZShapeQuad>(mesh,gel,index);
}
TPZCompEl *CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoTriangle,TPZShapeTriang>(mesh,gel,index);
}
TPZCompEl *CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoCube,TPZShapeCube>(mesh,gel,index);
}
TPZCompEl *CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoPrism,TPZShapePrism>(mesh,gel,index);
}
TPZCompEl *CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoPyramid,TPZShapePiram>(mesh,gel,index);
}
TPZCompEl *CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra>(mesh,gel,index);
}


template<>
TPZCompEl *(*TPZGeoElRefLess<TPZShapePoint,TPZGeoPoint>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreatePointEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZShapeLinear,TPZGeoLinear>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateLinearEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZShapeQuad,TPZGeoQuad>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateQuadEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZShapeTriang,TPZGeoTriangle>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateTriangleEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZShapeCube,TPZGeoCube>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateCubeEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZShapePrism,TPZGeoPrism>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreatePrismEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZShapeTetra,TPZGeoTetrahedra>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateTetraEl;
template<>
TPZCompEl *(*TPZGeoElRefLess<TPZShapePiram,TPZGeoPyramid>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreatePyramEl;


template class TPZGeoElRefLess<TPZShapeCube,TPZGeoCube>;
template class TPZGeoElRefLess<TPZShapeLinear,TPZGeoLinear>;
template class TPZGeoElRefLess<TPZShapeQuad,TPZGeoQuad>;
template class TPZGeoElRefLess<TPZShapeTriang,TPZGeoTriangle>;
template class TPZGeoElRefLess<TPZShapePrism,TPZGeoPrism>;
template class TPZGeoElRefLess<TPZShapeTetra,TPZGeoTetrahedra>;
template class TPZGeoElRefLess<TPZShapePiram,TPZGeoPyramid>;
template class TPZGeoElRefLess<TPZShapePoint,TPZGeoPoint>;

int main_refless(){
//TPZGeoElRefLess <TPZShapeCube,TPZGeoCube> el1;
//TPZGeoElRefLess <TPZShapeLinear,TPZGeoLinear> el2;
//TPZGeoElRefLess <TPZShapeQuad,TPZGeoQuad> el3;
//TPZGeoElRefLess <TPZShapeTriang,TPZGeoTriangle> el4;
//TPZGeoElRefLess <TPZShapePrism,TPZGeoPrism> el5;
//TPZGeoElRefLess <TPZShapeTetra,TPZGeoTetrahedra> el6;
//TPZGeoElRefLess <TPZShapePiram,TPZGeoPyramid> el7;
//TPZGeoElRefLess <TPZShapePoint,TPZGeoPoint> el8;

return 0;
}







