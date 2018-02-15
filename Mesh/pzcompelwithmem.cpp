/**
 * @file
 * @brief Contains the declaration of the TPZCompElWithMem class, it is as TPZCompEl with enable material memory feature.
 */

#include "pzcompelwithmem.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"
#include "pzshapepiram.h"

template<class TBASE>
TPZCompElWithMem<TBASE>::~TPZCompElWithMem() {
	SetFreeIntPtIndices();
}


#ifndef BORLAND
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePoint> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeLinear> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeTriang> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeQuad> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeCube> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeTetra> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePrism> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePiram> >>;

#include "pzgeopoint.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "pzgeopyramid.h"
#include "TPZMultiphysicsInterfaceEl.h"


template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoCube> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsInterfaceElement>>;

#endif
/*
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoCube>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid>;*/