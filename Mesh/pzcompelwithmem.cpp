/**
 * @file
 * @brief Contains the declaration of the TPZCompElWithMem class, it is as TPZCompEl with enable material memory feature.
 */

#include "pzcompelwithmem.h"

template<class TBASE>
TPZCompElWithMem<TBASE>::~TPZCompElWithMem() {
	SetFreeIntPtIndices();
}


#ifndef BORLAND
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePoint> >, TPZCOMPELWITHMEMPOINTID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeLinear> >, TPZCOMPELWITHMEMLINEARID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeTriang> >, TPZCOMPELWITHMEMTRIANGID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeQuad> >, TPZCOMPELWITHMEMQUADID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeCube> >, TPZCOMPELWITHMEMCUBEID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeTetra> >, TPZCOMPELWITHMEMTETRAID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePrism> >, TPZCOMPELWITHMEMPRISMID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePiram> >, TPZCOMPELWITHMEMPIRAMID>;

#include "pzgeopoint.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "pzgeopyramid.h"
#include "TPZMultiphysicsInterfaceEl.h"


template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint> >, TPZMPCOMPELWITHMEMPOINTID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear> >, TPZMPCOMPELWITHMEMLINEARID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle> >, TPZMPCOMPELWITHMEMTRIANGID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad> >, TPZMPCOMPELWITHMEMQUADID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoCube> >, TPZMPCOMPELWITHMEMCUBEID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra> >, TPZMPCOMPELWITHMEMTETRAID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism> >, TPZMPCOMPELWITHMEMPRISMID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid> >, TPZMPCOMPELWITHMEMPIRAMID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsInterfaceElement>, TPZMPCOMPELWITHMEMINTERFACE>;

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