
#include "TPZGeoElement.h.h"


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
#include "pzstream.h"
#include "pzmeshid.h"
//#include "pzstack.h"



using namespace pzgeom;
using namespace pzrefine;
using namespace pzshape;

template class TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>;
template class TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>;
template class TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>;
template class TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>;
template class TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>;
template class TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>;
template class TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>;
template class TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>;

#ifndef BORLAND

/** ClassId method for each instantiation followed by the registration of the class in the TPZRestoreClass */
template<>
int TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>::ClassId() const {
  return TPZFGEOELEMENTPOINTID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>, TPZFGEOELEMENTPOINTID>;

template<>
int TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::ClassId() const {
  return TPZFGEOELEMENTLINEARID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>, TPZFGEOELEMENTLINEARID>;

template<>
int TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::ClassId() const {
  return TPZFGEOELEMENTQUADID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>, TPZFGEOELEMENTQUADID>;

template<>
int TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::ClassId() const {
  return TPZFGEOELEMENTRIANGLEID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>, TPZFGEOELEMENTRIANGLEID>;

template<>
int TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::ClassId() const {
  return TPZFGEOELEMENTCUBEID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>, TPZFGEOELEMENTCUBEID>;

template<>
int TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::ClassId() const {
  return TPZFGEOELEMENTPRISMID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>, TPZFGEOELEMENTPRISMID>;

template<>
int TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::ClassId() const {
  return TPZFGEOELEMENTTETRAID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>, TPZFGEOELEMENTTETRAID>;

template<>
int TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::ClassId() const {
  return TPZFGEOELEMENTPYRAMID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>, TPZFGEOELEMENTPYRAMID>;

#else

int
TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>::ClassId() const {
  return TPZFGEOELEMENTPOINTID;
}

int
TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::ClassId() const {
  return TPZFGEOELEMENTLINEARID;
}

int
TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::ClassId() const {
  return TPZFGEOELEMENTQUADID;
}

int
TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::ClassId() const {
  return TPZFGEOELEMENTRIANGLEID;
}

int
TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::ClassId() const {
  return TPZFGEOELEMENTCUBEID;
}

int
TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::ClassId() const {
  return TPZFGEOELEMENTPRISMID;
}

int
TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::ClassId() const {
  return TPZFGEOELEMENTTETRAID;
}

int
TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::ClassId() const {
  return TPZFGEOELEMENTPYRAMID;
}
#endif
