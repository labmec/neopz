/**
 * @file
 * @brief Contains the instantiation of the TPZGeoElement classes from template.
 */

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


#ifndef BORLAND

/** ClassId method for each instantiation followed by the registration of the class in the TPZRestoreClass */
template<>
int TPZGeoElement<TPZGeoPoint,TPZRefPoint>::ClassId() const {
	return TPZFGEOELEMENTPOINTID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZGeoPoint,TPZRefPoint>, TPZFGEOELEMENTPOINTID>;

template<>
int TPZGeoElement<TPZGeoLinear,TPZRefLinear>::ClassId() const {
	return TPZFGEOELEMENTLINEARID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZGeoLinear,TPZRefLinear>, TPZFGEOELEMENTLINEARID>;

template<>
int TPZGeoElement<TPZGeoQuad,TPZRefQuad>::ClassId() const {
	return TPZFGEOELEMENTQUADID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZGeoQuad,TPZRefQuad>, TPZFGEOELEMENTQUADID>;

template<>
int TPZGeoElement<TPZGeoTriangle,TPZRefTriangle>::ClassId() const {
	return TPZFGEOELEMENTRIANGLEID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZGeoTriangle,TPZRefTriangle>, TPZFGEOELEMENTRIANGLEID>;

template<>
int TPZGeoElement<TPZGeoCube,TPZRefCube>::ClassId() const {
	return TPZFGEOELEMENTCUBEID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZGeoCube,TPZRefCube>, TPZFGEOELEMENTCUBEID>;

template<>
int TPZGeoElement<TPZGeoPrism,TPZRefPrism>::ClassId() const {
	return TPZFGEOELEMENTPRISMID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZGeoPrism,TPZRefPrism>, TPZFGEOELEMENTPRISMID>;

template<>
int TPZGeoElement<TPZGeoTetrahedra,TPZRefTetrahedra>::ClassId() const {
	return TPZFGEOELEMENTTETRAID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZGeoTetrahedra,TPZRefTetrahedra>, TPZFGEOELEMENTTETRAID>;

template<>
int TPZGeoElement<TPZGeoPyramid,TPZRefPyramid>::ClassId() const {
	return TPZFGEOELEMENTPYRAMID;
}
template class
TPZRestoreClass< TPZGeoElement<TPZGeoPyramid,TPZRefPyramid>, TPZFGEOELEMENTPYRAMID>;

#else

template<>
int
TPZGeoElement<TPZGeoPoint,TPZRefPoint>::ClassId() const {
	return TPZFGEOELEMENTPOINTID;
}

template<>
int
TPZGeoElement<TPZGeoLinear,TPZRefLinear>::ClassId() const {
	return TPZFGEOELEMENTLINEARID;
}

template<>
int
TPZGeoElement<TPZGeoQuad,TPZRefQuad>::ClassId() const {
	return TPZFGEOELEMENTQUADID;
}

template<>
int
TPZGeoElement<TPZGeoTriangle,TPZRefTriangle>::ClassId() const {
	return TPZFGEOELEMENTRIANGLEID;
}

template<>
int
TPZGeoElement<TPZGeoCube,TPZRefCube>::ClassId() const {
	return TPZFGEOELEMENTCUBEID;
}

template<>
int
TPZGeoElement<TPZGeoPrism,TPZRefPrism>::ClassId() const {
	return TPZFGEOELEMENTPRISMID;
}

template<>
int
TPZGeoElement<TPZGeoTetrahedra,TPZRefTetrahedra>::ClassId() const {
	return TPZFGEOELEMENTTETRAID;
}

template<>
int
TPZGeoElement<TPZGeoPyramid,TPZRefPyramid>::ClassId() const {
	return TPZFGEOELEMENTPYRAMID;
}

#endif

template class TPZGeoElement<TPZGeoPoint,TPZRefPoint>;
template class TPZGeoElement<TPZGeoLinear,TPZRefLinear>;
template class TPZGeoElement<TPZGeoTriangle,TPZRefTriangle>;
template class TPZGeoElement<TPZGeoQuad,TPZRefQuad>;
template class TPZGeoElement<TPZGeoCube,TPZRefCube>;
template class TPZGeoElement<TPZGeoPrism,TPZRefPrism>;
template class TPZGeoElement<TPZGeoTetrahedra,TPZRefTetrahedra>;
template class TPZGeoElement<TPZGeoPyramid,TPZRefPyramid>;
