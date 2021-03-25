/**
 * @file
 * @brief Contains the instantiation of the TPZGeoElement classes from template.
 */

#include "TPZGeoElement.h"

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
#include "TPZStream.h"
//#include "pzstack.h"

using namespace pzgeom;
using namespace pzrefine;
using namespace pzshape;


/** Registration of the class in the TPZRestoreClass */

template class TPZRestoreClass< TPZGeoElement<TPZGeoPoint,TPZRefPoint>>;
template class TPZRestoreClass< TPZGeoElement<TPZGeoLinear,TPZRefLinear>>;
template class TPZRestoreClass< TPZGeoElement<TPZGeoQuad,TPZRefQuad>>;
template class TPZRestoreClass< TPZGeoElement<TPZGeoTriangle,TPZRefTriangle>>;
template class TPZRestoreClass< TPZGeoElement<TPZGeoCube,TPZRefCube>>;
template class TPZRestoreClass< TPZGeoElement<TPZGeoPrism,TPZRefPrism>>;
template class TPZRestoreClass< TPZGeoElement<TPZGeoTetrahedra,TPZRefTetrahedra>>;
template class TPZRestoreClass< TPZGeoElement<TPZGeoPyramid,TPZRefPyramid>>;

template class TPZGeoElement<TPZGeoPoint,TPZRefPoint>;
template class TPZGeoElement<TPZGeoLinear,TPZRefLinear>;
template class TPZGeoElement<TPZGeoTriangle,TPZRefTriangle>;
template class TPZGeoElement<TPZGeoQuad,TPZRefQuad>;
template class TPZGeoElement<TPZGeoCube,TPZRefCube>;
template class TPZGeoElement<TPZGeoPrism,TPZRefPrism>;
template class TPZGeoElement<TPZGeoTetrahedra,TPZRefTetrahedra>;
template class TPZGeoElement<TPZGeoPyramid,TPZRefPyramid>;
