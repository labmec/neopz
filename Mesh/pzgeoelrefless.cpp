/**
 * @file
 * @brief Creates TPZGeoElRefLess classes for all topological master elements.
 */

#include "pzgeoelrefless.h"
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
#include "tpzellipse3d.h"
#include "tpzarc3d.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pzvec.h"
#include "pzmanvector.h"

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
#include "tpzquadraticquad.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticline.h"
#include "tpzquadraticcube.h"
#include "tpzquadraticprism.h"
#include "tpzquadraticpyramid.h"
#include "tpzquadratictetra.h"


template class TPZGeoElRefLess<TPZGeoCube>;
template class TPZGeoElRefLess<TPZGeoLinear>;
template class TPZGeoElRefLess<TPZGeoQuad>;
template class TPZGeoElRefLess<TPZGeoTriangle>;
template class TPZGeoElRefLess<TPZGeoPrism>;
template class TPZGeoElRefLess<TPZGeoTetrahedra>;
template class TPZGeoElRefLess<TPZGeoPyramid>;
template class TPZGeoElRefLess<TPZGeoPoint>;
template class TPZGeoElRefLess<TPZQuadraticQuad>;
template class TPZGeoElRefLess<TPZQuadraticTrig>;
template class TPZGeoElRefLess<TPZQuadraticLine>;
template class TPZGeoElRefLess<TPZQuadraticCube>;
template class TPZGeoElRefLess<TPZQuadraticPrism>;
template class TPZGeoElRefLess<TPZQuadraticPyramid>;
template class TPZGeoElRefLess<TPZQuadraticTetra>;

//static int main_refless()
//{
	
//	TPZGeoEl * teste = new TPZGeoElRefLess<TPZGeoTriangle>;
//	if(teste)
//		return 0;
//	return 1;
//}
