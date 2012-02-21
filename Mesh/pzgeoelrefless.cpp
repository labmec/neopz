/**
 * @file
 * @brief Creates TPZGeoElRefLess classes for all topological master elements.
 */
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

template class TPZGeoElRefLess<TPZGeoCube>;
template class TPZGeoElRefLess<TPZGeoLinear>;
template class TPZGeoElRefLess<TPZGeoQuad>;
template class TPZGeoElRefLess<TPZGeoTriangle>;
template class TPZGeoElRefLess<TPZGeoPrism>;
template class TPZGeoElRefLess<TPZGeoTetrahedra>;
template class TPZGeoElRefLess<TPZGeoPyramid>;
template class TPZGeoElRefLess<TPZGeoPoint>;

static int main_refless()
{
	
	TPZGeoEl * teste = new TPZGeoElRefLess<TPZGeoTriangle>;
	if(teste)
		return 0;
	return 1;
}

