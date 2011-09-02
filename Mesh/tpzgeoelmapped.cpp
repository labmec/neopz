/**
 * @file
 * @brief Contains the implementation of the TPZGeoElMapped methods. 
 */
//
// C++ Implementation: tpzgeoelmapped
//
// Description: 
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
#include "tpzgeoelmapped.h"

#include "tpzgeoelrefpattern.h"
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
#include "pzgeopoint.h"
#include "pzrefpoint.h"
#include "pzshapepoint.h"

using namespace pzgeom;
using namespace pzrefine;
using namespace pzshape;


/**
 * Creates a geometric element according to the type of the father element
 */
template<class TBase>
TPZGeoEl *TPZGeoElMapped<TBase>::CreateGeoElement(MElementType type,
												  TPZVec<int>& nodeindexes,
												  int matid,
												  int& index)
{
	TPZGeoMesh &mesh = *(this->Mesh());
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}



#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgeoprism.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "pzgeopoint.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "TPZGeoElement.h"

#include "pzvec.h"
#include "pzmanvector.h"
#include "tpzgeoelrefpattern.h"


TPZGeoEl *CreateGeoElementMapped(TPZGeoMesh &mesh,
								 MElementType type,
								 TPZVec<int>& nodeindexes,
								 int matid,
								 int& index)
{
	{
		if(!&mesh) return 0;
		switch( type ){
			case 0://point
			{
				TPZGeoEl * gel =
				new TPZGeoElMapped<TPZGeoElRefPattern<TPZGeoPoint> > (nodeindexes, matid, mesh, index);
				return gel;
			}
			case 1://line
			{
				TPZGeoEl *gel =
				new TPZGeoElMapped<TPZGeoElRefPattern < TPZGeoLinear > >
				(nodeindexes, matid, mesh, index);
				return gel;
			}
			case 2://triangle
			{
				TPZGeoEl *gel =
				new TPZGeoElMapped<TPZGeoElRefPattern < TPZGeoTriangle > >
				(nodeindexes, matid, mesh, index);
				return gel;
			}
			case 3://quadrilatera
			{
				TPZGeoEl* gel =
				new TPZGeoElMapped<TPZGeoElRefPattern < TPZGeoQuad > >
				(nodeindexes, matid, mesh, index);
				return gel;
			}
			case 4://tetraedra
			{
				TPZGeoEl*gel =
				new TPZGeoElMapped<TPZGeoElRefPattern < TPZGeoTetrahedra > >
				(nodeindexes, matid, mesh, index);
				return gel;
			}
			case 5://pyramid
			{
				TPZGeoEl *gel =
				new TPZGeoElMapped<TPZGeoElRefPattern < TPZGeoPyramid > >
				(nodeindexes, matid, mesh, index);
				return gel;
			}
			case 6://prism
			{
				TPZGeoEl*gel =
				new TPZGeoElMapped<TPZGeoElRefPattern < TPZGeoPrism > >
				(nodeindexes, matid, mesh, index);
				return gel;
			}
			case 7://cube
			{
				TPZGeoEl*gel =
				new TPZGeoElMapped<TPZGeoElRefPattern < TPZGeoCube > >
				(nodeindexes, matid, mesh, index);
				return gel;
			}
			default:
			{
				PZError << "TPZGeoMesh::CreateGeoElement type element not exists:"
				<< " type = " << type << std::endl;
				return NULL;
			}
		}
	}
}

using namespace pzgeom;

/// Macro to define templates to TPZGeoElMapped for all the geometric element types
#define INSERTCLASS(TCL,CLID) \
template<> \
int TPZGeoElMapped<TPZGeoElRefPattern< TCL > >::ClassId() const \
{ \
return CLID; \
} \
template class \
TPZRestoreClass< TPZGeoElMapped<TPZGeoElRefPattern<TCL > >, CLID>; \
template class TPZGeoElMapped< TPZGeoElRefPattern<TCL> >;

INSERTCLASS(TPZGeoPoint,TPZGEOELREFPATMAPPEDPOINTID)
INSERTCLASS(TPZGeoLinear,TPZGEOELREFPATMAPPEDLINEID)
INSERTCLASS(TPZGeoTriangle,TPZGEOELREFPATMAPPEDTRIANGLEID)
INSERTCLASS(TPZGeoQuad,TPZGEOELREFPATMAPPEDQUADRILATERALID)
INSERTCLASS(TPZGeoCube,TPZGEOELREFPATMAPPEDCUBEID)
INSERTCLASS(TPZGeoPrism,TPZGEOELREFPATMAPPEDPRISMID)
INSERTCLASS(TPZGeoTetrahedra,TPZGEOELREFPATMAPPEDTETRAHEDRAID)
INSERTCLASS(TPZGeoPyramid,TPZGEOELREFPATMAPPEDPYRAMIDID)

//template class TPZGeoElMapped< TPZGeoElRefPattern<TPZGeoCube> >;
// template class TPZGeoElMapped< TPZGeoElRefPattern<TPZGeoPoint> >;
// template class TPZGeoElMapped< TPZGeoElRefPattern<TPZGeoQuad> >;
// template class TPZGeoElMapped< TPZGeoElRefPattern<TPZGeoTriangle> >;
// template class TPZGeoElMapped< TPZGeoElRefPattern<TPZGeoPrism> >;
// template class TPZGeoElMapped< TPZGeoElRefPattern<TPZGeoTetrahedra> >;
// template class TPZGeoElMapped< TPZGeoElRefPattern<TPZGeoPyramid> >;
