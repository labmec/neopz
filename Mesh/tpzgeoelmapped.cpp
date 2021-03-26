/**
 * @file
 * @brief Contains the implementation of the TPZGeoElMapped methods. 
 */

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


/** Creates a geometric element according to the type of the father element */
template<class TBase>
TPZGeoEl *TPZGeoElMapped<TBase>::CreateGeoElement(MElementType type,
												  TPZVec<int64_t>& nodeindexes,
												  int matid,
												  int64_t& index)
{
	TPZGeoMesh &mesh = *(this->Mesh());
	return mesh.CreateGeoElementMapped(type,nodeindexes,matid,index);
}

template<class TBase>
TPZGeoEl * TPZGeoElMapped<TBase>::Clone(TPZGeoMesh &DestMesh) const
{
    if(&DestMesh == this->Mesh())
    {
        return new TPZGeoElMapped<TBase>(*this);
    }
    else
    {
        return new TPZGeoElMapped<TBase>(DestMesh,*this);
    }
}

/** @} */

template<class TBase>
TPZGeoEl * TPZGeoElMapped<TBase>::ClonePatchEl(TPZGeoMesh &DestMesh,
                                std::map<int64_t,int64_t> &gl2lcNdIdx,
                                std::map<int64_t,int64_t> &gl2lcElIdx) const
{
    return new TPZGeoElMapped<TBase>(DestMesh,*this,gl2lcNdIdx,gl2lcElIdx);
}

template <class TBase>
TPZGeoElMapped<TBase>::TPZGeoElMapped(TPZGeoMesh &destmesh, const TPZGeoElMapped<TBase> &copy) : TPZRegisterClassId(&TPZGeoElMapped::ClassId),
TBase(destmesh,copy), fCornerCo(copy.fCornerCo)
{
    
}

template <class TBase>
TPZGeoElMapped<TBase>::TPZGeoElMapped(TPZGeoMesh &destmesh, const TPZGeoElMapped<TBase> &copy, std::map<int64_t,int64_t> &gl2lcNdIdx,
                                      std::map<int64_t,int64_t> &gl2lcElIdx) : 
TPZRegisterClassId(&TPZGeoElMapped::ClassId),TBase(destmesh,copy,gl2lcNdIdx,gl2lcElIdx),
    fCornerCo(copy.fCornerCo)
{
    
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




using namespace pzgeom;

/// Macro to define templates to TPZGeoElMapped for all the geometric element types
#define INSERTCLASS2(TCL) \
template class \
TPZRestoreClass< TPZGeoElMapped<TPZGeoElRefPattern<TCL > >>; \
template class TPZGeoElMapped< TPZGeoElRefPattern<TCL> >;

INSERTCLASS2(TPZGeoPoint)
INSERTCLASS2(TPZGeoLinear)
INSERTCLASS2(TPZGeoTriangle)
INSERTCLASS2(TPZGeoQuad)
INSERTCLASS2(TPZGeoCube)
INSERTCLASS2(TPZGeoPrism)
INSERTCLASS2(TPZGeoTetrahedra)
INSERTCLASS2(TPZGeoPyramid)
