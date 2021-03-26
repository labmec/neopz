/**
 * @file
 * @brief Contains the implementation of the TPZGeoElRefPattern methods.
 */

#ifndef TPZGEOELREFPATTERN_H_H
#define TPZGEOELREFPATTERN_H_H

#include "pzlog.h"

#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"

template <class TGeo>
void TPZGeoElRefPattern<TGeo>::Read(TPZStream &str, void *context) {
    TPZGeoElRefLess<TGeo>::Read(str, context);
    str.Read(this->fSubEl);
    int refpatternindex;
    str.Read(&refpatternindex, 1);
    if (refpatternindex != -1) {
        const std::list< TPZAutoPointer<TPZRefPattern> > &RefPatternList = gRefDBase.RefPatternList(this->Type());
        std::list< TPZAutoPointer<TPZRefPattern> >::const_iterator it;

        for (it = RefPatternList.begin(); it != RefPatternList.end(); it++) {
            if ((*it)->Id() == refpatternindex) {
                break;
            }
        }

        if (it != RefPatternList.end()) {
            fRefPattern = *it;
        } else {
            DebugStop();
        }
    }
}

template <class TGeo>
void TPZGeoElRefPattern<TGeo>::Write(TPZStream &str, int withclassid) const {
    TPZGeoElRefLess<TGeo>::Write(str, withclassid);
    str.Write(this->fSubEl);
    int refpatternindex = -1;
    if (fRefPattern) refpatternindex = fRefPattern->Id();
    str.Write(&refpatternindex, 1);
}

template<class TGeo>
TPZGeoElRefPattern<TGeo>::TPZGeoElRefPattern():TPZRegisterClassId(&TPZGeoElRefPattern<TGeo>::ClassId),
TPZGeoElRefLess<TGeo>(), fSubEl(0), fRefPattern(0)
{
}

template<class TGeo>
TPZGeoElRefPattern<TGeo>::~TPZGeoElRefPattern() {
	//#warning "Acredito que os objetos deste tipo devem ser administrados pela malha"
	//  if (fRefPattern) delete fRefPattern;
}

template<class TGeo>
TPZGeoElRefPattern<TGeo>::TPZGeoElRefPattern(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh) :
TPZRegisterClassId(&TPZGeoElRefPattern<TGeo>::ClassId),TPZGeoElRefLess<TGeo>(nodeindices,matind,mesh) {
}
template<class TGeo>
TPZGeoElRefPattern<TGeo>::TPZGeoElRefPattern(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh, int64_t &index) :
TPZRegisterClassId(&TPZGeoElRefPattern<TGeo>::ClassId),TPZGeoElRefLess<TGeo>(nodeindices,matind,mesh,index)
{
}

template<class TGeo>
TPZGeoElRefPattern<TGeo>::TPZGeoElRefPattern(int64_t id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) :
TPZRegisterClassId(&TPZGeoElRefPattern<TGeo>::ClassId),TPZGeoElRefLess<TGeo>(id,nodeindexes,matind,mesh) {
}


template <class TGeo>
TPZGeoElRefPattern<TGeo>::TPZGeoElRefPattern(TPZGeoMesh &DestMesh, const TPZGeoElRefPattern<TGeo> &cp):
TPZRegisterClassId(&TPZGeoElRefPattern<TGeo>::ClassId),TPZGeoElRefLess<TGeo>(DestMesh,cp),
fRefPattern(cp.fRefPattern) {

	this->fSubEl = cp.fSubEl;
}

template <class TGeo>
TPZGeoEl * TPZGeoElRefPattern<TGeo>::Clone(TPZGeoMesh &DestMesh) const{
    if(&DestMesh == this->Mesh())
    {
        return new TPZGeoElRefPattern<TGeo>(*this);
    }
    else
    {
        return new TPZGeoElRefPattern<TGeo>(DestMesh, *this);
    }
}


template <class TGeo>
TPZGeoElRefPattern<TGeo>::TPZGeoElRefPattern(TPZGeoMesh &DestMesh,
											 const TPZGeoElRefPattern<TGeo> &cp,
											 std::map<int64_t,int64_t> &gl2lcNdMap,
											 std::map<int64_t,int64_t> &gl2lcElMap):
TPZRegisterClassId(&TPZGeoElRefPattern<TGeo>::ClassId),TPZGeoElRefLess<TGeo>(DestMesh,cp,gl2lcNdMap,gl2lcElMap),
fRefPattern ( cp.fRefPattern )
{
	int i;
    this->fSubEl.Resize(cp.fSubEl.size(),0);
	for (i=0;i<cp.fSubEl.NElements();i++)
	{
		if (cp.fSubEl[i] == -1)
		{
			this->fSubEl[i] = -1;
			continue;
		}
		if (gl2lcElMap.find(cp.fSubEl[i]) == gl2lcElMap.end())
		{
			std::stringstream sout;
			sout << "ERROR in - " << __PRETTY_FUNCTION__
			<< " subelement " << i << " index = " << cp.fSubEl[i] << " is not in the map.";
#ifdef PZ_LOG
            TPZLogger loggerrefpattern("pz.mesh.tpzgeoelrefpattern");
			LOGPZ_ERROR (loggerrefpattern,sout.str().c_str());
#endif
			DebugStop();
		}
		this->fSubEl[i] = gl2lcElMap[cp.fSubEl[i]];
	}
}

template <class TGeo>
TPZGeoEl * TPZGeoElRefPattern<TGeo>::ClonePatchEl(TPZGeoMesh &DestMesh,
												  std::map<int64_t,int64_t> &gl2lcNdMap,
												  std::map<int64_t,int64_t> &gl2lcElMap) const{
	return new TPZGeoElRefPattern<TGeo>(DestMesh, *this, gl2lcNdMap, gl2lcElMap);
}

#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgeoprism.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "pzgeopoint.h"


#define INSERTCLASS(TCL) \
template class \
TPZRestoreClass<TPZGeoElRefPattern<TCL > >; \
template class TPZGeoElRefPattern<TCL>;

INSERTCLASS(pzgeom::TPZGeoPoint)
INSERTCLASS(pzgeom::TPZGeoLinear)
INSERTCLASS(pzgeom::TPZGeoTriangle)
INSERTCLASS(pzgeom::TPZGeoQuad)
INSERTCLASS(pzgeom::TPZGeoCube)
INSERTCLASS(pzgeom::TPZGeoPrism)
INSERTCLASS(pzgeom::TPZGeoTetrahedra)
INSERTCLASS(pzgeom::TPZGeoPyramid)

#endif
