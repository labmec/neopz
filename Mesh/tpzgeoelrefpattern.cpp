/***************************************************************************
                          tpzgeoelrefpattern.cc  -  description
                             -------------------
    begin                : Tue Dec 23 2003
    copyright            : (C) 2003 by LabMeC - DES - FEC - UNICAMP (Edimar Cesar Rylo) & EMBRAER
    email                : cesar@labmec.fec.unicamp.br
 ***************************************************************************/

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
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "TPZGeoElement.h"
#include "TPZRefPattern.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzlog.h"
using namespace pzgeom;
using namespace pzshape;

// class TPZGeoElRefPattern<TPZShapeCube,TPZGeoCube>;
// class TPZGeoElRefPattern<TPZShapeLinear,TPZGeoLinear>;
// class TPZGeoElRefPattern<TPZShapeQuad,TPZGeoQuad>;
// class TPZGeoElRefPattern<TPZShapeTriang,TPZGeoTriangle>;
// class TPZGeoElRefPattern<TPZShapePrism,TPZGeoPrism>;
// class TPZGeoElRefPattern<TPZShapeTetra,TPZGeoTetrahedra>;
// class TPZGeoElRefPattern<TPZShapePiram,TPZGeoPyramid>;
// class TPZGeoElRefPattern<TPZShapePoint,TPZGeoPoint>;
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzgeoelrefpattern"));
#endif

/** ClassId method for each instantiation followed by the registration of the class in the TPZRestoreClass */
template < >
int TPZGeoElRefPattern<TPZGeoCube>::ClassId() const{
  return TPZGEOELREFPATCUBEID;
}
template class
TPZRestoreClass< TPZGeoElRefPattern<TPZGeoCube>, TPZGEOELREFPATCUBEID>;

template < >
int TPZGeoElRefPattern<TPZGeoLinear>::ClassId() const{
  return TPZGEOELREFPATLINEARID;
}
template class
TPZRestoreClass< TPZGeoElRefPattern<TPZGeoLinear>, TPZGEOELREFPATLINEARID>;

template < >
int TPZGeoElRefPattern<TPZGeoQuad>::ClassId() const{
  return TPZGEOELREFPATQUADID;
}
template class
TPZRestoreClass<TPZGeoElRefPattern<TPZGeoQuad>, TPZGEOELREFPATQUADID>;

template < >
int TPZGeoElRefPattern<TPZGeoTriangle>::ClassId() const{
  return TPZGEOELREFPATTRIANGLEID;
}
template class
TPZRestoreClass<TPZGeoElRefPattern<TPZGeoTriangle>, TPZGEOELREFPATTRIANGLEID>;

template < >
int TPZGeoElRefPattern<TPZGeoPrism>::ClassId() const{
  return TPZGEOELREFPATPRISMID;
}
template class
TPZRestoreClass< TPZGeoElRefPattern<TPZGeoPrism>, TPZGEOELREFPATPRISMID>;

template < >
int TPZGeoElRefPattern<TPZGeoTetrahedra>::ClassId() const{
  return TPZGEOELREFPATTETRAID;
}
template class
TPZRestoreClass< TPZGeoElRefPattern<TPZGeoTetrahedra>, TPZGEOELREFPATTETRAID>;

template < >
int TPZGeoElRefPattern<TPZGeoPyramid>::ClassId() const{
  return TPZGEOELREFPATPYRAMID;
}
template class
TPZRestoreClass< TPZGeoElRefPattern<TPZGeoPyramid>, TPZGEOELREFPATPYRAMID>;

template < >
int TPZGeoElRefPattern<TPZGeoPoint>::ClassId() const{
  return TPZGEOELREFPATPOINTID;
}
template class
TPZRestoreClass< TPZGeoElRefPattern<TPZGeoPoint>, TPZGEOELREFPATPOINTID>;

template <class TGeo>
void TPZGeoElRefPattern<TGeo>::Read(TPZStream &str, void *context){
  TPZGeoElRefLess<TGeo>::Read(str, context);
  TPZGeoMesh *gmesh = (TPZGeoMesh *) context;
  int refpatternindex;
  str.Read(&refpatternindex, 1);
  if(refpatternindex != -1)
  {
    const std::map<int, TPZAutoPointer<TPZRefPattern> > &RefPatternList = gmesh->RefPatternList(this->Type());
    std::map<int, TPZAutoPointer<TPZRefPattern> >::const_iterator it;
    it = RefPatternList.find(refpatternindex);
    if(it != RefPatternList.end()) fRefPattern = it->second;
  }
  TPZSaveable::ReadObjects(str, this->fSubEl);
}

template <class TGeo>
void TPZGeoElRefPattern<TGeo>::Write(TPZStream &str, int withclassid){
  TPZGeoElRefLess<TGeo>::Write(str, withclassid);
  int refpatternindex = -1;
  if(fRefPattern) refpatternindex = fRefPattern->Id();
  str.Write(&refpatternindex, 1);
  TPZSaveable::WriteObjects(str, this->fSubEl);
}

template <class TGeo>
TPZGeoElRefPattern<TGeo>::TPZGeoElRefPattern(TPZGeoMesh &DestMesh, const TPZGeoElRefPattern<TGeo> &cp):TPZGeoElRefLess<TGeo>(DestMesh,cp),
  fRefPattern(cp.fRefPattern) {
  this->fSubEl = cp.fSubEl;
}

template <class TGeo>
TPZGeoEl * TPZGeoElRefPattern<TGeo>::Clone(TPZGeoMesh &DestMesh) const{
  return new TPZGeoElRefPattern<TGeo>(DestMesh, *this);
}


template <class TGeo>
TPZGeoElRefPattern<TGeo>::TPZGeoElRefPattern(TPZGeoMesh &DestMesh,
                                                    const TPZGeoElRefPattern<TGeo> &cp,
                                                    std::map<int,int> &gl2lcNdMap,
                                                    std::map<int,int> &gl2lcElMap):
                                                    TPZGeoElRefLess<TGeo>(DestMesh,cp,gl2lcNdMap,gl2lcElMap),
                                                    fRefPattern ( cp.fRefPattern )
{
  int i;
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
      LOGPZ_ERROR (logger,sout.str().c_str());
      exit(-1);
    }
    this->fSubEl[i] = gl2lcElMap[cp.fSubEl[i]];
  }
}


template <class TGeo>
TPZGeoEl * TPZGeoElRefPattern<TGeo>::ClonePatchEl(TPZGeoMesh &DestMesh,
                                                        std::map<int,int> &gl2lcNdMap,
                                                        std::map<int,int> &gl2lcElMap) const{
  return new TPZGeoElRefPattern<TGeo>(DestMesh, *this, gl2lcNdMap, gl2lcElMap);
}

