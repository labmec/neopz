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


/** ClassId method for each instantiation followed by the registration of the class in the TPZRestoreClass */
template < >
int TPZGeoElRefPattern<TPZShapeCube,TPZGeoCube>::ClassId() const{
  return TPZGEOELREFPATCUBEID;
}
template class 
TPZRestoreClass< TPZGeoElRefPattern<TPZShapeCube,TPZGeoCube>, TPZGEOELREFPATCUBEID>;

template < >
int TPZGeoElRefPattern<TPZShapeLinear,TPZGeoLinear>::ClassId() const{
  return TPZGEOELREFPATLINEARID;
}
template class 
TPZRestoreClass< TPZGeoElRefPattern<TPZShapeLinear,TPZGeoLinear>, TPZGEOELREFPATLINEARID>;

template < >
int TPZGeoElRefPattern<TPZShapeQuad,TPZGeoQuad>::ClassId() const{
  return TPZGEOELREFPATQUADID;
}
template class 
TPZRestoreClass<TPZGeoElRefPattern<TPZShapeQuad,TPZGeoQuad>, TPZGEOELREFPATQUADID>;

template < >
int TPZGeoElRefPattern<TPZShapeTriang,TPZGeoTriangle>::ClassId() const{
  return TPZGEOELREFPATTRIANGLEID;
}
template class 
TPZRestoreClass<TPZGeoElRefPattern<TPZShapeTriang,TPZGeoTriangle>, TPZGEOELREFPATTRIANGLEID>;

template < >
int TPZGeoElRefPattern<TPZShapePrism,TPZGeoPrism>::ClassId() const{
  return TPZGEOELREFPATPRISMID;
}
template class 
TPZRestoreClass< TPZGeoElRefPattern<TPZShapePrism,TPZGeoPrism>, TPZGEOELREFPATPRISMID>;

template < >
int TPZGeoElRefPattern<TPZShapeTetra,TPZGeoTetrahedra>::ClassId() const{
  return TPZGEOELREFPATTETRAID;
}
template class 
TPZRestoreClass< TPZGeoElRefPattern<TPZShapeTetra,TPZGeoTetrahedra>, TPZGEOELREFPATTETRAID>;

template < >
int TPZGeoElRefPattern<TPZShapePiram,TPZGeoPyramid>::ClassId() const{
  return TPZGEOELREFPATPYRAMID;
}
template class 
TPZRestoreClass< TPZGeoElRefPattern<TPZShapePiram,TPZGeoPyramid>, TPZGEOELREFPATPYRAMID>;

template < >
int TPZGeoElRefPattern<TPZShapePoint,TPZGeoPoint>::ClassId() const{
  return TPZGEOELREFPATPOINTID;
}
template class 
TPZRestoreClass< TPZGeoElRefPattern<TPZShapePoint,TPZGeoPoint>, TPZGEOELREFPATPOINTID>;

template <class TShape, class TGeo>
void TPZGeoElRefPattern<TShape,TGeo>::Read(TPZStream &str, void *context){
  TPZGeoElRefLess<TShape,TGeo>::Read(str, context);
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

template <class TShape, class TGeo>
void TPZGeoElRefPattern<TShape,TGeo>::Write(TPZStream &str, int withclassid){
  TPZGeoElRefLess<TShape,TGeo>::Write(str, withclassid);
  int refpatternindex = -1;
  if(fRefPattern) refpatternindex = fRefPattern->Id();
  str.Write(&refpatternindex, 1);
  TPZSaveable::WriteObjects(str, this->fSubEl);
}

template <class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern(TPZGeoMesh &DestMesh, const TPZGeoElRefPattern<TShape,TGeo> &cp):TPZGeoElRefLess<TShape,TGeo>(DestMesh,cp),
  fRefPattern(cp.fRefPattern) {
  this->fSubEl = cp.fSubEl;
}

template <class TShape, class TGeo>
TPZGeoEl * TPZGeoElRefPattern<TShape,TGeo>::Clone(TPZGeoMesh &DestMesh) const{
  return new TPZGeoElRefPattern<TShape,TGeo>(DestMesh, *this);
}


