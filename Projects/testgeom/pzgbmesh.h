// -*- c++ -*-
#ifndef PZGBMESHH
#define PZGBMESHH

#include "pzgmesh.h"
#include "TPZGeoElement.h"

class TPZShapePoint;
class TPZGeoPoint;
class TPZRefPoint;
class TPZShapeLinear;
class TPZGeoLinear;
class TPZRefLinear;
class TPZShapeQuad;
class TPZGeoQuad;
class TPZRefQuad;
class TPZShapeTriang;
class TPZGeoTriangle;
class TPZRefTriangle;
class TPZShapeCube;
class TPZGeoCube;
class TPZRefCube;
class TPZShapePrism;
class TPZGeoPrism;
class TPZRefPrism;
class TPZShapeTetra;
class TPZGeoTetrahedra;
class TPZRefTetrahedra;
class TPZShapePiram;
class TPZGeoPyramid;
class TPZRefPyramid;

class TPZGeoEl;
//template<class TShape, class TGeo, class TRef>
//class TPZGeoElement<TShape,TGeo,TRef>;

/// Collection of typedefs for creation of geometric elements in groups
class GeoElTypes {

 public:

  typedef TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint> GPointType;
  typedef TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear> GLinearType;
  typedef TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad> GQuadType;
  typedef TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle> GTriangleType;
  typedef TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube> GHexahedronType;
  typedef TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism> GPrismType;
  typedef TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra> GTetrahedronType;
  typedef TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid> GPyramidType;

};

/// A geometric mesh which generates elements in blocks
/**
The idea behind this development is to reduce the number of dynamic memory allocations
*/
template<class Types>
class TPZGeoBMesh : public TPZGeoMesh {

  TPZAdmChunkVector<typename Types::GPointType> fPoint;
  TPZAdmChunkVector<typename Types::GLinearType> fLinear;
  TPZAdmChunkVector<typename Types::GQuadType> fQuad;
  TPZAdmChunkVector<typename Types::GTriangleType> fTriangle;
  TPZAdmChunkVector<typename Types::GHexahedronType> fHexahedron;
  TPZAdmChunkVector<typename Types::GPrismType> fPrism;
  TPZAdmChunkVector<typename Types::GTetrahedronType> fTetrahedron;
  TPZAdmChunkVector<typename Types::GPyramidType> fPyramid;

 public:

  TPZGeoBMesh() : TPZGeoMesh() {
  }

  virtual ~TPZGeoBMesh() {}

  virtual  TPZGeoEl *CreateGeoElement(MElementType type,TPZVec<int> &cornerindexes,int matid,int &index);

  virtual void DeleteElement(int gelindex);

  static int main();

};



#endif
