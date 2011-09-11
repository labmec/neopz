/**
 * @file
 * @brief Implements a test to reduce dynamic memory allocation
 */
#ifndef PZGBMESHH
#define PZGBMESHH

#include "pzgmesh.h"
#include "TPZGeoElement.h"

namespace pzshape {
class TPZShapePoint;
class TPZShapeLinear;
class TPZShapeQuad;
class TPZShapeTriang;
class TPZShapeCube;
class TPZShapePrism;
class TPZShapeTetra;
class TPZShapePiram;
}
namespace pzgeom {
class TPZGeoPoint;
class TPZGeoLinear;
class TPZGeoQuad;
class TPZGeoTriangle;
class TPZGeoCube;
class TPZGeoPrism;
class TPZGeoTetrahedra;
class TPZGeoPyramid;
}
namespace pzrefine {
class TPZRefPoint;
class TPZRefLinear;
class TPZRefQuad;
class TPZRefTriangle;
class TPZRefCube;
class TPZRefPrism;
class TPZRefTetrahedra;
class TPZRefPyramid;
}

class TPZGeoEl;
//template<class TShape, class TGeo, class TRef>
//class TPZGeoElement<TShape,TGeo,TRef>;

/// Collection of typedefs for creation of geometric elements in groups
class GeoElTypes {

 public:

  typedef TPZGeoElement<pzgeom::TPZGeoPoint,pzrefine::TPZRefPoint> GPointType;
  typedef TPZGeoElement<pzgeom::TPZGeoLinear,pzrefine::TPZRefLinear> GLinearType;
  typedef TPZGeoElement<pzgeom::TPZGeoQuad,pzrefine::TPZRefQuad> GQuadType;
  typedef TPZGeoElement<pzgeom::TPZGeoTriangle,pzrefine::TPZRefTriangle> GTriangleType;
  typedef TPZGeoElement<pzgeom::TPZGeoCube,pzrefine::TPZRefCube> GHexahedronType;
  typedef TPZGeoElement<pzgeom::TPZGeoPrism,pzrefine::TPZRefPrism> GPrismType;
  typedef TPZGeoElement<pzgeom::TPZGeoTetrahedra,pzrefine::TPZRefTetrahedra> GTetrahedronType;
  typedef TPZGeoElement<pzgeom::TPZGeoPyramid,pzrefine::TPZRefPyramid> GPyramidType;

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
