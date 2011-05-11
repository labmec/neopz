#include "pzgbmesh.h"
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
#include "pzshapepoint.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"

using namespace std;

/*
template<class T>
TPZGeoEl *TPZGeoBMesh<T>::CreateGeoElement(MElementType type,TPZVec<int> &nodeindexes,int matid,int &index) {

  int newelindex;
   switch( type )
   {
      case EPoint://point
	newelindex = fPoint.AllocateNewElement();
	fPoint[newelindex].Initialize(nodeindexes, matid, *this, index);
	return &fPoint[newelindex];

      case EOned://line
	newelindex = fLinear.AllocateNewElement();
	fLinear[newelindex].Initialize(nodeindexes, matid, *this, index);
	return &fLinear[newelindex];

      case ETriangle://triangle
	newelindex = fTriangle.AllocateNewElement();
	fTriangle[newelindex].Initialize(nodeindexes, matid, *this, index);
	return &fTriangle[newelindex];

      case EQuadrilateral://quadrilatera
	newelindex = fQuad.AllocateNewElement();
	fQuad[newelindex].Initialize(nodeindexes, matid, *this, index);
	return &fQuad[newelindex];

      case ETetraedro://tetraedra
	newelindex = fTetrahedron.AllocateNewElement();
	fTetrahedron[newelindex].Initialize(nodeindexes, matid, *this, index);
	return &fTetrahedron[newelindex];

      case EPiramide:
	newelindex = fPyramid.AllocateNewElement();
	fPyramid[newelindex].Initialize(nodeindexes, matid, *this, index);
	return &fPyramid[newelindex];

      case EPrisma:
	newelindex = fPrism.AllocateNewElement();
	fPrism[newelindex].Initialize(nodeindexes, matid, *this, index);
	return &fPrism[newelindex];

      case ECube:
	newelindex = fHexahedron.AllocateNewElement();
	fHexahedron[newelindex].Initialize(nodeindexes, matid, *this, index);
	return &fHexahedron[newelindex];

      default:
	 PZError << "TPZGeoBMesh::CreateGeoElement type element not exists:"
		 << " type = " << type << endl;
	 return NULL;
   }

   return NULL;
}

template<class T>
void TPZGeoBMesh<T>::DeleteElement(int elindex) {

  TPZGeoEl *gel = ElementVec()[elindex];
  if(!gel) return;
  ElementVec().SetFree(elindex);
  ElementVec()[elindex] = 0;
  int blindex = -1;
  switch(gel->Type()) {
  case EPoint: {//point
    typename T::GPointType *gstr = dynamic_cast<typename T::GPointType *>(gel);
    if(gstr) blindex = fPoint.FindObject(gstr);
    if(blindex != -1) {
      fPoint.SetFree(blindex);
    } else {
      delete gel;
    }
  }
    break;
    
  case EOned: {//line
    typename T::GLinearType *gstr = dynamic_cast<typename T::GLinearType *>(gel);
    if(gstr) blindex = fLinear.FindObject(gstr);
    if(blindex != -1) {
      fLinear.SetFree(blindex);
    } else {
      delete gel;
    }
  }
    break;

  case ETriangle: {//triangle
    typename T::GTriangleType *gstr = dynamic_cast<typename T::GTriangleType *>(gel);
    if(gstr) blindex = fTriangle.FindObject(gstr);
    if(blindex != -1) {
      fTriangle.SetFree(blindex);
    } else {
      delete gel;
    }
  }
    break;

  case EQuadrilateral: {//quadrilatera
    typename T::GQuadType *gstr = dynamic_cast<typename T::GQuadType *>(gel);
    if(gstr) blindex = fQuad.FindObject(gstr);
    if(blindex != -1) {
      fQuad.SetFree(blindex);
    } else {
      delete gel;
    }
  }
    break;

  case ETetraedro: {//tetraedra
    typename T::GTetrahedronType *gstr = dynamic_cast<typename T::GTetrahedronType *>(gel);
    if(gstr) blindex = fTetrahedron.FindObject(gstr);
    if(blindex != -1) {
      fTetrahedron.SetFree(blindex);
    } else {
      delete gel;
    }
  }
    break;

  case EPiramide: {
    typename T::GPyramidType *gstr = dynamic_cast<typename T::GPyramidType *>(gel);
    if(gstr) blindex = fPyramid.FindObject(gstr);
    if(blindex != -1) {
      fPyramid.SetFree(blindex);
    } else {
      delete gel;
    }
  }
    break;

  case EPrisma: {
    typename T::GPrismType *gstr = dynamic_cast<typename T::GPrismType *>(gel);
    if(gstr) blindex = fPrism.FindObject(gstr);
    if(blindex != -1) {
      fPrism.SetFree(blindex);
    } else {
      delete gel;
    }
  }
    break;

  case ECube: {
    typename T::GHexahedronType *gstr = dynamic_cast<typename T::GHexahedronType *>(gel);
    if(gstr) blindex = fHexahedron.FindObject(gstr);
    if(blindex != -1) {
      fHexahedron.SetFree(blindex);
    } else {
      delete gel;
    }
  }
    break;
    
  default:
    delete gel;
    PZError << "TPZGeoBMesh::CreateGeoElement type element not exists:"
	    << " type = " << gel->Type() << endl;
  }
}

template<class T>
int TPZGeoBMesh<T>::main() {
  return 0;
}

template class TPZGeoBMesh<GeoElTypes>;
*/