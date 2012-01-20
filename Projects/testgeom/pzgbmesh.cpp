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

#include "TPZRefPatternDataBase.h"
#include "pzgengrid.h"

using namespace std;

// Program to reproduce issue 1 into CodeGoogle NeoPZ

void UniformRefine(TPZAutoPointer<TPZGeoMesh> gmesh, int nDiv);

int main() {
    
#ifdef LOG4CXX
	InitializePZLOG();
#endif

    // First rectangular mesh
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;

	TPZManVector<int> nx(2,2);   // subdivisions in X and in Y
	TPZManVector<REAL> x0(3,0.), x1(3,1.);  // Corners of the rectangular mesh
	x1[2]=0.;

	TPZGenGrid gen(nx,x0,x1);    // mesh generator 
	gen.SetElementType(0);       // type = 0 means rectangular elements
	gen.Read(gmesh);            // generating mesh in gmesh

	ofstream saida("malhateste.txt");
	char namemesh[260];
	strncpy(namemesh,"Malha inicial",strlen("Malha inicial")+1);
	gen.Print(namemesh,saida);
//	gmesh->Print(saida);
	
	// Second rectangular domain - subdividions and corners of the second rectangular mesh
    TPZAutoPointer<TPZGeoMesh> gmesh2 = new TPZGeoMesh;
	nx[0] = nx[1] = 4;
	x0[1] = 1.;
	x1[0] = 2.;
	x1[1] = 3.;

	TPZGenGrid gen2(nx,x0,x1);   // second mesh generator
	gen2.SetElementType(0);

	// generating gmesh2 on data of the gen2 and merge gmesh into the gmesh2
	gen2.ReadAndMergeGeoMesh(gmesh2,gmesh);
	x0[1] = 0.;
	x1[1] = 1.;
	// setting bc condition -1 [no flux - is wall] from (0.,0.) until (2.,1.)
	gen2.SetBC(gmesh2,x0,x1,-1);
	x0[0] = 2.;
	x0[1] = 3.;
	x1[0] = 0.;
	x1[1] = 3.;
	// setting bc condition -1 from (2.,3.) until (0.,3.)
	gen2.SetBC(gmesh2,x0,x1,-1);
	x1[0] = 2.;
	x1[1] = 1.;
	// setting bc condition -2 from (2.,1.) until (2.,3.)
	gen2.SetBC(gmesh2,x1,x0,-2);
	x0[0] = 0.;
	x1[0] = x1[1] = 0.;
	// setting bc condition -3 from (0.,0.) until (0.,3.)
	gen2.SetBC(gmesh2, x0, x1, -3);
//	gmesh2->Print(saida);
	
	// Uniform refinement of the geometrical mesh, two level
	int nDiv = 2;
	UniformRefine(gmesh2, nDiv);
	gmesh2->Print(saida);

	saida.close();
	return 0;
}

void UniformRefine(TPZAutoPointer<TPZGeoMesh> gmesh, int nDiv)
{
    for(int D = 0; D < nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {    
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
}

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