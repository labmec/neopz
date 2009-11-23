//
// C++ Implementation: tpzquadrilateral
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpzquadrilateral.h"

#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"
#include "pzquad.h"
#include "tpzint1point.h"
#include "pzeltype.h"

#include "pzcreateapproxspace.h"

using namespace std;

namespace pztopology {

TPZCompEl *(*TPZQuadrilateral::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateQuadEl;	

TPZQuadrilateral::TPZQuadrilateral()
{
}


TPZQuadrilateral::~TPZQuadrilateral()
{
}

static int nsidenodes[9] = {
	1,1,1,1,2,2,2,2,4};

int TPZQuadrilateral::NSideNodes(int side)
{
	return nsidenodes[side];
}

static int nhighdimsides[9] = {3,3,3,3,1,1,1,1,0};


static int sidedimension[9] = {0,0,0,0,1,1,1,1,2};

static REAL sidetosidetransforms[9][3][4][3] = {
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}}
},
{
{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}}
},
{
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}}
},
{
{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}}
},
{
{{0,-1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
}
};


static int highsides[9][3] = {
{4,7,8},
{4,5,8},
{5,6,8},
{6,7,8},
{8},
{8},
{8},
{8},
{-999}
};

static REAL MidSideNode[9][3] = {
/*00*/{-1.,-1.},/*01*/{ 1.,-1.},/*02*/{1.,1.},
/*03*/{-1., 1.},/*04*/{ 0.,-1.},/*05*/{1.,0.},
/*06*/{ 0., 1.},/*07*/{-1., 0.},/*08*/{0.,0.} };

int TPZQuadrilateral::SideNodeLocId(int side, int node)
{
	if(side<4 && node==0) return side;
	if(side>=4 && side<8 && node <2) return (side+node)%4;
	if(side==8 && node <4) return node;
	PZError << "TPZQuadrilateral::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
	return -1;
}

void TPZQuadrilateral::LowerDimensionSides(int side,TPZStack<int> &smallsides)
{
     smallsides.Resize(0);
     int nsidecon = NSideConnects(side);
     int is;
     for(is=0; is<nsidecon-1; is++)
     smallsides.Push(SideConnectLocId(side,is));
}

void TPZQuadrilateral::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
{
     smallsides.Resize(0);
     int nsidecon = NSideConnects(side);
     for(int is = 0; is < nsidecon - 1; is++) {
     if (SideDimension(SideConnectLocId(side,is)) == DimTarget) smallsides.Push(SideConnectLocId(side,is));
  }
}

void TPZQuadrilateral::HigherDimensionSides(int side, TPZStack<int> &high)
{
	if(side <0 || side >= NSides) {
		PZError << "TPZQuadrilateral::HigherDimensionSides side "<< side << endl;
	}
	int is;
	for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
	
}

void TPZQuadrilateral::CenterPoint(int side, TPZVec<REAL> &center) {
  //center.Resize(Dimension);
  int i;
  for(i=0; i<Dimension; i++) {
    center[i] = MidSideNode[side][i];
  }
}

TPZIntPoints * TPZQuadrilateral::CreateSideIntegrationRule(int side, int order){
  if(side<0 || side>8) {
    PZError << "TPZQuadrilateral::CreateSideIntegrationRule wrong side " << side << endl;
    return 0;
  }
  if(side<4) return new TPZInt1Point();
  if(side<8) return new TPZInt1d(order);
  if(side==8) return new TPZIntQuad(order,order);
  return 0;
}


TPZTransform TPZQuadrilateral::TransformElementToSide(int side){

	if(side<0 || side>8){
  	PZError << "TPZShapeQuad::TransformElementToSide called with side error\n";
    return TPZTransform(0,0);
  }

  TPZTransform t(sidedimension[side],2);//t(dimto,2)
  t.Mult().Zero();	//TPZGeoElQ2d *gq;
  t.Sum().Zero();//int dimto = gq->SideDimension(side);

  switch(side){
    case 0:
    case 1:
    case 2:
    case 3:
      return t;
    case 4:
    	t.Mult()(0,0) = 1.0;//par. var.
      return t;
    case 5 :
      t.Mult()(0,1) = 1.0;
      return t;
    case 6:
    	t.Mult()(0,0) = -1.0;
      return t;
    case 7:
    	t.Mult()(0,1) = -1.0;
      return t;
    case 8:
    	t.Mult()(0,0) = 1.0;
      t.Mult()(1,1) = 1.0;
      return t;
  }
  return TPZTransform(0,0);
}

bool TPZQuadrilateral::IsInParametricDomain(TPZVec<REAL> &pt, REAL tol){
  const REAL qsi = pt[0];
  const REAL eta = pt[1];
  if( ( fabs(qsi) <= 1. + tol ) && ( fabs(eta) <= 1. + tol ) ){
    return true;
  }
  else{
    return false;
  }  
}///method

MElementType TPZQuadrilateral::Type()
{
  return EQuadrilateral;
}

MElementType TPZQuadrilateral::Type(int side)
{
  switch(side) {
    case 0:
    case 1:
    case 2:
    case 3:
      return EPoint;
    case 4:
    case 5:
    case 6:
    case 7:
      return EOned;
    case 8:
      return EQuadrilateral;
    default:
      return ENoType;
  }
}

int TPZQuadrilateral::NConnects() {
		return 9;
}

TPZTransform TPZQuadrilateral::SideToSideTransform(int sidefrom, int sideto)
{
	if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
		PZError << "TPZShapeQuad::HigherDimensionSides sidefrom "<< sidefrom << 
			' ' << sideto << endl;
		return TPZTransform(0);
	}
	if(sidefrom == sideto) {
		return TPZTransform(sidedimension[sidefrom]);
	}
	if(sidefrom == NSides-1) {
		return TransformElementToSide(sideto);
	}
	int nhigh = nhighdimsides[sidefrom];
	int is;
	for(is=0; is<nhigh; is++) {
		if(highsides[sidefrom][is] == sideto) {
			int dfr = sidedimension[sidefrom];
			int dto = sidedimension[sideto];
			TPZTransform trans(dto,dfr);
			int i,j;
			for(i=0; i<dto; i++) {
				for(j=0; j<dfr; j++) {
					trans.Mult()(i,j) = sidetosidetransforms[sidefrom][is][j][i];
				}
				trans.Sum()(i,0) = sidetosidetransforms[sidefrom][is][3][i];
			}
			return trans;
		}
	}
	PZError << "TPZShapeQuad::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
	return TPZTransform(0);
}


int TPZQuadrilateral::SideDimension(int side) {
	if(side<0 || side >= NSides) {
		PZError << "TPZShapeQuad::SideDimension side " << side << endl;
		return -1;
	}
	return sidedimension[side];
}

int TPZQuadrilateral::NSideConnects(int side) {
  if(side<0 || side>8) {
    PZError << "TPZShapeQuad::NSideConnects. Bad parameter i.\n";
    return 0;
  }
  if(side<4) return 1;
  if(side<8) return 3;
  return 9;//Cedric
}

/**It do not verify the values of the c*/
int TPZQuadrilateral::SideConnectLocId(int side,int c) {
  switch(side) {
  case 0:
  case 1:
  case 2:
  case 3:
    return side;
  case 4:
  case 5:
  case 6:
  case 7:
    if(!c) return side-4;
    if(c==1) return (side-3)%4;
    if(c==2) return side;
  case 8:
    return c;
  default:
    PZError << "TPZShapeQuad::SideConnectLocId, connect = " << c << endl;
    return -1;
  }
}


TPZTransform TPZQuadrilateral::TransformSideToElement(int side){

	if(side<0 || side>8){
  	PZError << "TPZShapeQuad::TransformSideToElement side out range\n";
    return TPZTransform(0,0);
  }
  TPZTransform t(2,sidedimension[side]);
  t.Mult().Zero();
  t.Sum().Zero();

  switch(side){
  	case 0:
    	t.Sum()(0,0) = -1.0;
      t.Sum()(1,0) = -1.0;
      return t;
    case 1:
    	t.Sum()(0,0) =  1.0;
      t.Sum()(1,0) = -1.0;
      return t;
    case 2:
    	t.Sum()(0,0) =  1.0;
      t.Sum()(1,0) =  1.0;
      return t;
    case 3:
    	t.Sum()(0,0) = -1.0;
      t.Sum()(1,0) =  1.0;
      return t;
    case 4:
      t.Mult()(0,0) =  1.0;
      t.Sum() (1,0) = -1.0;
      return t;
    case 5:
      t.Mult()(1,0) =  1.0;
      t.Sum() (0,0) =  1.0;
      return t;
    case 6:
      t.Mult()(0,0) = -1.0;
      t.Sum() (1,0) =  1.0;
      return t;
    case 7:
      t.Mult()(1,0) = -1.0;
      t.Sum() (0,0) = -1.0;
      return t;
    case 8:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      return t;
  }
	return TPZTransform(0,0);
}
}
