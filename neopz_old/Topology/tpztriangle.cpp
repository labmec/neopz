//
// C++ Implementation: tpztriangle
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpztriangle.h"
#include "pzquad.h"
#include "tpzint1point.h"

#include "pzcreateapproxspace.h"

using namespace std;

namespace pztopology {

TPZCompEl *(*TPZTriangle::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateTriangleEl;
	

static int sidedimension[7] = {0,0,0,1,1,1,2};

static int nhighdimsides[7] = {3,3,3,1,1,1,0};

static int highsides[7][3] = {
{3,5,6},
{3,4,6},
{4,5,6},
{6},
{6},
{6},
{-999}
};

static REAL sidetosidetransforms[7][3][4][3] = {
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}}
},
{
{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}}
},
{
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}}
},
{
{{0,-0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
}
};

static REAL MidSideNode[7][3] = {
/*00*/{.0,0.},/*01*/{1.0,.0},/*02*/{0.,1.0},
/*03*/{.5,0.},/*04*/{0.5,.5},/*05*/{0.,0.5},
/*06*/{ 1./3.,1./3.} };



static int nsidenodes[7] = {1,1,1,2,2,2,3};

int TPZTriangle::NSideNodes(int side)
{
	return nsidenodes[side];
}

bool TPZTriangle::IsInParametricDomain(TPZVec<REAL> &pt, REAL tol){
  const REAL qsi = pt[0];
  const REAL eta = pt[1];
  if( ( qsi <= 1. + tol ) && ( qsi >= 0. - tol ) &&
      ( eta <= 1. + tol ) && ( eta >= 0. - tol ) &&
      ( eta <= 1. - qsi + tol ) ){
    return true;
  }
  else{
    return false;
  }  
}///method

TPZTransform TPZTriangle::SideToSideTransform(int sidefrom, int sideto)
{
	if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
		PZError << "TPZTriangle::HigherDimensionSides sidefrom "<< sidefrom << 
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
	PZError << "TPZTriangle::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
	return TPZTransform(0);
}

int TPZTriangle::SideNodeLocId(int side, int node)
{
	if(side<3 && node == 0) return side;
	if(side>=3 && side<6 && node <2) return (side-3+node) %3;
	if(side==6 && node <3) return node;
	PZError << "TPZTriangle::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
	return -1;
}

void TPZTriangle::LowerDimensionSides(int side,TPZStack<int> &smallsides)
{
     smallsides.Resize(0);
     int nsidecon = NSideConnects(side);
     int is;
     for(is=0; is<nsidecon-1; is++)
     smallsides.Push(SideConnectLocId(side,is));
}

void TPZTriangle::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
{
  smallsides.Resize(0);
  int nsidecon = NSideConnects(side);
  for(int is = 0; is < nsidecon - 1; is++) {
    if (SideDimension(SideConnectLocId(side,is)) == DimTarget) smallsides.Push(SideConnectLocId(side,is));
  }
}

void TPZTriangle::HigherDimensionSides(int side, TPZStack<int> &high)
{
	if(side <0 || side >= NSides) {
		PZError << "TPZTriangle::HigherDimensionSides side "<< side << endl;
	}
	int is;
	for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
	
}


int TPZTriangle::SideDimension(int side) {
	if(side<0 || side >= NSides) {
		PZError << "TPZTriangle::SideDimension side " << side << endl;
		return -1;
	}
	return sidedimension[side];
}

void TPZTriangle::CenterPoint(int side, TPZVec<REAL> &center) {
  //center.Resize(Dimension);
  int i;
  for(i=0; i<Dimension; i++) {
    center[i] = MidSideNode[side][i];
  }
}

TPZTransform TPZTriangle::TransformElementToSide(int side){

	if(side<0 || side>6){
  	PZError << "TPZTriangle::TransformElementToSide called with side error\n";
    return TPZTransform(0,0);
  }

  TPZTransform t(sidedimension[side],2);//t(dimto,2)
  t.Mult().Zero();
  t.Sum().Zero();

  switch(side){
    case 0:
    case 1:
    case 2:
      return t;
    case 3:
    	t.Mult()(0,0) =  2.0;//par. var.
      t.Sum()(0,0)  = -1.0;
      return t;
    case 4:
      t.Mult()(0,0) = -1.0;
      t.Mult()(0,1) =  1.0;
      return t;
    case 5:
    	t.Mult()(0,1) = -2.0;
      t.Sum()(0,0)  =  1.0;
      return t;
    case 6:
      t.Mult()(0,0) =  1.0;
    	t.Mult()(1,1) =  1.0;
      return t;
  }
  return TPZTransform(0,0);
}

TPZTransform TPZTriangle::TransformSideToElement(int side){

	if(side<0 || side>6){
  	PZError << "TPZTriangle::TransformSideToElement side out range\n";
    return TPZTransform(0,0);
  }
  TPZTransform t(2,sidedimension[side]);
  t.Mult().Zero();
  t.Sum().Zero();

  switch(side){
  	case 0:
      return t;
    case 1:
    	t.Sum()(0,0) =  1.0;
      return t;
    case 2:
      t.Sum()(1,0) =  1.0;
      return t;
    case 3:
    	t.Mult()(0,0) =  0.5;
      t.Sum()(0,0)  =  0.5;
      return t;
    case 4:
      t.Mult()(0,0) = -0.5;
      t.Mult()(1,0) =  0.5;
      t.Sum() (0,0) =  0.5;
      t.Sum() (1,0) =  0.5;
      return t;
    case 5:
      t.Mult()(1,0) = -0.5;
      t.Sum() (1,0) =  0.5;
      return t;
    case 6:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      return t;
  }
	return TPZTransform(0,0);
}

TPZIntPoints * TPZTriangle::CreateSideIntegrationRule(int side, int order){
	if(side < 0 || side>6) {
		PZError << "TPZTriangle::CreateSideIntegrationRule wrong side " << side << endl;
		return 0;
	}
	if(side<3) return new TPZInt1Point();
	if(side<6) return new TPZInt1d(order);
	if(side==6)return new TPZIntTriang(order);
	return 0;
}


MElementType TPZTriangle::Type()
{
  return ETriangle;
}

MElementType TPZTriangle::Type(int side)
{
  switch(side) {
    case 0:
    case 1:
    case 2:
      return EPoint;
    case 3:
    case 4:
    case 5:
      return EOned;
    case 6:
      return ETriangle;
    default:
      return ENoType;
  }
}


int TPZTriangle::NConnects() {
	return NSides;
}



int TPZTriangle::NSideConnects(int side) {
  if(side<0 || side>6) {
    PZError << "TPZShapeTriang::NSideConnects. Bad parameter i.\n";
    return 0;
  }
  if(side<3) return 1;
  if(side<6) return 3;
  return 7;
}

/**It do not verify the values of the c*/
int TPZTriangle::SideConnectLocId(int side, int c) {
  switch(side) {
  case 0:
  case 1:
  case 2:
    return side;
  case 3:
  case 4:
  case 5:
    if(!c) return side-3;
    if(c==1) return (side-2)%3;
    if(c==2) return side;
  case 6:
    return c;
  default:
    PZError << "TPZShapeTriang::SideConnectLocId, connect = " << c << endl;
    return -1;
  }
}
}
