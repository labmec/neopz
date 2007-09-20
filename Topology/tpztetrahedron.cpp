//
// C++ Implementation: tpztetrahedron
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpztetrahedron.h"

#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"
#include "pzquad.h"
#include "pzeltype.h"

using namespace std;

namespace pztopology {

TPZTetrahedron::TPZTetrahedron()
{
}

TPZTetrahedron::~TPZTetrahedron()
{
}

static int nhighdimsides[15] = {7,7,7,7,3,3,3,3,3,3,1,1,1,1,0};

int TPZTetrahedron::FaceNodes[4][3]  = { {0,1,2},{0,1,3},{1,2,3},{0,2,3} };

int TPZTetrahedron::SideNodes[6][2]  = { {0,1},{1,2},{2,0},{0,3},{1,3},{2,3} };

int TPZTetrahedron::ShapeFaceId[4][3] = { {0,1,2},{0,1,3},{1,2,3},{0,2,3} };

static int sidedimension[15] = {0,0,0,0,1,1,1,1,1,1,2,2,2,2,3};


static int FaceConnectLocId[4][7] = { {0,1,2,4,5,6,10},{0,1,3,4,8,7,11},
					     {1,2,3,5,9,8,12},{0,2,3,6,9,7,13} };

static int highsides[15][7] = {
{4,6,7,10,11,13,14},
{4,5,8,10,11,12,14},
{5,6,9,10,12,13,14},
{7,8,9,11,12,13,14},
{10,11,14},
{10,12,14},
{10,13,14},
{11,13,14},
{11,12,14},
{12,13,14},
{14},
{14},
{14},
{14},
{-999}
};

static int nsidenodes[15] = 
{
	1,1,1,1,
	2,2,2,2,2,2,
	3,3,3,3,
	4};

static REAL sidetosidetransforms[15][7][4][3] = {
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,0}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,0}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,0}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,1}}
},
{
{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{0.5,0,0},{-99,-99,-99},{-99,-99,-99},{0.5,0,0}}
},
{
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{-0.5,0.5,0},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,0}}
},
{
{{0,-0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{-0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{0,-0.5,0},{-99,-99,-99},{-99,-99,-99},{0,0.5,0}}
},
{
{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{0,0,0.5},{-99,-99,-99},{-99,-99,-99},{0,0,0.5}}
},
{
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{-0.5,0,0.5},{-99,-99,-99},{-99,-99,-99},{0.5,0,0.5}}
},
{
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{0,-0.5,0.5},{-99,-99,-99},{-99,-99,-99},{0,0.5,0.5}}
},
{
{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,0}}
},
{
{{1,0,0},{0,0,1},{-99,-99,-99},{0,0,0}}
},
{
{{-1,1,0},{-1,0,1},{-99,-99,-99},{1,0,0}}
},
{
{{0,1,0},{0,0,1},{-99,-99,-99},{0,0,0}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
}
};

static REAL MidSideNode[15][3] = {
/*00*/{.0,.0},/*01*/{1.,.0},/*02*/{0.,1.,.0},/*03*/{.0,0.,1.0},/*04*/{.5,.0,.0},
/*05*/{.5,.5},/*06*/{0.,.5},/*07*/{0.,0.,.5},/*08*/{.5,0.,0.5},/*09*/{.0,.5,.5},
/*10*/{1./3.,1./3., 0.  }  ,/*11*/{1./3., .0  ,1./3.},
/*12*/{1./3.,1./3.,1./3.}  ,/*13*/{ 0.  ,1./3.,1./3.},/*14*/{1./4.,1./4.,1./4.} };

void TPZTetrahedron::LowerDimensionSides(int side,TPZStack<int> &smallsides)
{
     smallsides.Resize(0);
     int nsidecon = NSideConnects(side);
     int is;
     for(is=0; is<nsidecon-1; is++)
     smallsides.Push(SideConnectLocId(side,is));
}

void TPZTetrahedron::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
{
     smallsides.Resize(0);
     int nsidecon = NSideConnects(side);
     for(int is = 0; is < nsidecon - 1; is++) {
     if (SideDimension(SideConnectLocId(side,is)) == DimTarget) smallsides.Push(SideConnectLocId(side,is));
  }
}

void TPZTetrahedron::HigherDimensionSides(int side, TPZStack<int> &high)
{
	if(side <0 || side >= NSides) {
		PZError << "TPZTetrahedron::HigherDimensionSides side "<< side << endl;
	}
	int is;
	for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
	
}

int TPZTetrahedron::NSideNodes(int side)
{
	return nsidenodes[side];
}

int TPZTetrahedron::SideNodeLocId(int side, int node)
{
	if(side<4 && node == 0) return side;
	if(side>=4 && side < 10 && node <2) return SideNodes[side-4][node];
	if(side >= 10 && side < 14 && node <3) return FaceNodes[side-10][node];
	if(side ==14 && node < 4) return node;
	PZError << "TPZTetrahedron::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
	return -1;

}

void TPZTetrahedron::CenterPoint(int side, TPZVec<REAL> &center) {
  //center.Resize(Dimension);
  int i;
  for(i=0; i<Dimension; i++) {
    center[i] = MidSideNode[side][i];
  }
}

int TPZTetrahedron::SideDimension(int side) {
	if(side<0 || side >= NSides) {
		PZError << "TPZTetrahedron::SideDimension side " << side << endl;
		return -1;
	}
	return sidedimension[side];
}

TPZTransform TPZTetrahedron::SideToSideTransform(int sidefrom, int sideto)
{
	if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
		PZError << "TPZTetrahedron::HigherDimensionSides sidefrom "<< sidefrom << 
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
	PZError << "TPZTetrahedron::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
	return TPZTransform(0);
}

TPZTransform TPZTetrahedron::TransformElementToSide(int side){

	if(side<0 || side>14){
  	PZError << "TPZTetrahedron::TransformElementToSide called with side error\n";
    return TPZTransform(0,0);
  }

  TPZTransform t(sidedimension[side],3);
  t.Mult().Zero();
  t.Sum().Zero();

  switch(side){
    case 0:
    case 1:
    case 2:
    case 3:
      return t;
    case 4:
    	t.Mult()(0,0) =  2.0;
      t.Sum()(0,0)  = -1.0;
      return t;
    case 5:
      t.Mult()(0,0) = -1.0;
      t.Mult()(0,1) =  1.0;
      return t;
    case 6:
    	t.Mult()(0,1) = -2.0;
      t.Sum()(0,0)  =  1.0;
      return t;
    case 7:
    	t.Mult()(0,2) =  2.0;
      t.Sum()(0,0)  = -1.0;
      return t;
    case 8:
    	t.Mult()(0,0) = -1.0;
      t.Mult()(0,2) =  1.0;
      return t;
    case 9:
      t.Mult()(0,1) = -1.0;
      t.Mult()(0,2) =  1.0;
      return t;
    case 10:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      return t;
    case 11:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,2) =  1.0;
      return t;
    case 12:
    case 13:
      t.Mult()(0,1) =  1.0;
      t.Mult()(1,2) =  1.0;
      return t;
    case 14:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      t.Mult()(2,2) =  1.0;
      return t;
  }
  return TPZTransform(0,0);
}

TPZTransform TPZTetrahedron::TransformSideToElement(int side){

	if(side<0 || side>14){
  	PZError << "TPZTetrahedron::TransformSideToElement side out range\n";
    return TPZTransform(0,0);
  }
  TPZTransform t(3,sidedimension[side]);
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
    	t.Sum()(2,0) =  1.0;
      return t;
    case 4:
      t.Mult()(0,0) =  0.5;
      t.Sum() (0,0) =  0.5;
      return t;
    case 5:
      t.Mult()(0,0) = -0.5;
      t.Mult()(1,0) =  0.5;
      t.Sum() (0,0) =  0.5;
      t.Sum() (1,0) =  0.5;
      return t;
    case 6:
      t.Mult()(1,0) = -0.5;
      t.Sum() (1,0) =  0.5;
      return t;
    case 7:
      t.Mult()(2,0) =  0.5;
      t.Sum() (2,0) =  0.5;
      return t;
    case 8:
      t.Mult()(0,0) = -0.5;
      t.Mult()(2,0) =  0.5;
      t.Sum() (0,0) =  0.5;
      t.Sum() (2,0) =  0.5;
      return t;
    case 9:
      t.Mult()(1,0) = -0.5;
      t.Mult()(2,0) =  0.5;
      t.Sum() (1,0) =  0.5;
      t.Sum() (2,0) =  0.5;
      return t;
    case 10:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      return t;
    case 11:
      t.Mult()(0,0) =  1.0;
      t.Mult()(2,1) =  1.0;
      return t;
    case 12:
      t.Mult()(0,0) = -1.0;
      t.Mult()(0,1) = -1.0;
      t.Mult()(1,0) =  1.0;
      t.Mult()(2,1) =  1.0;
      t.Sum() (0,0) =  1.0;
      return t;
    case 13:
      t.Mult()(1,0) =  1.0;
      t.Mult()(2,1) =  1.0;
      return t;
    case 14:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      t.Mult()(2,2) =  1.0;
      return t;

  }
	return TPZTransform(0,0);
}

TPZIntPoints * TPZTetrahedron::CreateSideIntegrationRule(int side, int order) {
	if(side<0 || side>15) {
		PZError << "TPZTetrahedron::CreateSideIntegrationRule. bad side number.\n";
   		return 0;
	}
	//SideOrder corrige sides de 4 a 14 para 0 a 10
	if(side<4)   return new TPZInt1Point();//cantos 0 a 3 : cria regra com um ponto
	if(side<10)  return new TPZInt1d(order);//lados 4 a 9
	if(side<14)  {//faces : 10 a 13
   		return new TPZIntTriang(order);
	}
	if(side==14) {//integra�o do elemento
   		return new TPZIntTetra3D(order);
	}
	return 0;
}


MElementType TPZTetrahedron::Type()
{
  return ETetraedro;
}

MElementType TPZTetrahedron::Type(int side)
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
    case 8:
    case 9:
      return EOned;
    case 10:
    case 11:
    case 12:
    case 13:
      return ETriangle;
    case 14:
      return ETetraedro;
    default:
      return ENoType;
  }
}


int TPZTetrahedron::NConnects() {
	return NSides;
}


int TPZTetrahedron::NSideConnects(int side) {
	if(side<0)   return -1;
	if(side<4)   return 1;//cantos : 0 a 3
   if(side<10)  return 3;//lados : 4 a 9
   if(side<14)	 return 7;//faces : 10 a 13
   if(side==14) return 15;//centro : 15
   return -1;
}

int TPZTetrahedron::SideConnectLocId(int side, int node) {
	if(side<0 || side>15) return -1;
   if(side<4) {
   	if(node==0) return side;
   } else
   if(side<7) {//4,5,6
	   int s = side-4;//0,1,2
   	if(!node) return s;//0,1,2
      if(node==1) return (s+1)%3;//1,2,0
      if(node==2) return side;//4,5,6
   } else
   if(side<10) {//7,8,9
   	int s = side-7;//0,1,2
   	if(!node) return s;//0,1,2
      if(node==1) return 3;//4,4,4
      if(node==2) return side;//7,8,9
   } else
   if(side<14) {//10 a 13
   	int s = side-10;
   	if(node<7) return FaceConnectLocId[s][node];
   } else
   if(side==14 && node<15){
   	return node;
   }
	PZError << "TPZShapeTetra::SideConnectLocId called for node = "
   	      << node << " and side = " << side << "\n";
   return -1;
}


}
