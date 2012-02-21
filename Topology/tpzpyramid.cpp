/**
 * @file
 * @brief Contains the implementation of the TPZPyramid methods into the pztopology scope.
 */

#include "tpzpyramid.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"
#include "pzeltype.h"

#include "pzcreateapproxspace.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.topology.pzpyramid"));
#endif

using namespace std;

namespace pztopology {	

	static int FaceConnectLocId[5][9] = { {0,1,2,3,5,6,7,8,13},{0,1,4,5,10,9,14,-1,-1},
		{1,2,4,6,11,10,15,-1,-1},{3,2,4,7,11,12,16,-1,-1},{0,3,4,8,12,9,17,-1,-1} };
	
	
	static int nhighdimsides[19] = {7,7,7,7,9,3,3,3,3,3,3,3,3,1,1,1,1,1,0};
	
	int TPZPyramid::SideNodes[8][2]  = { {0,1},{1,2},{2,3},{3,0},{0,4},{1,4},{2,4},{3,4} };
	
	int TPZPyramid::FaceNodes[5][4]  = { {0,1,2,3},{0,1,4,-1},{1,2,4,-1},{3,2,4,-1},{0,3,4,-1} };
	
	int TPZPyramid::ShapeFaceId[5][4] = { {0,1,2,3},{0,1,4,-1},{1,2,4,-1},{3,2,4,-1},{0,3,4,-1} };
	
	static int sidedimension[19] = {0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,3};
	
	static int highsides[19][9] = {
		{5,8,9,13,14,17,18},
		{5,6,10,13,14,15,18},
		{6,7,11,13,15,16,18},
		{7,8,12,13,16,17,18},
		{9,10,11,12,14,15,16,17,18},
		{13,14,18},
		{13,15,18},
		{13,16,18},
		{13,17,18},
		{14,17,18},
		{14,15,18},
		{15,16,18},
		{16,17,18},
		{18},
		{18},
		{18},
		{18},
		{18},
		{-999}
	};
	
	int nsidenodes[19] = {1,1,1,1,1,
		2,2,2,2,2,2,2,2,
		4,3,3,3,3,
		5};
	
	int TPZPyramid::NSideNodes(int side)
	{
		return nsidenodes[side];
	}
	
	static REAL sidetosidetransforms[19][9][4][3] = {
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,0}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,0}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,0}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,0}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,1}}
		},
		{
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
			{{1,0,0},{-99,-99,-99},{-99,-99,-99},{0,-1,0}}
		},
		{
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
			{{0,1,0},{-99,-99,-99},{-99,-99,-99},{1,0,0}}
		},
		{
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
			{{-1,0,0},{-99,-99,-99},{-99,-99,-99},{0,1,0}}
		},
		{
			{{0,-1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{-0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
			{{0,-1,0},{-99,-99,-99},{-99,-99,-99},{-1,0,0}}
		},
		{
			{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
			{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
			{{0.5,0.5,0.5},{-99,-99,-99},{-99,-99,-99},{-0.5,-0.5,0.5}}
		},
		{
			{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
			{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
			{{-0.5,0.5,0.5},{-99,-99,-99},{-99,-99,-99},{0.5,-0.5,0.5}}
		},
		{
			{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
			{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
			{{-0.5,-0.5,0.5},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,0.5}}
		},
		{
			{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
			{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
			{{0.5,-0.5,0.5},{-99,-99,-99},{-99,-99,-99},{-0.5,0.5,0.5}}
		},
		{
			{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,0}}
		},
		{
			{{2,0,0},{1,1,1},{-99,-99,-99},{-1,-1,0}}
		},
		{
			{{0,2,0},{-1,1,1},{-99,-99,-99},{1,-1,0}}
		},
		{
			{{2,0,0},{1,-1,1},{-99,-99,-99},{-1,1,0}}
		},
		{
			{{0,2,0},{1,1,1},{-99,-99,-99},{-1,-1,0}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
		}
	};
	
	static REAL MidSideNode[19][3] = {
		/*00*/{-1.,-1.},   /*01*/{1.,-1.},   /*02*/{1.,1.},/*03*/{-1.,1.},/*04*/{0.,0.,1.},
		/*05*/{ 0.,-1.},   /*06*/{1., 0.},   /*07*/{0.,1.},/*08*/{-1.,0.},
		/*09*/{-.5,-.5,.5},/*10*/{.5,-.5,.5},/*11*/{.5,.5,.5},/*12*/{-.5,.5,.5},
		/*13*/{0.,  0. ,  0. },/*14*/{  0.  ,-2./3.,1./3.},/*15*/{2./3.,0.,1./3.},
		/*16*/{0.,2./3.,1./3.},/*17*/{-2./3.,  0.  ,1./3.},/*18*/{  0. ,0.,1./5.} };
	
	void TPZPyramid::LowerDimensionSides(int side,TPZStack<int> &smallsides)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		int is;
		for(is=0; is<nsidecon-1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	void TPZPyramid::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
		}
	}
	
	void TPZPyramid::HigherDimensionSides(int side, TPZStack<int> &high)
	{
		if(side <0 || side >= NSides) {
			PZError << "TPZPyramid::HigherDimensionSides side "<< side << endl;
		}
		int is;
		for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
		
	}
	
	int TPZPyramid::SideNodeLocId(int side, int node)
	{
		if(side <5 && node == 0) return side;
		if(side >= 5 && side <13 && node < 2) return SideNodes[side-5][node];
		if(side == 13 && node <4) return FaceNodes[side-13][node];
		if(side >13 && side < 18)
			if (node <3) return FaceNodes[side-13][node];
			else if (node==3) return -1;//Previsto receber pelas faces triangulares - Cesar 2003-01-02
		
		if(side == 18 && node < 5) return node;
		PZError << "TPZPyramid::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
		return -1;
	}
	
	void TPZPyramid::CenterPoint(int side, TPZVec<REAL> &center) {
		//center.Resize(Dimension);
		int i;
		for(i=0; i<Dimension; i++) {
			center[i] = MidSideNode[side][i];
		}
	}
	
	int TPZPyramid::SideDimension(int side) {
		if(side<0 || side >= NSides) {
			PZError << "TPZPyramid::SideDimension side " << side << endl;
			return -1;
		}
		return sidedimension[side];
	}
	
	TPZTransform TPZPyramid::SideToSideTransform(int sidefrom, int sideto)
	{
		if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
			PZError << "TPZPyramid::HigherDimensionSides sidefrom "<< sidefrom << 
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
		PZError << "TPZPyramid::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
		return TPZTransform(0);
	}
	
	TPZTransform TPZPyramid::TransformElementToSide(int side){
		
		if(side<0 || side>18){
			PZError << "TPZPyramid::TransformElementToSide called with side error\n";
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
			case 4:
				return t;
			case 5:
				t.Mult()(0,0) = 1.0;
				return t;
			case 6:
				t.Mult()(0,1) = 1.0;
				return t;
			case 7:
				t.Mult()(0,0) = -1.0;
				return t;
			case 8:
				t.Mult()(0,1) = -1.0;
				return t;
			case 9:
			case 12:
				t.Mult()(0,0) = 2.0;
				t.Sum()(0,0)  = 1.0;
				return t;
			case 10:
			case 11:
				t.Mult()(0,0) = -2.0;
				t.Sum()(0,0)  =  1.0;
				return t;
			case 13:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				return t;
			case 14:
				t.Mult()(0,0) =  0.5;
				t.Mult()(0,1) = -0.5;
				t.Mult()(1,2) =  1.0;
				return t;
			case 15:
			case 16:/** CONTEM ERRO AQUI */
				t.Mult()(0,0) =  0.5;
				t.Mult()(0,1) =  0.5;
				t.Mult()(1,2) =  1.0;
				return t;
			case 17:
				t.Mult()(0,0) = -0.5;
				t.Mult()(0,1) =  0.5;
				t.Mult()(1,2) =  1.0;
				return t;
			case 18:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,2) =  1.0;
				return t;
		}
		return TPZTransform(0,0);
	}
	
	TPZTransform TPZPyramid::TransformSideToElement(int side){
		
		if(side<0 || side>18){
			PZError << "TPZPyramid::TransformSideToElement side out range\n";
			return TPZTransform(0,0);
		}
		TPZTransform t(3,sidedimension[side]);
		t.Mult().Zero();
		t.Sum().Zero();
		
		switch(side){
			case 0:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) = -1.0;
				t.Sum()(2,0) =  0.0;
				return t;
			case 1:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) = -1.0;
				t.Sum()(2,0) =  0.0;
				return t;
			case 2:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) =  0.0;
				return t;
			case 3:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) =  0.0;
				return t;
			case 4:
				t.Sum()(0,0) =  0.0;
				t.Sum()(1,0) =  0.0;
				t.Sum()(2,0) =  1.0;
				return t;
			case 5:
				t.Mult()(0,0) =  1.0;
				t.Sum() (1,0) = -1.0;
				return t;
			case 6:
				t.Mult()(1,0) =  1.0;
				t.Sum() (0,0) =  1.0;
				return t;
			case 7:
				t.Mult()(0,0) = -1.0;
				t.Sum() (1,0) =  1.0;
				return t;
			case 8:
				t.Mult()(1,0) = -1.0;
				t.Sum()(0,0)  = -1.0;
				return t;
			case 9:
				t.Mult()(0,0) =  0.5;
				t.Mult()(1,0) =  0.5;
				t.Mult()(2,0) =  0.5;
				t.Sum()(0,0)  = -0.5;
				t.Sum()(1,0)  = -0.5;
				t.Sum()(2,0)  =  0.5;
				return t;
			case 10:
				t.Mult()(0,0) = -0.5;
				t.Mult()(1,0) =  0.5;
				t.Mult()(2,0) =  0.5;
				t.Sum()(0,0)  =  0.5;
				t.Sum()(1,0)  = -0.5;
				t.Sum()(2,0)  =  0.5;
				return t;
			case 11:
				t.Mult()(0,0) = -0.5;
				t.Mult()(1,0) = -0.5;
				t.Mult()(2,0) =  0.5;
				t.Sum()(0,0)  =  0.5;
				t.Sum()(1,0)  =  0.5;
				t.Sum()(2,0)  =  0.5;
				return t;
			case 12:
				t.Mult()(0,0) =  0.5;
				t.Mult()(1,0) = -0.5;
				t.Mult()(2,0) =  0.5;
				t.Sum()(0,0)  = -0.5;
				t.Sum()(1,0)  =  0.5;
				t.Sum()(2,0)  =  0.5;
				return t;
			case 13:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				return t;
			case 14:
				t.Mult()(0,0) =  2.0;
				t.Mult()(0,1) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(1,0)  = -1.0;
				return t;
			case 15:
				t.Mult()(1,0) =  2.0;
				t.Mult()(0,1) = -1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(0,0)  =  1.0;
				t.Sum()(1,0)  = -1.0;
				return t;
			case 16:
				t.Mult()(0,0) =  2.0;
				t.Mult()(0,1) =  1.0;
				t.Mult()(1,1) = -1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(1,0)  =  1.0;
				return t;
			case 17:
				t.Mult()(1,0) =  2.0;
				t.Mult()(0,1) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(1,0)  = -1.0;
				return t;
			case 18:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,2) =  1.0;
				return t;
		}
		return TPZTransform(0,0);
	}
	
	
	
	TPZIntPoints * TPZPyramid::CreateSideIntegrationRule(int side, int order){
		if(side<0 || side>18) {
			PZError << "TPZPyramid::CreateSideIntegrationRule. Bad side number.\n";
			return 0;
		}
		if(side<5)   return new TPZInt1Point(order);            // sides 0 to 4 are vertices
		if(side<13)  return new TPZInt1d(order);           // sides 5 to 12 are lines
		if(side==13) return new TPZIntQuad(order,order);   // side 13 are quadrilateral (pyramid base)
		if(side<18)  {
			return new TPZIntTriang(order);                // sides 14 to 17 are triangles
		}
		if(side==18) {
			return new IntruleType(order);               // integration of the element
		}
		return 0;
	}

	MElementType TPZPyramid::Type()
	{
		return EPiramide;
	}

	MElementType TPZPyramid::Type(int side)
	{
		switch(side) {
			case 0:
			case 1:
			case 2:
			case 3:
			case 4:
				return EPoint;
			case 5:
			case 6:
			case 7:
			case 8:
			case 9:
			case 10:
			case 11:
			case 12:
				return EOned;
			case 13:
				return EQuadrilateral;
			case 14:
			case 15:
			case 16:
			case 17:
				return ETriangle;
			case 18:
				return EPiramide;
			default:
				return ENoType;
		}
	}
	
	
	int TPZPyramid::NumSides() {
		return 19;
	}
	
	//Tentando criar o metodo
	int TPZPyramid::NumSides(int dimension) {
		if(dimension<0 || dimension> 3) {
			PZError << "TPZPyramid::NumSides. Bad parameter i.\n";
			return 0;
		}
		if(dimension==0) return 5;
		if(dimension==1) return 8;
		if(dimension==2) return 5;
		if(dimension==3) return 1;
		return -1;
	}
	int TPZPyramid::NContainedSides(int side) {
		if(side<0)   return -1;
		if(side<5)   return 1;//cantos : 0 a 4
		if(side<13)  return 3;//lados : 5 a 12
		if(side==13) return 9;//face : 13 , quadrilateral
		if(side<18)	 return 7;//faces : 14 a 17 , triangulares
		if(side==18) return 19;//centro : 18
		return -1;
	}
	
	int TPZPyramid::ContainedSideLocId(int side, int node) {
		if(side<0 || side>19 || node < 0) return -1;
		if(side<5) {
			if(node==0) return side;
		} 
		else if(side<9) {//5 a 8
			int s = side-5;//0,1,2
			if(!node) return s;//0,1,2,3
			if(node==1) return (s+1)%4;//1,2,0
			if(node==2) return side;//5,6,7,8
		}
		else if(side<13) {//9 a 12
			int s = side-9;//0,1,2,3
			if(!node) return s;//0,1,2,3
			if(node==1) return 4;//
			if(node==2) return side;//9,10,11,12
		} 
		else if(side<18) {//13 a 17
			int s = side-13;
			if(node<9) return FaceConnectLocId[s][node];
		} 
		else if(side==18 && node<19){
			return node;
		}
		PZError << "TPZShapePiram::ContainedSideLocId called for node = "
		<< node << " and side = " << side << "\n";
		return -1;
	}
	
	bool TPZPyramid::IsInParametricDomain(TPZVec<REAL> &pt, REAL tol){
		const REAL qsi = pt[0];
		const REAL eta = pt[1];
		const REAL zeta = pt[2];
		
		if( (qsi < -1. - tol) || (qsi > 1.+tol) ||
		   (eta < -1. - tol) || (eta > 1.+tol) || 
		   (zeta < 0. - tol) || (zeta > 1.+tol) ||
		   (fabs(qsi) > 1.-zeta + tol) || (fabs(eta) > 1.-zeta + tol)
		   ) {
			return false;
		}
		else{
			return true;
		}  
		
		
	}//method
	
	
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
	void TPZPyramid::GetSideHDivPermutation(int side, TPZVec<int> &id, TPZVec<int> &permgather)
	{
		std::cout << "Please implement me\n";
		DebugStop();
	}
}
