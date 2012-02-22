/**
 * @file
 * @brief Contains the implementation of the TPZPrism methods. 
 */

#include "tpzprism.h"

#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapelinear.h"
#include "pzshapepoint.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"
#include "pzquad.h"
#include "pzeltype.h"

#include "pzcreateapproxspace.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.topology.pzprism"));
#endif

using namespace std;

namespace pztopology {

	static int FaceConnectLocId[5][9] = { {0,1,2,6,7,8,15,-1,-1},{0,1,4,3,6,10,12,9,16},
		{1,2,5,4,7,11,13,10,17},{0,2,5,3,8,11,14,9,18},{3,4,5,12,13,14,19,-1,-1} };
	
	
	int TPZPrism::FaceNodes[5][4]  = { {0,1,2,-1},{0,1,4,3},{1,2,5,4},{0,2,5,3},{3,4,5,-1} };
	//F15        F16       F17       F18        F19
	
	int TPZPrism::SideNodes[9][2]  = { {0,1},{1,2},{2,0},{0,3},{1,4},{2,5},{3,4},{4,5},{5,3} };
	//arestas   6     7      8    9     10    11    12    13    14
	
	int TPZPrism::ShapeFaceId[5][4] = { {0,1,2,-1},{0,1,4,3},{1,2,5,4},{0,2,5,3},{3,4,5,-1} };
	//F15        F16       F17       F18       F19
	
	static int sidedimension[21] = {0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,3};
	
	static int nhighdimsides[21] = {7,7,7,7,7,7,3,3,3,3,3,3,3,3,3,1,1,1,1,1,0};
	
	static int highsides[21][7] = {
		{6,8,9,15,16,18,20},
		{6,7,10,15,16,17,20},
		{7,8,11,15,17,18,20},
		{9,12,14,16,18,19,20},
		{10,12,13,16,17,19,20},
		{11,13,14,17,18,19,20},
		{15,16,20},
		{15,17,20},
		{15,18,20},
		{16,18,20},
		{16,17,20},
		{17,18,20},
		{16,19,20},
		{17,19,20},
		{18,19,20},
		{20},
		{20},
		{20},
		{20},
		{20},
		{-999}
	};
	
	static REAL sidetosidetransforms[21][7][4][3] = {
		{
			//0
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-1}}
		},
		{
			//1
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			//{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},//est� errado deve ser {-1,-1,-99}
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},//CEDRIC
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-1}}
		},
		{
			//2
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-1}}
		},
		{
			//3
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,1}}
		},
		{
			//4
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,1}}
		},
		{
			//5
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,1}}
		},
		{
			//6
			{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{0.5,0,0},{-99,-99,-99},{-99,-99,-99},{0.5,0,-1}}
		},
		{
			//7
			{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{-0.5,0.5,0},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-1}}
		},
		{
			//8
			{{0,-0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{0,-0.5,0},{-99,-99,-99},{-99,-99,-99},{0,0.5,-1}}
		},
		{
			//9
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{0,0,0}}
		},
		{
			//10
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{1,0,0}}
		},
		{
			//11
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{0,1,0}}
		},
		{
			//12
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
			{{0.5,0,0},{-99,-99,-99},{-99,-99,-99},{0.5,0,1}}
		},
		{
			//13
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
			{{-0.5,0.5,0},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,1}}
		},
		{
			//14
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{0,-0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
			{{0,-0.5,0},{-99,-99,-99},{-99,-99,-99},{0,0.5,1}}
		},
		{
			//15
			{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,-1}}
		},
		{
			//16
			{{0.5,0,0},{0,0,1},{-99,-99,-99},{0.5,0,0}}
		},
		{
			//17
			{{-0.5,0.5,0},{0,0,1},{-99,-99,-99},{0.5,0.5,0}}
		},
		{
			//18
			{{0,0.5,0},{0,0,1},{-99,-99,-99},{0,0.5,0}}
		},
		{
			//19
			{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,1}}
		},
		{
			//20
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
		}
	};
	
	static REAL MidSideNode[21][3] = {
		/*00*/{0.,.0,-1.},/*01*/{1.,0.,-1.},/*02*/{.0,1.,-1.},/*03*/{.0,.0, 1.},
		/*04*/{1.,.0, 1.},/*05*/{0.,1., 1.},/*06*/{.5,.0,-1.},/*07*/{.5,.5,-1.},
		/*08*/{.0,.5,-1.},/*09*/{0.,.0, 0.},/*10*/{1.,.0, 0.},/*11*/{.0,1., 0.},
		/*12*/{.5,.0, 1.},/*13*/{.5,.5, 1.},/*14*/{.0,.5, 1.},/*15*/{1./3.,1./3.,-1.},
		/*16*/{.5,.0, 0.},/*17*/{.5,.5, 0.},/*18*/{0.,.5, 0.},/*19*/{1./3.,1./3., 1.},
		/*20*/{1./3.,1./3.,0.} };
	
	void TPZPrism::LowerDimensionSides(int side,TPZStack<int> &smallsides)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		int is;
		for(is=0; is<nsidecon-1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	void TPZPrism::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
		}
	}
	
	void TPZPrism::HigherDimensionSides(int side, TPZStack<int> &high)
	{
		if(side <0 || side >= NSides) {
			PZError << "TPZPrism::HigherDimensionSides side "<< side << endl;
		}
		int is;
		for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
		
	}
	
	static int nsidenodes[21] = {
		1,1,1,1,1,1,
		2,2,2,2,2,2,2,2,2,
		3,4,4,4,3,
		6};
	
	int TPZPrism::NSideNodes(int side)
	{
		return nsidenodes[side];
	}
	//Tentando criar o metodo
	int TPZPrism::NumSides(int dimension) {
		if(dimension<0 || dimension> 3) {
			PZError << "TPZPrism::NumSides. Bad parameter i.\n";
			return 0;
		}
		if(dimension==0) return 6;
		if(dimension==1) return 9;
		if(dimension==2) return 5;
		if(dimension==3) return 1;
		return -1;
	}
	
	int TPZPrism::SideNodeLocId(int side, int node)
	{
		if(side<6 && node == 0) return side;
		if(side>= 6 && side < 15 && node<2) return SideNodes[side-6][node];
		if(side==15)
			if (node < 3) return FaceNodes[side-15][node];
			else if (node == 3) return -1; //previsto para faces triangulares
		
		if(side>15 && side <19 && node <4) return FaceNodes[side-15][node];
		if(side==19)
			if (node<3) return FaceNodes[side-15][node];
			else if (node == 3) return -1; // Previsto p/ faces triangulares
		
		if(side==20 && node<6) return node;
		PZError << "TPZPrism::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
		return -1;
	}
	
	void TPZPrism::CenterPoint(int side, TPZVec<REAL> &center) {
		center.Resize(Dimension);
		int i;
		for(i=0; i<Dimension; i++) {
			center[i] = MidSideNode[side][i];
		}
	}
	
	int TPZPrism::SideDimension(int side) {
		if(side<0 || side >= NSides) {
			PZError << "TPZPrism::SideDimension side " << side << endl;
			return -1;
		}
		return sidedimension[side];
	}
	
	TPZTransform TPZPrism::SideToSideTransform(int sidefrom, int sideto)
	{
		if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
			PZError << "TPZPrism::HigherDimensionSides sidefrom "<< sidefrom << 
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
		PZError << "TPZPrism::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
		return TPZTransform(0);
	}
	
	
	TPZTransform TPZPrism::TransformElementToSide(int side){
		
		if(side<0 || side>20){
			PZError << "TPZPrism::TransformElementToSide called with side error\n";
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
			case 5:
				return t;
			case  6:
			case 12:
				t.Mult()(0,0) =  2.0;
				t.Sum()(0,0)  = -1.0;
				return t;
			case  7:
			case 13:
				t.Mult()(0,0) = -1.0;
				t.Mult()(0,1) =  1.0;
				return t;
			case  8:
			case 14:
				t.Mult()(0,1) = -2.0;
				t.Sum()(0,0)  =  1.0;
				return t;
			case  9:
			case 10:
			case 11:
				t.Mult()(0,2) =  1.0;
				return t;
			case 15:
			case 19:
				t.Mult()(0,0) = 1.0;
				t.Mult()(1,1) = 1.0;
				return t;
			case 16:
				t.Mult()(0,0) =  2.0;
				t.Mult()(1,2) =  1.0;
				t.Sum()(0,0)  = -1.0;
				return t;
			case 17:
				t.Mult()(0,0) = -1.0;
				t.Mult()(0,1) =  1.0;
				t.Mult()(1,2) =  1.0;
				return t;
			case 18:
				t.Mult()(0,1) =  2.0;
				t.Mult()(1,2) =  1.0;
				t.Sum()(0,0)  = -1.0;
				return t;
			case 20:
				t.Mult()(0,0) = 1.0;
				t.Mult()(1,1) = 1.0;
				t.Mult()(2,2) = 1.0;
				return t;
		}
		return TPZTransform(0,0);
	}
	
	TPZTransform TPZPrism::TransformSideToElement(int side){
		
		if(side<0 || side>20){
			PZError << "TPZPrism::TransformSideToElement side out range\n";
			return TPZTransform(0,0);
		}
		TPZTransform t(3,sidedimension[side]);
		t.Mult().Zero();
		t.Sum().Zero();
		
		switch(side){
			case 0:
				t.Sum()(0,0) =  0.0;
				t.Sum()(1,0) =  0.0;
				t.Sum()(2,0) = -1.0;
				return t;
			case 1:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) =  0.0;
				t.Sum()(2,0) = -1.0;
				return t;
			case 2:
				t.Sum()(0,0) =  0.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) = -1.0;
				return t;
			case 3:
				t.Sum()(0,0) =  0.0;
				t.Sum()(1,0) =  0.0;
				t.Sum()(2,0) =  1.0;
				return t;
			case 4:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) =  0.0;
				t.Sum()(2,0) =  1.0;
				return t;
			case 5:
				t.Sum()(0,0) =  0.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) =  1.0;
				return t;
			case 6:
				t.Mult()(0,0) =  0.5;
				t.Sum() (0,0) =  0.5;
				t.Sum() (2,0) = -1.0;
				return t;
			case 7:
				t.Mult()(0,0) = -0.5;
				t.Mult()(1,0) =  0.5;
				t.Sum() (0,0) =  0.5;
				t.Sum() (1,0) =  0.5;
				t.Sum() (2,0) = -1.0;
				return t;
			case 8:
				t.Mult()(1,0) = -0.5;
				t.Sum() (1,0) =  0.5;
				t.Sum() (2,0) = -1.0;
				return t;
			case 9:
				t.Mult()(2,0) =  1.0;
				return t;
			case 10:
				t.Mult()(2,0) =  1.0;
				t.Sum() (0,0) =  1.0;
				return t;
			case 11:
				t.Mult()(2,0) =  1.0;
				t.Sum() (1,0) =  1.0;
				return t;
			case 12:
				t.Mult()(0,0) =  0.5;
				t.Sum() (0,0) =  0.5;
				t.Sum() (2,0) =  1.0;
				return t;
			case 13:
				t.Mult()(0,0) = -0.5;
				t.Mult()(1,0) =  0.5;
				t.Sum() (0,0) =  0.5;
				t.Sum() (1,0) =  0.5;
				t.Sum() (2,0) =  1.0;
				return t;
			case 14:
				t.Mult()(1,0) = -0.5;
				t.Sum() (1,0) =  0.5;
				t.Sum() (2,0) =  1.0;
				return t;
			case 15:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Sum() (2,0) = -1.0;
				return t;
			case 16:
				t.Mult()(0,0) =  0.5;
				t.Mult()(2,1) =  1.0;
				t.Sum() (0,0) =  0.5;
				return t;
			case 17:
				t.Mult()(0,0) = -0.5;
				t.Mult()(1,0) =  0.5;
				t.Mult()(2,1) =  1.0;
				t.Sum() (0,0) =  0.5;
				t.Sum() (1,0) =  0.5;
				return t;
			case 18:
				t.Mult()(1,0) =  0.5;
				t.Mult()(2,1) =  1.0;
				t.Sum() (1,0) =  0.5;
				return t;
			case 19:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Sum() (2,0) =  1.0;
				return t;
			case 20:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,2) =  1.0;
				return t;
		}
		return TPZTransform(0,0);
	}
	
	
	TPZIntPoints * TPZPrism::CreateSideIntegrationRule(int side, int order){
		if(side<0 || side>20) {
			PZError << "TPZPrism::CreateSideIntegrationRule. bad side number.\n";
			return 0;
		}
		if(side<6)   return new TPZInt1Point(order);                   // sides 0 to 7 are vertices
		if(side<15)  return new TPZInt1d(order);                  // sides 7 to 14 are lines
		if(side==15 || side==19) return new TPZIntTriang(order);  // sides 15 and 19 are triangles
		if(side<20)  {
			return new TPZIntQuad(order,order);                   // sides 16 to 18 are quadrilaterals
		}
		if(side==20) {
			return new IntruleType(order,order);                  // integration of the element
		}
		return 0;
	}
	
	
	MElementType TPZPrism::Type()
	{
		return EPrisma;
	}
	
	MElementType TPZPrism::Type(int side)
	{
		switch(side) {
			case 0:
			case 1:
			case 2:
			case 3:
			case 4:
			case 5:
				return EPoint;
			case 6:
			case 7:
			case 8:
			case 9:
			case 10:
			case 11:
			case 12:
			case 13:
			case 14:
				return EOned;
			case 15:
				return ETriangle;
			case 16:
			case 17:
			case 18:
				return EQuadrilateral;
			case 19:
				return ETriangle;
			case 20:
				return EPrisma;
			default:
				return ENoType;
		}
	}
	
	int TPZPrism::NumSides() {
		return NSides;
	}
	
	
	int TPZPrism::NContainedSides(int side) {
		if(side<0)   return -1;
		if(side<6)   return 1;//cantos : 0 a 5
		if(side<15)  return 3;//arestas
		if(side==15 || side==19)  return 7;//faces : 15,19 , triangulares
		if(side<19) return 9;//faces : 16 a 18  quadrilaterais
		if(side==20) return 21;//centro : 20
		return -1;
	}
	
	int TPZPrism::ContainedSideLocId(int side, int node) {
		if(side<0 || side>20 || node < 0) return -1;
		if(side<6) {
			if(node==0) return side;
		} else
			if(side<9) {//6,7,8
				int s = side-6;//0,1,2
				if(!node) return s;//0,1,2
				if(node==1) return (s+1)%3;//1,2,0
				if(node==2) return side;//6,7,8
			} else
				if(side<12) {//9,10,11
					int s = side-9;//0,1,2
					if(!node) return s;//0,1,2,4
					if(node==1) return s+3;//3,4,5
					if(node==2) return side;//5,6,7
				} else
					if(side<15) {//12,13,14
						int s = side-9;//3,4,5
						if(!node) return s;//3,4,5
						if(node==1) return (s+1)%3+3;//4,5,3
						if(node==2) return side;//12,13,14
					} else
						if(side==15 || side==19) {
							int s = side-15;
							if(side==15 && node<7) return FaceConnectLocId[s][node];
							if(side==19 && node<7) return FaceConnectLocId[s][node];
						} else
							if(side<20) {//16,17,18
								int s = side-15;
								if(node<9) return FaceConnectLocId[s][node];
							} else
								if(side==20 && node<21){
									return node;
								}
		PZError << "TPZShapePrism::ContainedSideLocId called for node = "
		<< node << " and side = " << side << "\n";
		return -1;
	}
	
	bool TPZPrism::IsInParametricDomain(TPZVec<REAL> &pt, REAL tol){
		const REAL qsi = pt[0];
		const REAL eta = pt[1];
		const REAL zeta = pt[2];
		if( ( qsi <= 1. + tol ) && ( qsi >= 0. - tol ) &&
		   ( eta <= 1. + tol ) && ( eta >= 0. - tol ) &&
		   ( eta <= 1. - qsi + tol ) && 
		   ( zeta <= 1. + tol ) && ( zeta >= -1. - tol) ){
			return true;
		}
		else{
			return false;
		}  
		
	}//method
	
	/**
	 * Method which identifies the transformation based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */
	int TPZPrism::GetTransformId(TPZVec<int> &id)
	{
		LOGPZ_ERROR(logger,"Please implement me")
		return -1;
	}
	
	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	int TPZPrism::GetTransformId(int side, TPZVec<int> &id)
	{
		LOGPZ_ERROR(logger,"Please implement me")
		return -1;
	}
	
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
	void TPZPrism::GetSideHDivPermutation(int side, TPZVec<int> &id, TPZVec<int> &permgather)
	{
		LOGPZ_ERROR(logger,"Please implement me")
		int nel = permgather.NElements();
		int iel;
		for(iel=0; iel<nel; iel++)
			permgather[iel]=iel;
	}

}
