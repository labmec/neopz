/**
 * @file
 * @brief Contains the implementation of the TPZTetrahedron methods. 
 */

#include "tpztetrahedron.h"

#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"
#include "pzquad.h"
#include "pzeltype.h"
#include "tpztriangle.h"

#ifdef _AUTODIFF
#include "fad.h"
#endif

#include "pzcreateapproxspace.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.topology.pztetrahedron"));
#endif

using namespace std;

namespace pztopology {

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
		/*00*/{.0,.0,.0},/*01*/{1.,.0,.0},/*02*/{0.,1.,.0},/*03*/{.0,0.,1.0},/*04*/{.5,.0,.0},
		/*05*/{.5,.5,.0},/*06*/{0.,.5,.0},/*07*/{0.,0.,.5},/*08*/{.5,0.,0.5},/*09*/{.0,.5,.5},
		/*10*/{1./3.,1./3., 0.  }  ,/*11*/{1./3., .0  ,1./3.},
		/*12*/{1./3.,1./3.,1./3.}  ,/*13*/{ 0.  ,1./3.,1./3.},/*14*/{1./4.,1./4.,1./4.} };
    
    static REAL bTetra[45][3] = // direcao perpendicular ao lado
    {
        {0,0,-1}, {1,0,-1}, {0,1,-1}, {0,0,-1}, {0.5,0.5,-1}, {0,0,-1}, {0,0,-1},// face 0
        {0,-1,0}, {1,-1,0},  {0,-1,1}, {0,-1,0}, {0.5,-1,0.5}, {0,-1,0}, {0,-1,0},// face 1
        {1,0,0} , {0,1,0} ,  {0,0,1}, {1,1,0}, {0,0.5,0.5}, {0.5,0,0.5}, {1,1,1} ,// face 2
        {-1,0,0}, {-1,1,0}, {-1,0,1}, {-1,0,0}, {-1,0.5,0.5}, {-1,0,0} , {-1,0,0},// face 3
        //interior
        //aresta
        {1,0,0},{-1,1,0},{0,-1,0},  {0,0,1},  {-1,0,1},  {0,-1,1},
        //faces
        {-1,0,0}, {0,1,0},// face 0
        {1,0,0}, {0,0,1},// face 1
        {1,0,-1}, {-1,2,-1},// face 2
        {0,0,1}, {0,1,0} ,// face 3
        //interior
        {1,0,0} ,
        {0,1,0} ,
        {0,0,1}
    };
    static REAL t1Tetra[45][3] =
    {
        {-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},//face 0
        {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0}, //face 1
        {1/sqrt(2),0,-1/sqrt(2)},{1/sqrt(2),0,-1/sqrt(2)},{1/sqrt(2),0,-1/sqrt(2)},{1/sqrt(2),0,-1/sqrt(2)},{1/sqrt(2),0,-1/sqrt(2)},{1/sqrt(2),0,-1/sqrt(2)},{1/sqrt(2),0,-1/sqrt(2)},//face 2
        {0,0,1} ,{0,0,1} ,{0,0,1} ,{0,0,1} ,{0,0,1} ,{0,0,1} ,{0,0,1} ,//face 3
        //interior
        //aresta
        {0,-1,0},{1,1,0},{0,0,-1},  {0,-1,0},  {0,-1,0},  {1,1,1},
        //faces
        {0,1,0}, {1,0,0},// face 0
        {0,0,1}, {-1,0,0},// face 1
        {-1,2,-1}, {-1,0,1},// face 2
        {0,1,0}, {0,0,-1} ,// face 3
        //interior
        {0,1,0} ,
        {0,0,1} ,
        {1,0,0}
        
    };
    static REAL t2Tetra[45][3] =
    {
        {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, // face 0
        {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1},// face 1
        {-1/sqrt(6),2/sqrt(6),-1/sqrt(6)},{-1/sqrt(6),2/sqrt(6),-1/sqrt(6)},{-1/sqrt(6),2/sqrt(6),-1/sqrt(6)},{-1/sqrt(6),2/sqrt(6),-1/sqrt(6)},{-1/sqrt(6),2/sqrt(6),-1/sqrt(6)},{-1/sqrt(6),2/sqrt(6),-1/sqrt(6)},{-1/sqrt(6),2/sqrt(6),-1/sqrt(6)},// face 2
        {0,1,0},{0,1,0},{0,1,0},{0,1,0},{0,1,0},{0,1,0},{0,1,0},// face 3
        //interior
        //aresta
        {0,0,-1},{0,0,-1},{-1,0,0},  {-1,0,0},  {1,1,1},  {-1,0,0},
        //faces
        {0,0,-1}, {0,0,-1},// face 0
        {0,-1,0}, {0,-1,0},// face 1
        {1,1,1}, {1,1,1},// face 2
        {-1,0,0}, {-1,0,0} ,// face 3
        //interior
        {0,0,1} ,
        {1,0,0} ,
        {0,1,0}
    };

    static int vectorsideorderTe [45] =
    {
        0,1,2,4,5,6,10, //face 0
        0,1,3,4,8,7,11,//face 1
        1,2,3,5,9,8,12,//face 2
        0,2,3,6,9,7,13,//face 3
        4,5,6,7,
        8,9,
        10,10,//tg face 0
        11,11,//tg face 1
        12,12,//tg face 2
        13,13,//tg face 3
        14,14,14
    };
    
    static int bilinearounao [45] =
    {
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1
    }; //P*k Pk

//    static int bilinearounao [45] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; //Pk Pk-1

    static int direcaoksioueta [45] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,2};
    
	int TPZTetrahedron::NBilinearSides()
    {
        DebugStop();
        return 0;
    }
    
	void TPZTetrahedron::LowerDimensionSides(int side,TPZStack<int> &smallsides)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		int is;
		for(is=0; is<nsidecon-1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	void TPZTetrahedron::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
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
	//Tentando criar o metodo
	int TPZTetrahedron::NumSides(int dimension) {		if(dimension<0 || dimension> 3) {
		PZError << "TPZTetrahedron::NumSides. Bad parameter i.\n";
		return 0;
	}
		if(dimension==0) return 4;
		if(dimension==1) return 6;
		if(dimension==2) return 4;
		if(dimension==3) return 1;
		return -1;
	}
	
	int TPZTetrahedron::SideNodeLocId(int side, int node)
	{
		if(side<4 && node == 0) return side;
		if(side>=4 && side < 10 && node <2) return SideNodes[side-4][node];
		if(side >= 10 && side < 14 && node <3) return FaceNodes[side-10][node];
		if(side ==14 && node < 4) return node;
		PZError << "TPZTetrahedron::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
        DebugStop();
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
	
	TPZTransform<> TPZTetrahedron::SideToSideTransform(int sidefrom, int sideto)
	{
		if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
			PZError << "TPZTetrahedron::HigherDimensionSides sidefrom "<< sidefrom << 
			' ' << sideto << endl;
			return TPZTransform<>(0);
		}
		if(sidefrom == sideto) {
			return TPZTransform<>(sidedimension[sidefrom]);
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
				TPZTransform<> trans(dto,dfr);
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
		return TPZTransform<>(0);
	}
	
	TPZTransform<> TPZTetrahedron::TransformElementToSide(int side){
		
		if(side<0 || side>14){
			PZError << "TPZTetrahedron::TransformElementToSide called with side error\n";
			return TPZTransform<>(0,0);
		}
		
		TPZTransform<> t(sidedimension[side],3);
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
		return TPZTransform<>(0,0);
	}
	
	TPZTransform<> TPZTetrahedron::TransformSideToElement(int side){
		
		if(side<0 || side>14){
			PZError << "TPZTetrahedron::TransformSideToElement side out range\n";
			return TPZTransform<>(0,0);
		}
		TPZTransform<> t(3,sidedimension[side]);
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
		return TPZTransform<>(0,0);
	}
	
	TPZIntPoints * TPZTetrahedron::CreateSideIntegrationRule(int side, int order) {
		if(side<0 || side>15) {
			PZError << "TPZTetrahedron::CreateSideIntegrationRule. bad side number.\n";
			return 0;
		}
		if(side<4)   return new TPZInt1Point(order);    // sides 0 to 3 are points (vertices)
		if(side<10)  return new TPZInt1d(order);   // sides 4 to 9 are lines
		if(side<14)  {                             // sides 10 to 13 are triangles
			return new TPZIntTriang(order);
		}
		if(side==14) {                            // integration of the element
			return new IntruleType(order);
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
	
	
	int TPZTetrahedron::NumSides() {
		return NSides;
	}
	
	
	int TPZTetrahedron::NContainedSides(int side) {
		if(side<0)   return -1;
		if(side<4)   return 1;//cantos : 0 a 3
		if(side<10)  return 3;//lados : 4 a 9
		if(side<14)	 return 7;//faces : 10 a 13
		if(side==14) return 15;//centro : 14
		return -1;
	}
	
	int TPZTetrahedron::ContainedSideLocId(int side, int node) {
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
		PZError << "TPZShapeTetra::ContainedSideLocId called for node = "
		<< node << " and side = " << side << "\n";
		return -1;
	}
	
	bool TPZTetrahedron::IsInParametricDomain(TPZVec<REAL> &pt, REAL tol){
		const REAL qsi = pt[0];
		const REAL eta = pt[1];
		const REAL zeta = pt[2];
		if( (qsi < 0. - tol) || (qsi > 1. + tol) ||
		   (eta < 0. - tol) || (eta > 1. + tol) || 
		   (zeta < 0. -tol) || (zeta > 1. +tol) || 
		   (qsi+eta+zeta > 1.+tol) ) {
			return false;
		}
		else{
			return true;
		}
		
	}//method
    
    
    /** @brief Generates a random point in the master domain */
    void TPZTetrahedron::RandomPoint(TPZVec<REAL> &pt)
    {
        REAL val = (REAL) rand() / (RAND_MAX);
        pt[0] = val;
        val = (1.-pt[0]) * (REAL) rand() / (RAND_MAX);
        pt[1] = val;
        val = (1.-pt[0]-pt[1]) * (REAL) rand() / (RAND_MAX);
        pt[2] = val;
    }
    
    template<class T>
    bool TPZTetrahedron::MapToSide(int side, TPZVec<T> &InternalPar, TPZVec<T> &SidePar, TPZFMatrix<T> &JacToSide) {
		
		T zero = 1.E-5;
		
		T qsi = InternalPar[0], eta = InternalPar[1], zeta = InternalPar[2];
		bool regularmap = true;
		switch(side)
		{
            case 0:
            case 1:
            case 2:
            case 3:
            {
                SidePar.Resize(0); JacToSide.Resize(0,0);
                break;
            }
			case 4://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs((T)(eta + zeta - 1.)) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = -(2.*qsi + eta + zeta -1.)/(eta + zeta -1.);
                    JacToSide(0,0) = -2./(eta + zeta - 1.); JacToSide(0,1) = 2.*qsi/((eta + zeta - 1.)*(eta + zeta - 1.)); JacToSide(0,2) = 2.*qsi/((eta + zeta - 1.)*(eta + zeta - 1.));
				}
				break;
				
			case 5://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs((T)(qsi + eta)) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = (eta - qsi)/(eta + qsi);
                    JacToSide(0,0) = -2.*eta/((qsi + eta)*(qsi + eta)); JacToSide(0,1) = 2.*qsi/((qsi + eta)*(qsi + eta)); JacToSide(0,2) = 0.;
				}
				break;
				
			case 6://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs((T)(qsi + zeta - 1.)) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = (qsi + 2.*eta + zeta -1.)/(qsi + zeta -1.);
                    JacToSide(0,0) = -2.*eta/((qsi + zeta - 1.)*(qsi + zeta - 1.)); JacToSide(0,1) = 2./(qsi + zeta - 1.); JacToSide(0,2) = -2.*eta/((qsi + zeta - 1.)*(qsi + zeta - 1.));
				}
				break;
				
			case 7://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs((T)(qsi + eta - 1.)) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = -(qsi + eta + 2.*zeta -1.)/(qsi + eta -1.);
                    JacToSide(0,0) = 2.*zeta/((qsi + eta - 1.)*(qsi + eta - 1.)); JacToSide(0,1) = 2.*zeta/((qsi + eta - 1.)*(qsi + eta - 1.)); JacToSide(0,2) = -2./(qsi + eta -1.);
				}
				break;
				
			case 8://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs((T)(qsi + zeta)) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = (zeta - qsi)/(zeta + qsi);
                    JacToSide(0,0) = -2.*zeta/(qsi + zeta); JacToSide(0,1) = 0.; JacToSide(0,2) = 2.*qsi/((qsi + zeta)*(qsi + zeta));
				}
				break;
				
			case 9://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs((T)(eta + zeta)) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = (zeta - eta)/(zeta + eta);
                    JacToSide(0,0) = -2./(eta + zeta - 1.); JacToSide(0,1) = 2.*qsi/((eta + zeta - 1.)*(eta + zeta - 1.)); JacToSide(0,2) = 2.*qsi/((eta + zeta - 1.)*(eta + zeta - 1.));
				}
				break;
				
			case 10://2D
				SidePar.Resize(2); JacToSide.Resize(2,3);
				if(fabs((T)(zeta - 1.)) < zero)
				{
                    SidePar[0] = 1./3.; SidePar[1] = 1./3.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = qsi/(1.-zeta); SidePar[1] = eta/(1.-zeta);
                    JacToSide(0,0) = 1/(1.-zeta); JacToSide(0,1) = 0.; JacToSide(0,2) = qsi/((zeta-1.)*(zeta-1.));
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 1/(1.-zeta); JacToSide(1,2) = eta/((zeta-1.)*(zeta-1.));
				}
				break;
				
			case 11://2D
				SidePar.Resize(2); JacToSide.Resize(2,3);
				if(fabs((T)(eta - 1.)) < zero)
				{
                    SidePar[0] = 1./3.; SidePar[1] = 1./3.;
                    JacToSide(0,0) = 1./(1.-eta); JacToSide(0,1) = qsi/((eta-1.)*(eta-1.)); JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = zeta/((eta-1.)*(eta-1.)); JacToSide(1,2) = 1./(1.-eta);
					regularmap = false;
				}
				else
				{
                    SidePar[0] = qsi/(1.-eta); SidePar[1] = zeta/(1.-eta);
                    JacToSide(0,0) = 1/(1.-zeta); JacToSide(0,1) = 1/(1.-zeta); JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 1/(1.-zeta); JacToSide(1,1) = 1/(1.-zeta); JacToSide(1,2) = 0.;
				}
				break;
				
			case 12://2D
				SidePar.Resize(2); JacToSide.Resize(2,3);
				if(fabs((T)(qsi+eta+zeta)) < zero)
				{
                    SidePar[0] = 1./6.; SidePar[1] = 1./3.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = eta/(qsi+eta+zeta); SidePar[1] = zeta/(qsi+eta+zeta);
                    JacToSide(0,0) = -eta/((qsi+eta+zeta)*(qsi+eta+zeta)); JacToSide(0,1) = (qsi+zeta)/((qsi+eta+zeta)*(qsi+eta+zeta)); JacToSide(0,2) = -eta/((qsi+eta+zeta)*(qsi+eta+zeta));
                    JacToSide(1,0) = -zeta/((qsi+eta+zeta)*(qsi+eta+zeta)); JacToSide(1,1) = -zeta/((qsi+eta+zeta)*(qsi+eta+zeta)); JacToSide(1,2) = (qsi+eta)/((qsi+eta+zeta)*(qsi+eta+zeta));
				}
				break;
				
			case 13://2D
				SidePar.Resize(2); JacToSide.Resize(2,3);
				if(fabs((T)(qsi - 1.)) < zero)
				{
                    SidePar[0] = 1./3.; SidePar[1] = 1./3.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = eta/(1.-qsi); SidePar[1] = zeta/(1.-qsi);
                    JacToSide(0,0) = eta/((qsi-1.)*(qsi-1.)); JacToSide(0,1) = 1/(1.-qsi); JacToSide(0,2) = 0.;
                    JacToSide(1,0) = zeta/((qsi-1.)*(qsi-1.)); JacToSide(1,1) = 0.; JacToSide(1,2) = 1/(1.-qsi);
				}
				break;
            case 14: //Interno para interno ? aquidouglas aqui
                SidePar = InternalPar;
                JacToSide.Resize(3, 3);
                JacToSide.Identity();
                break;
		}
		if(side > 14)
		{
			cout << "Cant compute MapToSide method in TPZGeoTetrahedra class!\nParameter (SIDE) must be between 4 and 13!\nMethod Aborted!\n";
			DebugStop();
		}
		return regularmap;
	}
    
    void TPZTetrahedron::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
    {
        if(node > NCornerNodes)
        {
            DebugStop();
        }
        nodeCoord.Resize(Dimension, 0.);
        switch (node) {
            case (0):
            {
                nodeCoord[0] = 0.;
                nodeCoord[1] = 0.;
                nodeCoord[2] = 0.;
                break;
            }
            case (1):
            {
                nodeCoord[0] = 1.;
                nodeCoord[1] = 0.;
                nodeCoord[2] = 0.;
                break;
            }
            case (2):
            {
                nodeCoord[0] = 0.;
                nodeCoord[1] = 1.;
                nodeCoord[2] = 0.;
                break;
            }
            case (3):
            {
                nodeCoord[0] = 0.;
                nodeCoord[1] = 0.;
                nodeCoord[2] = 1.;
                break;
            }
            default:
            {
                DebugStop();
                break;
            }
        }
    }
	
	/**
	 * Method which identifies the transformation based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */
	int TPZTetrahedron::GetTransformId(TPZVec<int64_t> &id)
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
	int TPZTetrahedron::GetTransformId(int side, TPZVec<int64_t> &id)
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
	void TPZTetrahedron::GetSideHDivPermutation(int transformationid, TPZVec<int> &permgather)
	{
	// Not complete
        DebugStop();
        
//#ifdef PZDEBUG
//        if (SideDimension(side) != 2) {
//            DebugStop();
//        }
//#endif
//        permgather.Resize(7);
//        TPZManVector<int64_t,7> locids(3);
//        for (int in=0; in<3; in++) {
//            locids[in] = id[SideNodeLocId(side, in)];
//        }
//        int transformid = pztopology::TPZTriangle::GetTransformId(locids);
//        pztopology::TPZTriangle::GetHDivGatherPermute(transformid,permgather);
    }
    
    
    void computedirectionsT3(int inicio, int fim, TPZFMatrix<REAL> &bvec, TPZFMatrix<REAL> &t1vec,
                           TPZFMatrix<REAL> &t2vec, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions);
    
    void computedirectionsT3(int inicio, int fim, TPZFMatrix<REAL> &bvec, TPZFMatrix<REAL> &t1vec,
                           TPZFMatrix<REAL> &t2vec, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions)
    {
        REAL detgrad = 0.0;
        TPZVec<REAL> u(3);
        TPZVec<REAL> v(3);
        TPZVec<REAL> uxv(3);// result
        int cont = 0;
        
        for (int ivet=inicio; ivet<=fim; ivet++)
        {
            for (int ilin=0; ilin<3; ilin++)
            {
                u[ilin] = t1vec(ilin,ivet);
                v[ilin] = t2vec(ilin,ivet);
            }
            TPZVec<REAL> e2(3);
            detgrad = 0.0;
            REAL normaX0xX1 = 0.0;
            //TPZNumeric::ProdVetorial(u,v,e2);
            e2[0] = u[1]*v[2]-u[2]*v[1];
            e2[1] = -(u[0]*v[2]-v[0]*u[2]);
            e2[2] = u[0]*v[1]-v[0]*u[1];
            
            // calc do v gradx*b
            TPZManVector<REAL,3> dxt1(3,0.),dxt2(3,0.),dxt3(3,0.),Vvec(3,0.);
            REAL be2 = 0.0, ne2 = 0.0;
            for(int i=0;i<3;i++)
            {
                ne2 += e2[i]*e2[i];
            }
            ne2 = sqrt(fabs(ne2));
            for (int il=0; il<3; il++)
            {
                for (int i = 0 ; i<3; i++)
                {
                    dxt1[il] += gradx(il,i) * t1vec(i,ivet);
                    dxt2[il] += gradx(il,i) * t2vec(i,ivet);
                    dxt3[il] += gradx(il,i) * e2[i]/ne2;
                    Vvec[il] += gradx(il,i) * bvec(i,ivet);
                }
                be2 += bvec(il,ivet)*e2[il]/ne2;
            }
            TPZManVector<REAL,3> normal(3,0.);
            //TPZNumeric::ProdVetorial(dxt1,dxt2,normal);
            normal[0] = dxt1[1]*dxt2[2]-dxt1[2]*dxt2[1];
            normal[1] = -(dxt1[0]*dxt2[2]-dxt2[0]*dxt1[2]);
            normal[2] = dxt1[0]*dxt2[1]-dxt2[0]*dxt1[1];
            
            for (int pos=0; pos<3; pos++)
            {
                detgrad += normal[pos]*dxt3[pos];//uxv[pos]*gradx.GetVal(pos, 2);
                normaX0xX1 += normal[pos]*normal[pos]; //uxv[pos]*uxv[pos];
            }
            TPZFMatrix<REAL> Wvec(3,1);
            detgrad = fabs(detgrad);
            normaX0xX1 = sqrt(normaX0xX1);
            
            for (int il=0; il<3; il++)
            {
                Wvec(il,0) = Vvec[il]*normaX0xX1/(detgrad*be2);
                directions(il,cont) = Wvec(il,0);
            }
            cont++;
        }
        
        
    }

    
    void TPZTetrahedron::ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors)
    {
        if(gradx.Cols()!=3)
        { std::cout << "Gradient dimensions are not compatible with this topology" << std::endl;
            DebugStop();
        }
        TPZFMatrix<REAL> bvec(3, 45);
        int numvec = bvec.Cols();
        TPZFMatrix<REAL> t1vec(3, numvec);
        TPZFMatrix<REAL> t2vec(3,numvec);
        
        directions.Redim(3, numvec);
        
        for (int lin = 0; lin<numvec ; lin++)
        {
            for(int col = 0;col<3;col++)
            {
                bvec.PutVal(col,  lin, bTetra[lin][col]);
                t1vec.PutVal(col, lin, t1Tetra[lin][col]);
                t2vec.PutVal(col, lin, t2Tetra[lin][col]);
            }
        }
        
        // calcula os vetores
        switch (side) {
            case 10:
            {
                directions.Resize(3, 7);
                sidevectors.Resize(7);
                int inicio = 0, fim = 6;
                computedirectionsT3( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderTe[ip+inicio];
                }
            }
                break;
            case 11:
            {
                directions.Resize(3, 7);
                sidevectors.Resize(7);
                int inicio = 7, fim = 13;
                computedirectionsT3( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderTe[ip+inicio];
                }
            }
                break;
            case 12:
            {
                directions.Resize(3, 7);
                sidevectors.Resize(7);
                int inicio = 14, fim = 20;
                computedirectionsT3( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderTe[ip+inicio];
                }
            }
                break;
            case 13:
            {
                directions.Resize(3, 7);
                sidevectors.Resize(7);
                int inicio = 21, fim = 27;
                computedirectionsT3( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderTe[ip+inicio];
                }
            }
                break;
            case 14:
            {
                directions.Resize(3, 17);
                sidevectors.Resize(17);
                int inicio = 28, fim =  44;
                computedirectionsT3( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderTe[ip+inicio];
                }
            }
                break;
                
            default:
                break;
        }
                
                
	}
    
    void TPZTetrahedron::ComputeDirections(TPZFMatrix<REAL> &gradx, REAL detjac, TPZFMatrix<REAL> &directions)
    {
        REAL detgrad = gradx(0,0)*gradx(1,1)*gradx(2,2) + gradx(0,1)*gradx(1,2)*gradx(2,0) + gradx(0,2)*gradx(1,0)*gradx(2,1) - gradx(0,2)*gradx(1,1)*gradx(2,0) - gradx(0,0)*gradx(1,2)*gradx(2,1) - gradx(0,1)*gradx(1,0)*gradx(2,2);
        detgrad = fabs(detgrad);
        
        TPZManVector<REAL,3> v1(3),v2(3),v3(3),v1v2(3),v3v1(3),v2v3(3),vdiagxy(3),vi(3),vivdiagxy(3);
        
        for (int i=0; i<3; i++) {
            v1[i] = gradx(i,0);
            v2[i] = gradx(i,1);
            v3[i] = gradx(i,2);
            vdiagxy[i] = (gradx(i,0)-gradx(i,1));
            vi[i] = gradx(i,2)-gradx(i,0);//gradx(i,2)-0.5*(gradx(i,0)+gradx(i,1));
        }
        
        
        /**
         * @file
         * @brief Computing mapped vector with scaling factor equal 1.0.
         * using contravariant piola mapping.
         */
        
        REAL Nv1v2 = 1.0;
        REAL Nv2v3 = 1.0;
        REAL Nv3v1 = 1.0;
        REAL Nvivdiagb = 1.0;

        if (HDivPiola == 1) {
            // the above constants are wrong
            for (int i=0; i<3; i++) {
                v1[i] /= detgrad;
                v2[i] /= detgrad;
                v3[i] /= detgrad;
            }
            for (int i=0; i<3; i++)
            {
                
                //face 0
                directions(i,0) = -v3[i];
                directions(i,1) = (v1[i]-v3[i]);
                directions(i,2) = (v2[i]-v3[i]);
                directions(i,3) = (directions(i,0)+directions(i,1))/2.;
                directions(i,4) = (directions(i,1)+directions(i,2))/2.;
                directions(i,5) = (directions(i,0)+directions(i,2))/2.;
                directions(i,6) = (directions(i,3)+directions(i,4)+directions(i,5))/3.;
                //face 1
                directions(i,7) = -v2[i];
                directions(i,8) = (v1[i]-v2[i]);
                directions(i,9) = (v3[i]-v2[i]);
                directions(i,10) = (directions(i,7)+directions(i,8))/2.;
                directions(i,11) = (directions(i,8)+directions(i,9))/2.;
                directions(i,12) = (directions(i,7)+directions(i,9))/2.;
                directions(i,13) = (directions(i,10)+directions(i,11)+directions(i,12))/3.;
                //face 2 face diagonal
                
                directions(i,14) = v1[i];
                directions(i,15) = v2[i];
                directions(i,16) = v3[i];
                directions(i,17) = (directions(i,14)+directions(i,15))/2.;
                directions(i,18) = (directions(i,15)+directions(i,16))/2.;
                directions(i,19) = (directions(i,14)+directions(i,16))/2.;
                directions(i,20) = (directions(i,17)+directions(i,18)+directions(i,19))/3.;
                //face 3
                directions(i,21) = -v1[i];
                directions(i,22) = (v2[i]-v1[i]);
                directions(i,23) = (v3[i]-v1[i]);
                directions(i,24) = (directions(i,21)+directions(i,22))/2.;
                directions(i,25) = (directions(i,22)+directions(i,23))/2.;
                directions(i,26) = (directions(i,21)+directions(i,23))/2.;
                directions(i,27) = (directions(i,24)+directions(i,25)+directions(i,26))/3.;
                
                //arestas
                directions(i,28) = v1[i];
                directions(i,29) = (v2[i]-v1[i]);
                directions(i,30) = -v2[i];
                directions(i,31) = v3[i];
                directions(i,32) = (v3[i]-v1[i]);
                directions(i,33) = (v3[i]-v2[i]);
                
                //faces
                directions(i,34) = v1[i];
                directions(i,35) = v2[i];
                directions(i,36) = v1[i];
                directions(i,37) = v3[i];
                directions(i,38) = (v2[i]-v1[i]);
                directions(i,39) = (v3[i]-v1[i]);//v3[i]-0.5*(v1[i]+v2[i]);//
                directions(i,40) = v2[i];
                directions(i,41) = v3[i];
                
                directions(i,42) = v1[i];
                directions(i,43) = v2[i];
                directions(i,44) = v3[i];
                
            }        

        }
        else
        {
            TPZNumeric::ProdVetorial(v1,v2,v1v2);
            TPZNumeric::ProdVetorial(v2,v3,v2v3);
            TPZNumeric::ProdVetorial(v3,v1,v3v1);
            TPZNumeric::ProdVetorial(vi,vdiagxy,vivdiagxy);
            
            Nv1v2 = TPZNumeric::Norma(v1v2);
            Nv2v3 = TPZNumeric::Norma(v2v3);
            Nv3v1 = TPZNumeric::Norma(v3v1);
            Nvivdiagb = TPZNumeric::Norma(vivdiagxy);


        
            for (int i=0; i<3; i++) {
                v1[i] /= detgrad;
                v2[i] /= detgrad;
                v3[i] /= detgrad;
            }
            for (int i=0; i<3; i++)
            {
                
                //face 0
                directions(i,0) = -v3[i]*Nv1v2;
                directions(i,1) = (v1[i]-v3[i])*Nv1v2;                
                directions(i,2) = (v2[i]-v3[i])*Nv1v2;                
                directions(i,3) = (directions(i,0)+directions(i,1))/2.;
                directions(i,4) = (directions(i,1)+directions(i,2))/2.;
                directions(i,5) = (directions(i,0)+directions(i,2))/2.;
                directions(i,6) = (directions(i,3)+directions(i,4)+directions(i,5))/3.;
                //face 1
                directions(i,7) = -v2[i]*Nv3v1;
                directions(i,8) = (v1[i]-v2[i])*Nv3v1;                
                directions(i,9) = (v3[i]-v2[i])*Nv3v1;                
                directions(i,10) = (directions(i,7)+directions(i,8))/2.;
                directions(i,11) = (directions(i,8)+directions(i,9))/2.;     
                directions(i,12) = (directions(i,7)+directions(i,9))/2.;    
                directions(i,13) = (directions(i,10)+directions(i,11)+directions(i,12))/3.;
                //face 2
                
                directions(i,14) = v1[i]*Nvivdiagb; 
                directions(i,15) = v2[i]*Nvivdiagb; 
                directions(i,16) = v3[i]*Nvivdiagb; 
                directions(i,17) = (directions(i,14)+directions(i,15))/2.;              
                directions(i,18) = (directions(i,15)+directions(i,16))/2.;               
                directions(i,19) = (directions(i,14)+directions(i,16))/2.;              
                directions(i,20) = (directions(i,17)+directions(i,18)+directions(i,19))/3.; 
                //face 3
                directions(i,21) = -v1[i]*Nv2v3;
                directions(i,22) = (v2[i]-v1[i])*Nv2v3;               
                directions(i,23) = (v3[i]-v1[i])*Nv2v3;
                directions(i,24) = (directions(i,21)+directions(i,22))/2.;
                directions(i,25) = (directions(i,22)+directions(i,23))/2.;
                directions(i,26) = (directions(i,21)+directions(i,23))/2.;
                directions(i,27) = (directions(i,24)+directions(i,25)+directions(i,26))/3.;
                
                //arestas
                directions(i,28) = v1[i];
                directions(i,29) = (v2[i]-v1[i]);                           
                directions(i,30) = -v2[i];
                directions(i,31) = v3[i];
                directions(i,32) = (v3[i]-v1[i]);               
                directions(i,33) = (v3[i]-v2[i]);               
                
                //faces
                directions(i,34) = v1[i];
                directions(i,35) = v2[i];
                directions(i,36) = v1[i];
                directions(i,37) = v3[i];
                directions(i,38) = (v2[i]-v1[i]);                
                directions(i,39) = (v3[i]-v1[i]);//v3[i]-0.5*(v1[i]+v2[i]);//
                directions(i,40) = v2[i];
                directions(i,41) = v3[i];
                
                directions(i,42) = v1[i];
                directions(i,43) = v2[i];
                directions(i,44) = v3[i];
                
            }        
        }
    }

    
    void TPZTetrahedron::GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao)
    {
        int nsides = NumSides()*3;
        
        sides.Resize(nsides);
        dir.Resize(nsides);
        bilounao.Resize(nsides);
        
        for (int is = 0; is<nsides; is++)
        {
            sides[is] = vectorsideorderTe[is];
            dir[is] = direcaoksioueta[is];
            bilounao[is] = bilinearounao[is];
        }
    }
    
    void TPZTetrahedron::GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao, TPZVec<int> &sidevectors)
    {
        int nsides = NumSides()*3;
        
        sides.Resize(nsides);
        dir.Resize(nsides);
        bilounao.Resize(nsides);
        
        for (int is = 0; is<nsides; is++)
        {
            sides[is] = vectorsideorderTe[is];
            dir[is] = direcaoksioueta[is];
            bilounao[is] = bilinearounao[is];
        }
        for (int i=0; i<Dimension*NumSides(); i++) {
            sidevectors[i] = vectorsideorderTe[i];
        }
    }
    
    
    int TPZTetrahedron::ClassId() const{
        return Hash("TPZTetrahedron");
    }

    void TPZTetrahedron::Read(TPZStream& buf, void* context) {

    }

    void TPZTetrahedron::Write(TPZStream& buf, int withclassid) const {

    }


}

template
bool pztopology::TPZTetrahedron::MapToSide<REAL>(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide);

#ifdef _AUTODIFF
template
bool pztopology::TPZTetrahedron::MapToSide<Fad<REAL> >(int side, TPZVec<Fad<REAL> > &InternalPar, TPZVec<Fad<REAL> > &SidePar, TPZFMatrix<Fad<REAL> > &JacToSide);
#endif
