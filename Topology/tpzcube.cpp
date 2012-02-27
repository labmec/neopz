/**
 * @file
 * @brief Contains the implementation of the TPZCube methods. 
 */

#include "tpzcube.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"
#include "pzquad.h"
#include "pzeltype.h"
#include "tpzquadrilateral.h"

#include "pzcreateapproxspace.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.topology.pzcube"));
#endif
using namespace std;

namespace pztopology {

	/**
	 * @brief For each face (quadrilateral sides) was enumerated in sequence the sides contained in the closure of them.
	 * For example: First face (side 20) in its closure contains side 0, 1, 2, 3 (vertices), 8, 9, 10, 11, (edges) and it self 
	 */
	static int FaceConnectLocId[6][9] = { {0,1,2,3,8,9,10,11,20},{0,1,5,4,8,13,16,12,21},
		{1,2,6,5,9,14,17,13,22},{3,2,6,7,10,14,18,15,23},//{2,3,7,6,10,15,18,14,23}
		{0,3,7,4,11,15,19,12,24},{4,5,6,7,16,17,18,19,25} };
	
	/** @brief For each face was enumerated the pontoal sides (vertices) */
	int TPZCube::FaceNodes[6][4]  = { {0,1,2,3},{0,1,5,4},{1,2,6,5},{3,2,6,7},{0,3,7,4},{4,5,6,7} };
	
	/** @brief For each edge was enumerated the pontoal sides (vertices) */
	int TPZCube::SideNodes[12][2]  = { {0,1},{1,2},{2,3},{3,0},{0,4},{1,5},
		{2,6},{3,7},{4,5},{5,6},{6,7},{7,4} };
	
	/** @brief For each face was enumerated the vertice sides on its main diagonal */
	int TPZCube::ShapeFaceId[6][2] = { {0,2},{0,5},{1,6},{3,6},{0,7},{4,6} };
	
	/** @brief Vector of the dimension for each side */
	static int sidedimension[27] = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3};
	
	/** @brief Vector with the number of vertices contained in the closure of the side */
	static int nsidenodes[27] = {1,1,1,1,1,1,1,1,
		2,2,2,2,2,2,2,2,2,2,2,2,
		4,4,4,4,4,4,
		8};
	

	/** @brief For each side was stored the sides connected with it but of the higher dimension */ 
	static int highsides[27][7] = {
		{8,11,12,20,21,24,26},
		{8,9,13,20,21,22,26},
		{9,10,14,20,22,23,26},
		{10,11,15,20,23,24,26},
		{12,16,19,21,24,25,26},
		{13,16,17,21,22,25,26},
		{14,17,18,22,23,25,26},
		{15,18,19,23,24,25,26},
		{20,21,26},
		{20,22,26},
		{20,23,26},
		{20,24,26},
		{21,24,26},
		{21,22,26},
		{22,23,26},
		{23,24,26},
		{21,25,26},
		{22,25,26},
		{23,25,26},
		{24,25,26},
		{26},
		{26},
		{26},
		{26},
		{26},
		{26},
		{-999}
	};

	/**
	 * @brief For each side was stored the number of sides connected of the higher dimension than it self 
	 * For example: First side (side 0 - vertice) was connected with the sides 8, 11, 12 (edges) 20, 21, 24 (faces) and 26 (the hexahedra). At total 7 sides */
	static int nhighdimsides[27] = {7,7,7,7,7,7,7,7,3,3,3,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,0};

	/** @brief The transformations for each side over neighboard side with higher dimension */
	static REAL sidetosidetransforms[27][7][4][3] = {
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,1}}
		},
		{
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{1,0,0},{-99,-99,-99},{-99,-99,-99},{0,-1,-1}}
		},
		{
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{0,1,0},{-99,-99,-99},{-99,-99,-99},{1,0,-1}}
		},
		{
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{-1,0,0},{-99,-99,-99},{-99,-99,-99},{0,1,-1}}
		},
		{
			{{0,-1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{0,-1,0},{-99,-99,-99},{-99,-99,-99},{-1,0,-1}}
		},
		{
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{-1,-1,0}}
		},
		{
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{1,-1,0}}
		},
		{
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{1,1,0}}
		},
		{
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{-1,1,0}}
		},
		{
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{1,0,0},{-99,-99,-99},{-99,-99,-99},{0,-1,1}}
		},
		{
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,1,0},{-99,-99,-99},{-99,-99,-99},{1,0,1}}
		},
		{
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-1,0,0},{-99,-99,-99},{-99,-99,-99},{0,1,1}}
		},
		{
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{0,-1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,-1,0},{-99,-99,-99},{-99,-99,-99},{-1,0,1}}
		},
		{
			{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,-1}}
		},
		{
			{{1,0,0},{0,0,1},{-99,-99,-99},{0,-1,0}}
		},
		{
			{{0,1,0},{0,0,1},{-99,-99,-99},{1,0,0}}
		},
		{
			{{1,0,0},{0,0,1},{-99,-99,-99},{0,1,0}}
		},
		{
			{{0,1,0},{0,0,1},{-99,-99,-99},{-1,0,0}}
		},
		{
			{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
		}
	};
	
	/** @brief For each side the vector related the coordinates of the point considered middle of the side */
	static REAL MidSideNode[27][3] = {
		/*00*/{-1.,-1.,-1.},/*01*/{1.,-1.,-1.},/*02*/{1.,1.,-1.},/*03*/{-1.,1.,-1.},
		/*04*/{-1.,-1., 1.},/*05*/{1.,-1., 1.},/*06*/{1.,1., 1.},/*07*/{-1.,1., 1.},
		/*08*/{ 0.,-1.,-1.},/*09*/{1., 0.,-1.},/*10*/{0.,1.,-1.},/*11*/{-1.,0.,-1.},
		/*12*/{-1.,-1., 0.},/*13*/{1.,-1., 0.},/*14*/{1.,1., 0.},/*15*/{-1.,1., 0.},
		/*16*/{ 0.,-1., 1.},/*17*/{1., 0., 1.},/*18*/{0.,1., 1.},/*19*/{-1.,0., 1.},
		/*20*/{ 0., 0.,-1.},/*21*/{0.,-1., 0.},/*22*/{1.,0., 0.},/*23*/{ 0.,1., 0.},
		/*24*/{-1., 0., 0.},/*25*/{0., 0., 1.},/*26*/{0.,0., 0.} };
	
	void TPZCube::LowerDimensionSides(int side,TPZStack<int> &smallsides)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	void TPZCube::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
		}
	}
	
	void TPZCube::HigherDimensionSides(int side, TPZStack<int> &high)
	{
		if(side <0 || side >= NSides) {
			PZError << "TPZCube::HigherDimensionSides side "<< side << endl;
		}
		int is;
		for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
		
	}
	
	int TPZCube::NSideNodes(int side)
	{
		return nsidenodes[side];
	}
	
	int TPZCube::SideNodeLocId(int side, int node)
	{
		if(side<8 && node == 0) return side;
		if(side>=8 && side < 20 && node < 2) return SideNodes[side-8][node];
		if(side>=20 && side < 26 && node < 4) return FaceNodes[side-20][node];
		if(side == 26 && node < 8) return node;
		PZError << "TPZCube::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
		return -1;
	}
	
	void TPZCube::CenterPoint(int side, TPZVec<REAL> &center) {
		center.Resize(Dimension);
		int i;
		for(i=0; i<Dimension; i++) {
			center[i] = MidSideNode[side][i];
		}
	}
	
	int TPZCube::SideDimension(int side) {
		if(side<0 || side >= NSides) {
			PZError << "TPZCube::SideDimension side " << side << endl;
			return -1;
		}
		return sidedimension[side];
	}
	
	TPZTransform TPZCube::SideToSideTransform(int sidefrom, int sideto)
	{
		if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
			PZError << "TPZCube::HigherDimensionSides sidefrom "<< sidefrom << 
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
		PZError << "TPZCube::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
		return TPZTransform(0);
	}
	
	TPZTransform TPZCube::TransformElementToSide(int side){
		
		if(side<0 || side>26){
			PZError << "TPZCube::TransformElementToSide called with side error\n";
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
			case 5 :
			case 6:
			case 7:
				return t;
			case  8:
			case 16:
				t.Mult()(0,0) =  1.0;
				return t;
			case  9:
			case 17:
				t.Mult()(0,1) =  1.0;
				return t;
			case 10:
			case 18:
				t.Mult()(0,0) = -1.0;
				return t;
			case 11:
			case 19:
				t.Mult()(0,1) = -1.0;
				return t;
			case 12:
			case 13:
			case 14:
			case 15:
				t.Mult()(0,2) = 1.0;
				return t;
			case 20:
			case 25:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				return t;
			case 21:
			case 23:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,2) =  1.0;
				return t;
			case 22:
			case 24:
				t.Mult()(0,1) =  1.0;
				t.Mult()(1,2) =  1.0;
				return t;
			case 26:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,2) =  1.0;
				return t;
		}
		return TPZTransform(0,0);
	}
	
	TPZTransform TPZCube::TransformSideToElement(int side){
		
		if(side<0 || side>26){
			PZError << "TPZCube::TransformSideToElement side out range\n";
			return TPZTransform(0,0);
		}
		TPZTransform t(3,sidedimension[side]);
		t.Mult().Zero();
		t.Sum().Zero();
		
		switch(side){
			case 0:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) = -1.0;
				t.Sum()(2,0) = -1.0;
				return t;
			case 1:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) = -1.0;
				t.Sum()(2,0) = -1.0;
				return t;
			case 2:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) = -1.0;
				return t;
			case 3:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) = -1.0;
				return t;
			case 4:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) = -1.0;
				t.Sum()(2,0) =  1.0;
				return t;
			case 5:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) = -1.0;
				t.Sum()(2,0) =  1.0;
				return t;
			case 6:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) =  1.0;
				return t;
			case 7:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) =  1.0;
				return t;
			case 8:
				t.Mult()(0,0) =  1.0;
				t.Sum()(1,0)  = -1.0;
				t.Sum()(2,0)  = -1.0;
				return t;
			case 9:
				t.Mult()(1,0) =  1.0;
				t.Sum()(0,0)  =  1.0;
				t.Sum()(2,0)  = -1.0;
				return t;
			case 10:
				t.Mult()(0,0) = -1.0;
				t.Sum()(1,0)  =  1.0;
				t.Sum()(2,0)  = -1.0;
				return t;
			case 11:
				t.Mult()(1,0) = -1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(2,0)  = -1.0;
				return t;
			case 12:
				t.Mult()(2,0) =  1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(1,0)  = -1.0;
				return t;
			case 13:
				t.Mult()(2,0) =  1.0;
				t.Sum()(0,0)  =  1.0;
				t.Sum()(1,0)  = -1.0;
				return t;
			case 14:
				t.Mult()(2,0) =  1.0;
				t.Sum()(0,0)  =  1.0;
				t.Sum()(1,0)  =  1.0;
				return t;
			case 15:
				t.Mult()(2,0) =  1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(1,0)  =  1.0;
				return t;
			case 16:
				t.Mult()(0,0) =  1.0;
				t.Sum()(1,0)  = -1.0;
				t.Sum()(2,0)  =  1.0;
				return t;
			case 17:
				t.Mult()(1,0) =  1.0;
				t.Sum()(0,0)  =  1.0;
				t.Sum()(2,0)  =  1.0;
				return t;
			case 18:
				t.Mult()(0,0) = -1.0;
				t.Sum()(1,0)  =  1.0;
				t.Sum()(2,0)  =  1.0;
				return t;
			case 19:
				t.Mult()(1,0) = -1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(2,0)  =  1.0;
				return t;
			case 20:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Sum()(2,0)  = -1.0;
				return t;
			case 21:
				t.Mult()(0,0) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(1,0)  = -1.0;
				return t;
			case 22:
				t.Mult()(1,0) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(0,0)  =  1.0;
				return t;
			case 23:
				t.Mult()(0,0) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(1,0)  =  1.0;
				return t;
			case 24:
				t.Mult()(1,0) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(0,0)  = -1.0;
				return t;
			case 25:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Sum()(2,0)  =  1.0;
				return t;
			case 26:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,2) =  1.0;
				return t;
		}
		return TPZTransform(0,0);
	}
	
	
	TPZIntPoints *TPZCube::CreateSideIntegrationRule(int side, int order){
		
		if(side<0 || side>26) {
			PZError << "TPZCube::CreateSideIntegrationRule. bad side number.\n";
			return 0;
		}
		if(side<8)   return new TPZInt1Point(order);            // sides 0 to 7 are vertices (corners)
		if(side<20)  return new TPZInt1d(order);           // sides 8 to 19 are lines
		if(side<26)  {                                     // sides 20 to 25 are quadrilaterals
			return new TPZIntQuad(order,order);
		}
		if(side==26) {
			return new IntruleType(order,order,order);     // integration of the element
		}
		return 0;
		
	}
	
	
	MElementType TPZCube::Type()
	{
		return ECube;
	}
	
	MElementType TPZCube::Type(int side)
	{
		switch(side) {
			case 0:
			case 1:
			case 2:
			case 3:
			case 4:
			case 5:
			case 6:
			case 7:
				return EPoint;
			case 8:
			case 9:
			case 10:
			case 11:
			case 12:
			case 13:
			case 14:
			case 15:
			case 16:
			case 17:
			case 18:
			case 19:        
				return EOned;
			case 20:
			case 21:
			case 22:
			case 23:
			case 24:
			case 25:
				return EQuadrilateral;
			case 26:
				return ECube;
			default:
				return ENoType;
		}
	}
	
	
	int TPZCube::NumSides() {
		return 27;
	}
	
	
	int TPZCube::NContainedSides(int side) {
		if(side<0) return -1;
		if(side<8) return 1;//cantos : 0 a 7
		if(side<20)	return 3;//lados : 8 a 19
		if(side<26)	return 9;//faces : 20 a 25
		if(side==26)	return 27;//centro : 26
		return -1;
	}
	/**
	 return number of sides of dimension dimension
	 **/
	int TPZCube::NumSides(int dimension) {
		if(dimension<0 || dimension> 3) {
			PZError << "TPZCube::NumSides. Bad parameter i.\n";
			return 0;
		}
		if(dimension==0) return 8;
		if(dimension==1) return 12;
		if(dimension==2) return 6;
		if(dimension==3) return 1;
		return -1;
	}
	
	// Pronto 23/04/98
	int TPZCube::ContainedSideLocId(int side, int node) {
		if(side<0 || side>26) return -1;
		if(side<8) {
			if(node==0) return side;
		} 
		else if(side<12) {//8,9,10,11
			int s = side-8;//0,1,2,3
			if(!node) return s;
			if(node==1) return (s+1)%4;
			if(node==2) return side;
		} 
		else if(side<16) {//12,13,14,15
			int s = side-12;//0,1,2,3
			if(!node) return s;
			if(node==1) return s+4;
			if(node==2) return side;
		} 
		else if(side<20) {//16,17,18,19
			int s = side-16;//0,1,2,3
			if(!node) return s+4;
			if(node==1) return (s+1)%4+4;
			if(node==2) return side;
		} 
		else if(side<26) {//20 a 25
			int s = side-20;
			if(node<9) return FaceConnectLocId[s][node];
		} 
		else if(side==26){
			return node;
		}
		PZError << "TPZShapeCube::ContainedSideLocId called for node = "
		<< node << " and side = " << side << "\n";
		return -1;
	}
	
	bool TPZCube::IsInParametricDomain(TPZVec<REAL> &pt, REAL tol){
		const REAL qsi = pt[0];
		const REAL eta = pt[1];
		const REAL zeta = pt[2];
		if( ( fabs(qsi) <= 1. + tol ) && ( fabs(eta) <= 1. + tol ) && ( fabs(zeta) <= 1. + tol ) ){
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
	int TPZCube::GetTransformId(TPZVec<int> &id)
	{
		LOGPZ_ERROR(logger,"GetTransformId not implemented")
		return -1;
	}
	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	int TPZCube::GetTransformId(int side, TPZVec<int> &id)
	{
		switch (side) {
			case 0:
			case 1:
			case 2:
			case 3:
			case 4:
			case 5:
			case 6:
			case 7:
				return 0;
				break;
			case 8:
			case 9:
			case 10:
			case 11:
			case 12:
			case 13:
			case 14:
			case 15:
			case 16:
			case 17:
			case 18:
			case 19:
			{
				int in1 = ContainedSideLocId(side,0);
				int in2 = ContainedSideLocId(side,1);
				return id[in1]<id[in2] ? 0 : 1;
			}
				break;
			case 20:
			case 21:
			case 22:
			case 23:
			case 24:
			case 25:
			{
				TPZManVector<int,4> locid;
				int i;
				for(i=0; i<4; i++) locid[i] = id[ContainedSideLocId(side,i)];
				return pztopology::TPZQuadrilateral::GetTransformId(locid);
				//			return pzshape::TPZShapeQuad::GetTransformId2dQ(locid);
			}
				break;			
			case 26:
				LOGPZ_ERROR(logger,"Please Implement me")
				return -1;
			default:
				break;
		}
		return -1;
	}
	
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
	void TPZCube::GetSideHDivPermutation(int side, TPZVec<int> &id, TPZVec<int> &permgather)
	{
	}

}
