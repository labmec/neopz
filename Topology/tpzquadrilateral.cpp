/**
 * @file
 * @brief Contains the implementation of the TPZQuadrilateral methods. 
 */

#include "tpzquadrilateral.h"

#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"
#include "pzquad.h"
#include "pzeltype.h"

#include "pzcreateapproxspace.h"
#include "pzshapequad.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.topology.pzquadrilateral"));
#endif

using namespace std;

namespace pztopology {

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
		int nsidecon = NContainedSides(side);
		int is;
		for(is=0; is<nsidecon-1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	void TPZQuadrilateral::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
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
		if(side<4) return new TPZInt1Point(order);              // sides 0 to 3 are vertices (corners)
		if(side<8) return new TPZInt1d(order);             // sides 4 to 7 are lines
		if(side==8) return new IntruleType(order,order);   // integration of the element
		return 0;
	}
	
	
	TPZTransform TPZQuadrilateral::TransformElementToSide(int side){
		
		if(side<0 || side>8){
			PZError << "TPZShapeQuad::TransformElementToSide called with side error\n";
			return TPZTransform(0,0);
		}
		
		TPZTransform t(sidedimension[side],2);//t(dimto,2)
		t.Mult().Zero();
		t.Sum().Zero();
		
		switch(side){
			case 0:
			case 1:
			case 2:
			case 3:
				return t;
			case 4:
				t.Mult()(0,0) = 1.0;
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
	}//method
	
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
	
	int TPZQuadrilateral::NContainedSides(int side) {
		if(side<0 || side>8) {
			PZError << "TPZShapeQuad::NContainedSides. Bad parameter i.\n";
			return 0;
		}
		if(side<4) return 1;
		if(side<8) return 3;
		return 9;//Cedric
	}
	
	int TPZQuadrilateral::NumSides(int dimension) {
		if(dimension<0 || dimension> 2) {
			PZError << "TPZShapeQuad::NumSides. Bad parameter i.\n";
			return 0;
		}
		if(dimension==0) return 4;
		if(dimension==1) return 4;
		if(dimension==2) return 1;
		return -1;
		
	}
	
	/**It do not verify the values of the c*/
	// side é o lado do elemento, c é o noh do lado
	int TPZQuadrilateral::ContainedSideLocId(int side,int c) {
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
				PZError << "TPZShapeQuad::ContainedSideLocId, connect = " << c << endl;
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
	
	
	/**
	 * Method which identifies the transformation based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */
	int TPZQuadrilateral::GetTransformId(TPZVec<int> &id)
	{
		return pzshape::TPZShapeQuad::GetTransformId2dQ(id);
	}
	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	int TPZQuadrilateral::GetTransformId(int side, TPZVec<int> &id)
	{
		switch (side) {
			case 0:
			case 1:
			case 2:
			case 3:
				return 0;
				break;
			case 4:
			case 5:
			case 6:
			case 7:
			{
				int in1 = ContainedSideLocId(side,0);
				int in2 = ContainedSideLocId(side,1);
				return id[in1]<id[in2] ? 0 : 1;
			}
				break;
			case 8:
			{
				return pzshape::TPZShapeQuad::GetTransformId2dQ(id);
			}
				break;
			default:
				break;
		}
		LOGPZ_ERROR(logger,"Wrong side parameter")
		return -1;
	}
	
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
	void TPZQuadrilateral::GetSideHDivPermutation(int side, TPZVec<int> &id, TPZVec<int> &permgather)
	{
		switch (side) {
			case 0:
			case 1:
			case 2:
			case 3:
				permgather[0] = 0;
				break;
			case 4:
			case 5:
			case 6:
			case 7:
			{
				int in1 = ContainedSideLocId(side,0);
				int in2 = ContainedSideLocId(side,1);
				if(in1<in2)
				{
					permgather[0] = 0;
					permgather[1] = 1;
					permgather[2] = 2;
				} 
				else 
				{
					permgather[0] = 1;
					permgather[1] = 0;
					permgather[2] = 2;
				}
			}
				break;
			case 8:
			{
				int i;
				int tid = pzshape::TPZShapeQuad::GetTransformId2dQ(id);
				if(tid%2 == 0)
				{
					switch (tid)
					{
						case 0:
							for(i=0; i<9; i++) permgather[i] = i;
							break;
						case 2:
							for(i=0; i<4; i++) permgather[i] = (i+1)%4;
							for(i=4; i<8; i++) permgather[i] = 4+(i+1)%4;
							permgather[8] = 8;
							break;
						case 4:
							for(i=0; i<4; i++) permgather[i] = (i+2)%4;
							for(i=4; i<8; i++) permgather[i] = 4+(i+2)%4;
							permgather[8] = 8;
							break;
						case 6:
							for(i=0; i<4; i++) permgather[i] = (i+3)%4;
							for(i=4; i<8; i++) permgather[i] = 4+(i+3)%4;
							permgather[8] = 8;
							break;
					}
				}
				else
				{
					TPZManVector<int,4> invid(4);
					invid[0] = 0;
					invid[1] = 3;
					invid[2] = 2;
					invid[3] = 1;
					switch (tid) {
						case 1:
							for(i=0; i<4; i++) permgather[i] = invid[i];
							for(i=4; i<8; i++) permgather[i] = 4+invid[(i+0)%4];
							permgather[8] = 8;
							break;
						case 3:
							for(i=0; i<4; i++) permgather[i] = invid[(i+1)%4];
							for(i=4; i<8; i++) permgather[i] = 4+invid[(i+1)%4];
							permgather[8] = 8;
							break;
						case 5:
							for(i=0; i<4; i++) permgather[i] = invid[(i+2)%4];
							for(i=4; i<8; i++) permgather[i] = 4+invid[(i+2)%4];
							permgather[8] = 8;
							break;
						case 7:
							for(i=0; i<4; i++) permgather[i] = invid[(i+3)%4];
							for(i=4; i<8; i++) permgather[i] = 4+invid[(i+3)%4];
							permgather[8] = 8;
							break;
						default:
							break;
					}
				}
			}
				break;
			default:
				break;
		}
		LOGPZ_ERROR(logger,"Wrong side parameter")
	}
	
}
