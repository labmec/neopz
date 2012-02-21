/**
 * @file
 * @brief Contains the implementation of the TPZTriangle methods. 
 */

#include "tpztriangle.h"
#include "pzquad.h"

#include "pzshapetriang.h"
#include "pzcreateapproxspace.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.topology.pztriangle"));
#endif

using namespace std;

namespace pztopology {
	
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
	}//method
	
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
		int nsidecon = NContainedSides(side);
		int is;
		for(is=0; is<nsidecon-1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	void TPZTriangle::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
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
	
	//Tentando criar o metodo
	int TPZTriangle::NumSides(int dimension) {	
		if(dimension<0 || dimension> 2) {
			PZError << "TPZTriangle::NumSides. Bad parameter i.\n";
			return 0;
		}
		if(dimension==0) return 3;
		if(dimension==1) return 3;
		if(dimension==2) return 1;
		return -1;
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
	
	TPZTransform TPZTriangle::TransformElementToSide(int side) {
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
		if(side<3) return new TPZInt1Point(order);     // sides 0 to 2 are vertices (corners)
		if(side<6) return new TPZInt1d(order);    // sides 3 to 5 are lines
		if(side==6)return new IntruleType(order); // integration of the element
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
	
	
	int TPZTriangle::NumSides() {
		return NSides;
	}
	
	
	
	int TPZTriangle::NContainedSides(int side) {
		if(side<0 || side>6) {
			PZError << "TPZShapeTriang::NContainedSides. Bad parameter i.\n";
			return 0;
		}
		if(side<3) return 1;
		if(side<6) return 3;
		return 7;
	}
	
	/**It do not verify the values of the c*/
	// side Ž o lado do elemento, c Ž o noh do lado
	int TPZTriangle::ContainedSideLocId(int side, int c) {
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
				PZError << "TPZShapeTriang::ContainedSideLocId, connect = " << c << endl;
				return -1;
		}
	}
	
	/**
	 * Method which identifies the transformation based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */
	int TPZTriangle::GetTransformId(TPZVec<int> &id)
	{
		return pzshape::TPZShapeTriang::GetTransformId2dT(id);
	}
	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	int TPZTriangle::GetTransformId(int side, TPZVec<int> &id)
	{
		switch (side) {
			case 0:
			case 1:
			case 2:
				return 0;
				break;
			case 3:
			case 4:
			case 5:
			{
				int in1 = ContainedSideLocId(side,0);
				int in2 = ContainedSideLocId(side,1);
				return id[in1]<id[in2] ? 0 : 1;
			}
				break;
			case 6:
			{
				return pzshape::TPZShapeTriang::GetTransformId2dT(id);
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
	void TPZTriangle::GetSideHDivPermutation(int side, TPZVec<int> &id, TPZVec<int> &permgather)
	{
		switch (side) {
			case 0:
			case 1:
			case 2:
				permgather[0] = 0;
				break;
			case 3:
			case 4:
			case 5:
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
			case 6:
			{
				int i;
				int tid = pzshape::TPZShapeTriang::GetTransformId2dT(id);
				if(tid%2 == 0)
				{
					switch (tid)
					{
						case 0:
							for(i=0; i<7; i++) permgather[i] = i;
							break;
						case 2:
							for(i=0; i<3; i++) permgather[i] = (i+1)%3;
							for(i=4; i<6; i++) permgather[i] = 3+(i+1)%3;
							permgather[6] = 6;
							break;
						case 4:
							for(i=0; i<3; i++) permgather[i] = (i+2)%3;
							for(i=4; i<6; i++) permgather[i] = 3+(i+2)%3;
							permgather[6] = 6;
							break;
					}
				}
				else
				{
					TPZManVector<int,3> invid(3);
					invid[0] = 0;
					invid[1] = 2;
					invid[2] = 1;
					switch (tid) {
						case 1:
							for(i=0; i<3; i++) permgather[i] = invid[i];
							for(i=4; i<6; i++) permgather[i] = 3+invid[(i+0)%3];
							permgather[6] = 6;
							break;
						case 3:
							for(i=0; i<3; i++) permgather[i] = invid[(i+1)%3];
							for(i=4; i<6; i++) permgather[i] = 3+invid[(i+1)%3];
							permgather[6] = 6;
							break;
						case 5:
							for(i=0; i<3; i++) permgather[i] = invid[(i+2)%3];
							for(i=4; i<6; i++) permgather[i] = 3+invid[(i+2)%3];
							permgather[6] = 6;
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
