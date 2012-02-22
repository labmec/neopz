/**
 * @file
 * @brief Contains the implementation of the TPZLine methods. 
 */

#include "tpzline.h"

#include "pzerror.h"
#include "pzreal.h"
#include "pztrnsform.h"
#include "pzquad.h"
#include "pzeltype.h"

#include "pzlog.h"
#include "pzcreateapproxspace.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.topology.pzline"));
#endif

using namespace std;

namespace pztopology {

	static int nhighdimsides[3] = {1,1,0};
	
	static int sidedimension[3] = {0,0,1};
	
	static int highsides[3][1] = {
		{2},
		{2},
		{0}
	};
	
	static REAL sidetosidetransforms[3][1][4][3] = {
		{
			{{0,0,0},{0,0,0},{0,0,0},{-1,0,0}}
		},
		{
			{{0,0,0},{0,0,0},{0,0,0},{1,0,0}}
		},
		{
			{{0,0,0},{0,0,0},{0,0,0},{0,0,0}}}
	};
	
	static REAL MidSideNode[3][1] = {{-1.},{1.},{0.}};
	
	static int nsidenodes[3] = {1,1,2};
	
	int TPZLine::NSideNodes(int side)
	{
		return nsidenodes[side];
	}
	
	int TPZLine::NumSides(int dimension) {
		if(dimension<0 || dimension> 1) {
			PZError << "TPZLine::NumSides. Bad parameter i.\n";
			return 0;
		}
		if(dimension==0) return 2;
		if(dimension==1) return 1;
		return -1;
	}
	
	
	void TPZLine::LowerDimensionSides(int side,TPZStack<int> &smallsides)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		int is;
		for(is=0; is<nsidecon-1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	void TPZLine::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
		}
	}
	
	bool TPZLine::IsInParametricDomain(TPZVec<REAL> &pt, REAL tol) {
		const REAL qsi = pt[0];
		if( fabs(qsi) <= 1. + tol){
			return true;
		}
		else{
			return false;
		}  
	}//method
	
	void TPZLine::HigherDimensionSides(int side, TPZStack<int> &high)
	{
		if(side <0 || side >= NSides) {
			PZError << "TPZLine::HigherDimensionSides side "<< side << endl;
		}
		int is;
		for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
		
	}
	
	int TPZLine::SideNodeLocId(int side, int node)
	{
		if(side <2 && node == 0) return side;
		if(side == 2 && node <2) return node;
		PZError << "TPZLine::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
		return -1;
	}
	void TPZLine::CenterPoint(int side, TPZVec<REAL> &center) {
		//center.Resize(Dimension);
		int i;
		for(i=0; i<Dimension; i++) {
			center[i] = MidSideNode[side][i];
		}
	}
	
	int TPZLine::SideDimension(int side) {
		if(side<0 || side >= NSides) {
			PZError << "TPZLine::SideDimension side " << side << endl;
			return -1;
		}
		return sidedimension[side];
	}
	
	TPZTransform TPZLine::SideToSideTransform(int sidefrom, int sideto)
	{
		if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
			PZError << "TPZLine::HigherDimensionSides sidefrom "<< sidefrom << 
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
		PZError << "TPZLine::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
		return TPZTransform(0);
	}
	
	TPZTransform TPZLine::TransformElementToSide(int side){
		
		if(side<0 || side>2){
			PZError << "TPZLine::TransformElementToSide called with side error\n";
			return TPZTransform(0,0);
		}
		
		//  int sidedim = SideDimension(side);
		TPZTransform t(SideDimension(side),1);//t(dimto,2)
		t.Mult().Zero();	//TPZGeoElQ2d *gq;
		t.Sum().Zero();//int dimto = gq->SideDimension(side);
		
		switch(side){
			case 0:
			case 1:
				return t;
			case 2:
				t.Mult()(0,0) = 1.0;//par. var.
				return t;
		}
		return TPZTransform(0,0);
	}
	
	TPZTransform TPZLine::TransformSideToElement(int side){
		
		if(side<0 || side>2){
			PZError << "TPZLine::TransformSideToElement side out range\n";
			return TPZTransform(0,0);
		}
		int sidedim = 1;
		if(side <2) sidedim = 0;
		
		TPZTransform t(1,sidedim);
		t.Mult().Zero();
		t.Sum().Zero();
		
		switch(side){
			case 0:
				t.Sum()(0,0) = -1.0;
				return t;
			case 1:
				t.Sum()(0,0) =  1.0;
				return t;
			case 2:
				t.Mult()(0,0) =  1.0;
				return t;
		}
		return TPZTransform(0,0);
	}
	
	
	TPZIntPoints *TPZLine::CreateSideIntegrationRule(int side, int order) {
		
		if(side<0 || side>2) {
			PZError << "TPZLine::CreateSideIntegrationRule wrong side " << side << endl;
			return 0;
		}
		if(side != 2) return new TPZInt1Point(order);   // sides 0 and 1 are vertices (corners)
		return new IntruleType(order);
		
		
	}
	
	
	MElementType TPZLine::Type()
	{
		return EOned;
	}
	
	MElementType TPZLine::Type(int side)
	{
		switch(side) {
			case 0:
			case 1:
				return EPoint;
			case 2:
				return EOned;
			default:
				return ENoType;
		}
	}
	
	
	int TPZLine::NumSides() {
		return 3;
	}
	
	
	int TPZLine::NContainedSides(int i) {
		if(i==0 || i==1) return 1;
		else if(i==2) return 3;
		PZError << "TPZLine::NContainedSides. Bad parameter i = " << i << " .\n";
		return 0;
	}
	
	int TPZLine::ContainedSideLocId(int side,int c) {
		switch(side) {
			case 0:
			case 1:
				if(c != 0)
				{
					PZError << "TPZLine::ContainedSideLocId, connect = " << c << endl;
				}
				return side;
			case 2:
				return c;
			default:
				PZError << "TPZLine::ContainedSideLocId called with side = " << side << endl;
				return 0;
		}
	}
	
	/**
	 * Method which identifies the transformation based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */
	int TPZLine::GetTransformId(TPZVec<int> &id)
	{
		return id[0] < id[1] ? 0 : 1;
	}
	
	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	int TPZLine::GetTransformId(int side, TPZVec<int> &id)
	{
		switch (side) {
			case 0:
			case 1:
				return 0;
				break;
			case 2:
				return id[0] < id[1] ? 0 : 1;
			default:
				break;
		}
		LOGPZ_ERROR(logger,"Wrong input parameter")
		return -1;
	}
	
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
	void TPZLine::GetSideHDivPermutation(int side, TPZVec<int> &id, TPZVec<int> &permgather)
	{
		switch (side) {
			case 0:
			case 1:
				permgather[0] = 0;
				break;
			case 2:
				if(id[0]<id[1])
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
				
			default:
				break;
		}
		LOGPZ_ERROR(logger,"Wrong input parameter")
	}

}
