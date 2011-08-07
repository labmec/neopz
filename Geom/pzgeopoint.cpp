/**
 * @file
 * @brief Contains the implementation of the TPZGeoPoint methods. 
 */

#include "pzgeopoint.h"
//#include "pzelgpoint.h"
//#include "pzelg1d.h"
#include "pzquad.h"
#include "pzfmatrix.h"
#include "pzstack.h"
#include "pzgeoel.h"
#include "tpzgeoelrefpattern.h"

using namespace std;

namespace pzgeom {
	
	void TPZGeoPoint::X(TPZFMatrix &coord,TPZVec<REAL> &loc,TPZVec<REAL> &result){
		int i;
		for (i=0;i<coord.Rows();i++){
			result[i] = coord(i,0);
		}
	}
	
	bool TPZGeoPoint::MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix &JacToSide) {
		SidePar.Resize(0); JacToSide.Resize(0,0);
		return true;
	}
	
	void TPZGeoPoint::Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {
		phi(0,0) = 1.;
	}
	
	void TPZGeoPoint::Jacobian(TPZFMatrix coord,TPZVec<REAL> &param,TPZFMatrix &jacobian,
							   TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) {
		/*result*/jacobian.Redim(0,0);
		jacinv.Redim(0,0);
		detjac = 1.;
		axes.Zero();
		axes.Redim(3,3);
		axes.Zero();
		axes(0,0) = 1.;
		axes(1,1) = 1.;
		axes(2,2) = 1.;
	}
	
	TPZGeoEl *TPZGeoPoint::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc){
		if(side==0) {
			TPZManVector<int> nodeindexes(1);
			nodeindexes[0] = orig->NodeIndex(0);
			//TPZGeoElPoint *gel = CreateGeoEl(nodes,bc,*Mesh());
			int index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			
			//    TPZGeoElPoint *gel = new TPZGeoElPoint(nodeindexes,bc,*orig->Mesh());
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,0));
			return gel;
		}
		else PZError << "TPZGeoPoint::CreateBCGeoEl. Side = " << side << endl;
		return 0;
	}
	
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	TPZGeoEl *TPZGeoPoint::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
											TPZVec<int>& nodeindexes,
											int matid,
											int& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
	
	
};
