/**
 * @file
 * @brief Contains the implementation of the TPZGeoPoint methods. 
 */

#include "pzgeopoint.h"
#include "pzquad.h"
#include "pzfmatrix.h"
#include "pzstack.h"
#include "pzgeoel.h"
#include "tpzgeoelrefpattern.h"

using namespace std;

namespace pzgeom {
	
	void TPZGeoPoint::X(const TPZFMatrix<REAL> &coord,TPZVec<REAL> &loc,TPZVec<REAL> &result){
		int i;
		for (i=0;i<coord.Rows();i++){
			result[i] = coord.GetVal(i,0);
		}
	}
	
	void TPZGeoPoint::Shape(TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		phi(0,0) = 1.;
	}
	
	void TPZGeoPoint::Jacobian(const TPZFMatrix<REAL> &coord,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,
							   TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) {
		jacobian.Redim(0,0);
		jacinv.Redim(0,0);
		detjac = 1.;
		axes.Zero();
		axes.Redim(0,3);
		axes.Zero();
		/*axes(0,0) = 1.;
		axes(1,1) = 1.;
		axes(2,2) = 1.;*/
	}
	
	TPZGeoEl *TPZGeoPoint::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc){
		if(side==0) {
			TPZManVector<long> nodeindexes(1);
			nodeindexes[0] = orig->NodeIndex(0);
			long index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,0));
			return gel;
		}
		else PZError << "TPZGeoPoint::CreateBCGeoEl. Side = " << side << endl;
		return 0;
	}
	
	/** Creates a geometric element according to the type of the father element */
	TPZGeoEl *TPZGeoPoint::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
											TPZVec<long>& nodeindexes,
											int matid,
											long& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
};
