// TPZGeoTriangle.c: implementation of the TPZGeoTriangle class.
//
//////////////////////////////////////////////////////////////////////

#include "pzgeotriangle.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "pzshapetriang.h"


using namespace pzshape;
using namespace std;

namespace pzgeom {

MElementType TPZGeoTriangle::Type()
{
  return ETriangle;
}

MElementType TPZGeoTriangle::Type(int side)
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

void TPZGeoTriangle::Shape(TPZVec<REAL> &param,TPZFMatrix &phi,TPZFMatrix &dphi) {
	REAL x = param[0], y = param[1];
	phi(0,0) = 1.-x-y;
	phi(1,0) = x;
	phi(2,0) = y;
	dphi(0,0) = dphi(1,0) = -1.;
	dphi(0,1) = dphi(1,2) = 1.;
	dphi(1,1) = dphi(0,2) = 0.;
}


void TPZGeoTriangle::Jacobian(TPZFMatrix & coord, TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){

	REAL spacephi[3],spacedphi[6];
	int i,j;
	TPZFMatrix phi(3,1,spacephi,3);
	TPZFMatrix dphi(2,3,spacedphi,6);
	jacobian.Zero();
	Shape(param,phi,dphi);

	TPZVec<REAL> V1(3,0.),V2(3,0.),V2til(3,0.),V3(3,0.);
	REAL V1Norm=0.,V1V2=0.,V2tilNorm=0.;
//	TPZGeoNode *np;

	for(i=0;i<3;i++) {
//		np = NodePtr(i);
		for(j=0;j<3;j++) {
			V1[j] += coord(j,i)*dphi(0,i);
			V2[j] += coord(j,i)*dphi(1,i);
		}
	}
	for(i=0;i<3;i++) {
		V1Norm += V1[i]*V1[i];
		V1V2 += V1[i]*V2[i];
	}
	V1Norm = sqrt(V1Norm);
	for(i=0;i<3;i++) {
		V1[i] /= V1Norm;
		V2til[i] = V2[i] - V1V2*V1[i]/V1Norm;
		V2tilNorm += V2til[i]*V2til[i];
	}
	V2tilNorm = sqrt(V2tilNorm);
	jacobian(0,0) = V1Norm;
	jacobian(0,1) = V1V2/V1Norm;
	jacobian(1,1) = V2tilNorm;
	for(i=0;i<3;i++) {
		axes(0,i) = V1[i];
		axes(1,i) = V2til[i]/V2tilNorm;
	}
	detjac = jacobian(0,0)*jacobian(1,1)-jacobian(1,0)*jacobian(0,1);
	jacinv(0,0) = jacobian(1,1)/detjac;
	jacinv(1,1) = jacobian(0,0)/detjac;
	jacinv(0,1) = -jacobian(0,1)/detjac;
	jacinv(1,0) = -jacobian(1,0)/detjac;

	axes(2,0) = axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1);
	axes(2,1) = -axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0);
	axes(2,2) = axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0);
}

void TPZGeoTriangle::X(TPZFMatrix & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result){

	REAL spacephi[3],spacedphi[6];
	TPZFMatrix phi(3,1,spacephi,3);
	TPZFMatrix dphi(2,3,spacedphi,6);
	//Shape(par,phi,dphi);
	Shape(loc,phi,dphi);

	int i,j;
	for(i=0;i<3;i++) {
		result[i] = 0.0;
		for(j=0;j<3;j++)
			//result[i] += phi(j,0)*NodePtr(j)->Coord(i);
			result[i] += phi(j,0)*coord(i,j);
	}
}

TPZGeoEl *TPZGeoTriangle::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
  if(side==6) {
		TPZManVector<int> nodes(3);
		int i;
		for (i=0;i<3;i++){
			nodes[i] = orig->SideNodeIndex(side,i);
		}
		//TPZGeoElT2d *gel = CreateGeoEl(nodes,bc,*Mesh());
		//		TPZGeoElT2d *gel = new TPZGeoElT2d (nodes,bc,*orig->Mesh());
		int index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(ETriangle,nodes,bc,index);
		int iside;
		for (iside = 0; iside <6; iside++){
			TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::SideConnectLocId(side,iside)));
		}	    
		TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(orig,side));
		//TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(this,6));
		return gel;
	}
	else if(side>-1 && side<3) {
		TPZManVector<int> nodeindexes(1);
		//		TPZGeoElPoint *gel;
		nodeindexes[0] = orig->SideNodeIndex(side,0);
		int index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
		//		gel = new TPZGeoElPoint(nodeindexes,bc,*(orig->Mesh()));
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	}
	else if(side > 2 && side < 6) {
		TPZManVector<int> nodes(2);
		nodes[0] = orig->SideNodeIndex(side,0);
		nodes[1] = orig->SideNodeIndex(side,1);
		int index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
		//		TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
		//Cesar 2003-01-07
		//		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::SideConnectLocId(side,0)));
		//		TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::SideConnectLocId(side,1)));
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::SideConnectLocId(side,0)));
		TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::SideConnectLocId(side,1)));
		TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	}
	else PZError << "TPZGeoTriangle::CreateBCGeoEl has no bc.\n";
	return 0;
}


TPZIntPoints * TPZGeoTriangle::CreateSideIntegrationRule(int side, int order){
	if(side < 0 || side>6) {
		PZError << "TPZGeoTriangle::CreateSideIntegrationRule wrong side " << side << endl;
		return 0;
	}
	if(side<3) return new TPZInt1Point();
	if(side<6) return new TPZInt1d(order);
	if(side==6)return new TPZIntTriang(order);
	return 0;
}

};
