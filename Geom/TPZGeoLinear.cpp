
#include "TPZGeoLinear.h"
//#include "pzelg1d.h"
//#include "pzelgpoint.h"
#include "pzquad.h"
#include "pzshapelinear.h"
#include "pzgeoel.h"



void TPZGeoLinear::Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {
	REAL x = pt[0];
    phi(0,0) = (1-x)/2.;
    phi(1,0) = (1+x)/2.;
    dphi(0,0) = -0.5;
    dphi(0,1) = 0.5;
}


TPZGeoEl *TPZGeoLinear::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc){
  if(side==2) {
      TPZManVector<int> nodes(2);
      nodes[0] = orig->SideNodeIndex(side,0);
      nodes[1] = orig->SideNodeIndex(side,1);
      int index;
      TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
      //      TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::SideConnectLocId(side,0)));
      TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::SideConnectLocId(side,1)));
      TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
      return gel;
  }
  else if(side==0 || side==1) {
      TPZManVector<int> nodeindexes(1);
      //      TPZGeoElPoint *gel;
      nodeindexes[0] = orig->SideNodeIndex(side,0);
      int index;
      TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
      //      gel = new TPZGeoElPoint(nodeindexes,bc,*(orig->Mesh()));
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
	  return gel;
  }
  else PZError << "TPZGeoLinear::CreateBCGeoEl. Side = " << side << endl;
  return 0;
}

TPZIntPoints *TPZGeoLinear::CreateSideIntegrationRule(int side, int order) {

	if(side<0 || side>2) {
		PZError << "TPZGeoLinear::CreateSideIntegrationRule wrong side " << side << endl;
		return 0;
	}
	if(side != 2) return new TPZInt1Point();
	return new TPZInt1d(order);


}

