
#include "TPZGeoLinear.h"
#include "pzelg1d.h"
#include "pzelgpoint.h"
#include "pzquad.h"
#include "pzshapelinear.h"
#include "pzgeoel.h"


void TPZGeoLinear::X(TPZFMatrix &coord,TPZVec<REAL> &loc,TPZVec<REAL> &result){

  int numnodes = NNodes;
  REAL spacephi[9], spacedphi[18];
  TPZFMatrix phi(numnodes,1,spacephi,9),dphi(1,numnodes,spacedphi,18);
  Shape/*1d*/(loc/*[0]*/,/*numnodes,*/phi,dphi);
  int in;
  for(in=0; in<3; in++) result[in] = 0.;
  for(in = 0; in < numnodes; in++) {
    int ic;
    for(ic=0; ic<3 ; ic++) {
      result[ic] += coord(ic,in)*phi(in,0);
    }
  }
}

void TPZGeoLinear::Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {
	REAL x = pt[0];
    phi(0,0) = (1-x)/2.;
    phi(1,0) = (1+x)/2.;
    dphi(0,0) = -0.5;
    dphi(0,1) = 0.5;
}

void TPZGeoLinear::Jacobian(TPZFMatrix coord,TPZVec<REAL> &param,TPZFMatrix &jacobian,
							TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) {


  REAL spacephi[9], spacedphi[18];
  TPZFMatrix phi(NNodes,1,spacephi,9),dphi(1,NNodes,spacedphi,18);
  Shape(param,phi,dphi);

  int ic;
  TPZVec<REAL> v1(3,0.);
  REAL mod1 = 0;

  for(int i=0; i < NNodes; i++) {
    for(ic = 0; ic < 3; ic++) {
      v1[ic] += coord(ic,i)*dphi(0,i);
    }
  }

  for(ic=0; ic<3; ic++) {
    mod1 += v1[ic]*v1[ic];
  }
  mod1 = sqrt(mod1);
  jacobian(0,0) = mod1;
  detjac = mod1;
  jacinv(0,0) = 1./detjac;

  axes.Zero();
  for(ic=0; ic<3; ic++) {
	  axes(0,ic) = v1[ic]/mod1;
  }
}

TPZGeoEl *TPZGeoLinear::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc){
  if(side==2) {
      TPZManVector<int> nodes(2);
      nodes[0] = orig->SideNodeIndex(side,0);
      nodes[1] = orig->SideNodeIndex(side,1);
      TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::SideConnectLocId(side,0)));
      TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::SideConnectLocId(side,1)));
      TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
      return gel;
  }
  else if(side==0 || side==1) {
      TPZManVector<int> nodeindexes(1);
      TPZGeoElPoint *gel;
      nodeindexes[0] = orig->SideNodeIndex(side,0);
      gel = new TPZGeoElPoint(nodeindexes,bc,*(orig->Mesh()));
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

