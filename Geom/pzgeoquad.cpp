// TPZGeoQuad.cpp: implementation of the TPZGeoQuad class.
//
//////////////////////////////////////////////////////////////////////

#include "pzgeoquad.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
//#include "pzelgpoint.h"
//#include "pzelg1d.h"
//#include "pzelgq2d.h"
#include "pzshapequad.h"


void TPZGeoQuad::Shape(TPZVec<REAL> &param,TPZFMatrix &phi,TPZFMatrix &dphi) {

   REAL x=param[0], y=param[1];

   phi(0,0) = .25*(1.-x)*(1.-y);
   phi(1,0) = .25*(1.+x)*(1.-y);
   phi(2,0) = .25*(1.+x)*(1.+y);
   phi(3,0) = .25*(1.-x)*(1.+y);

   dphi(0,0) = .25*(y-1.);
   dphi(1,0) = .25*(x-1.);

   dphi(0,1) = .25*(1.-y);
   dphi(1,1) =-.25*(1.+x);

   dphi(0,2) = .25*(1.+y);
   dphi(1,2) = .25*(1.+x);

   dphi(0,3) =-.25*(1.+y);
   dphi(1,3) = .25*(1.-x);
}


void TPZGeoQuad::Jacobian(TPZFMatrix & coord, TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){

  int nnodes = NNodes;
#ifdef DEBUG
  if (nnodes != 4) {
    PZError << "TPZGeoQuad.jacobian only implemented for"
      " 4, 8 or 9 nodes, NumberOfNodes = " << nnodes << "\n";
  }
  if( param[0] < -1. || param[0] > 1. ||  //param.NElements() != 3 ||
     param[1] < -1. || param[1] > 1.) {
    PZError << "TPZGeoQuad.jacobian. param out of range : "
      " param.NElements() = " << param.NElements() <<
      "\nparam[0] = " << param[0] << " param[1] = " << param[1] << "\n";
    return;
  }
#endif
  REAL spacephi[4];
  TPZFMatrix phi(4,1,spacephi,4);
  REAL spacedphi[8];
  TPZFMatrix dphi(2,4,spacedphi,8);
  Shape(param,phi,dphi);
  int i,j;
  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      jacobian(i,j)=0.;

//  TPZGeoNode *np;
  TPZVec<REAL> V1(3,0.),V2(3,0.),V2til(3,0.),V3(3,0.);
  REAL V1Norm=0.,V1V2=0.,V2tilNorm=0.;
  for(i=0;i<nnodes;i++) {
    //np = NodePtr(i);
    for(j=0;j<3;j++) {
      //V1[j] += np->Coord(j)*dphi(0,i);
		V1[j] += coord(j,i)*dphi(0,i);
      //V2[j] += np->Coord(j)*dphi(1,i);
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

void TPZGeoQuad::X(TPZFMatrix & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result){
  REAL spacephi[4],spacedphi[8];
  int i,j;
  TPZFMatrix phi(4,1,spacephi,4);
  TPZFMatrix dphi(2,4,spacedphi,8);
  Shape(loc,phi,dphi);
  for(i=0;i<3;i++) {
    result[i] = 0.0;
    for(j=0;j<4;j++)
      //result[i] += phi(j,0)*NodePtr(j)->Coord(i);
	  result[i] += phi(j,0)*coord(i,j);
  }
}

TPZGeoEl *TPZGeoQuad::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
  if(side==8) {//8
    TPZManVector<int> nodes(4);
    int i; 
    for (i=0;i<4;i++) {
      nodes[i] = orig->SideNodeIndex(side,i);
    }
    int index;
    TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EQuadrilateral,nodes,bc,index);
    //    TPZGeoElQ2d *gel = new TPZGeoElQ2d(nodes,bc,*orig->Mesh());
    int iside;
    for (iside = 0; iside <8; iside++){
      TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapeQuad::SideConnectLocId(side,iside)));
    }	    
    TPZGeoElSide(gel,8).SetConnectivity(TPZGeoElSide(orig,side));
    return gel;
  }
  else if(side>-1 && side<4) {//side = 0,1,2,3
    TPZManVector<int> nodeindexes(1);
    //    TPZGeoElPoint *gel;
    nodeindexes[0] = orig->SideNodeIndex(side,0);
    int index;
    TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
    //    gel = new TPZGeoElPoint(nodeindexes,bc,*(orig->Mesh()));
    TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
    return gel;
  }
  else if(side>3 && side<8) {
    TPZManVector<int> nodes(2);
    nodes[0] = orig->SideNodeIndex(side,0);
    nodes[1] = orig->SideNodeIndex(side,1);
    //    TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
    int index;
    TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
    TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeQuad::SideConnectLocId(side,0)));
    TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeQuad::SideConnectLocId(side,1)));
    TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
    return gel;
  }
  else PZError << "TPZGeoQuad::CreateBCCompEl has no bc.\n";
  return 0;
}


TPZIntPoints * TPZGeoQuad::CreateSideIntegrationRule(int side, int order){
  if(side<0 || side>8) {
    PZError << "TPZGeoQuad::CreateSideIntegrationRule wrong side " << side << endl;
    return 0;
  }
  if(side<4) return new TPZInt1Point();
  if(side<8) return new TPZInt1d(order);
  if(side==8) return new TPZIntQuad(order,order);
  return 0;
}
