#include "TPZGeoCube.h"
#include "pzelgpoint.h"
#include "pzelg1d.h"
#include "pzelgq2d.h"
#include "pzshapecube.h"
#include "pzquad.h"

void TPZGeoCube::X(TPZFMatrix &nodes,TPZVec<REAL> & loc,TPZVec<REAL> &result){

  int nrow = nodes.Rows();
  int ncol = nodes.Cols();
  if(nrow != 3 || ncol != 8){//8x3 nós por linhas
	cout << "TPZGeoCube::X nodes matrix error size\n";
  }
  REAL spacephi[12],spacedphi[30];
  int i,j;
  TPZFMatrix phi(8,1,spacephi,12);
  TPZFMatrix dphi(3,8,spacedphi,30);
  Shape(loc,phi,dphi);
  for(j=0;j<3;j++) {
    result[j] = 0.0;
    for(i=0;i<8;i++) result[j] += nodes(j,i)*phi(i,0);
  }
}

void TPZGeoCube::Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {

  REAL x[2],dx[2],y[2],dy[2],z[2],dz[2];
  x[0] = (1.-pt[0])/2.;
  x[1] = (1.+pt[0])/2.;
  dx[0] = -0.5;
  dx[1] = 0.5;
  y[0] = (1.-pt[1])/2.;
  y[1] = (1.+pt[1])/2.;
  dy[0] = -0.5;
  dy[1] = 0.5;
  z[0] = (1.-pt[2])/2.;
  z[1] = (1.+pt[2])/2.;
  dz[0] = -0.5;
  dz[1] = 0.5;

  phi(0,0) = x[0]*y[0]*z[0];
  phi(1,0) = x[1]*y[0]*z[0];
  phi(2,0) = x[1]*y[1]*z[0];
  phi(3,0) = x[0]*y[1]*z[0];
  phi(4,0) = x[0]*y[0]*z[1];
  phi(5,0) = x[1]*y[0]*z[1];
  phi(6,0) = x[1]*y[1]*z[1];
  phi(7,0) = x[0]*y[1]*z[1];
  dphi(0,0) = dx[0]*y[0]*z[0];
  dphi(1,0) = x[0]*dy[0]*z[0];
  dphi(2,0) = x[0]*y[0]*dz[0];
  dphi(0,1) = dx[1]*y[0]*z[0];
  dphi(1,1) = x[1]*dy[0]*z[0];
  dphi(2,1) = x[1]*y[0]*dz[0];
  dphi(0,2) = dx[1]*y[1]*z[0];
  dphi(1,2) = x[1]*dy[1]*z[0];
  dphi(2,2) = x[1]*y[1]*dz[0];
  dphi(0,3) = dx[0]*y[1]*z[0];
  dphi(1,3) = x[0]*dy[1]*z[0];
  dphi(2,3) = x[0]*y[1]*dz[0];
  dphi(0,4) = dx[0]*y[0]*z[1];
  dphi(1,4) = x[0]*dy[0]*z[1];
  dphi(2,4) = x[0]*y[0]*dz[1];
  dphi(0,5) = dx[1]*y[0]*z[1];
  dphi(1,5) = x[1]*dy[0]*z[1];
  dphi(2,5) = x[1]*y[0]*dz[1];
  dphi(0,6) = dx[1]*y[1]*z[1];
  dphi(1,6) = x[1]*dy[1]*z[1];
  dphi(2,6) = x[1]*y[1]*dz[1];
  dphi(0,7) = dx[0]*y[1]*z[1];
  dphi(1,7) = x[0]*dy[1]*z[1];
  dphi(2,7) = x[0]*y[1]*dz[1];
}

void TPZGeoCube::Jacobian(TPZFMatrix &nodes,TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){

#ifdef DEBUG
  if (NNodes != 8) {
    PZError << "TPZGeoCube.jacobian only implemented for"
      " 8 nodes, NumberOfNodes = " << NNodes << "\n";
  }
  if(param.NElements() != 3 || param[0] < -1. || param[0] > 1. ||
     param[1] < -1. || param[1] > 1. || param[2] < -1. || param[2] > 1.) {
    PZError << "TPZGeoCube.jacobian. param out of range : "
      " param.NElements() = " << param.NElements() <<
      "\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2] << "\n";
    return;
  }
#endif

  int nrow = nodes.Rows();
  int ncol = nodes.Cols();
  if(nrow != 3 || ncol != 8){//8x3 nós por linhas
	cout << "TPZGeoCube::X nodes matrix error size\n";
  }

  REAL spacephi[12];
  TPZFMatrix phi(8,1,spacephi,12);
  REAL spacedphi[30];
  TPZFMatrix dphi(3,8,spacedphi,30);
  Shape(param,phi,dphi);
  jacobian.Zero();
  REAL coor;
  int i,j;
  for(i=0;i<8;i++) {
    for(j=0;j<3;j++) {
      coor = nodes(j,i); // ??
      jacobian(j,0) += coor*dphi(0,i);
      jacobian(j,1) += coor*dphi(1,i);
      jacobian(j,2) += coor*dphi(2,i);
    }
  }

  detjac = -jacobian(0,2)*jacobian(1,1)*jacobian(2,0);//- a02 a11 a20
  detjac += jacobian(0,1)*jacobian(1,2)*jacobian(2,0);//+ a01 a12 a20
  detjac += jacobian(0,2)*jacobian(1,0)*jacobian(2,1);//+ a02 a10 a21
  detjac -= jacobian(0,0)*jacobian(1,2)*jacobian(2,1);//- a00 a12 a21
  detjac -= jacobian(0,1)*jacobian(1,0)*jacobian(2,2);//- a01 a10 a22
  detjac += jacobian(0,0)*jacobian(1,1)*jacobian(2,2);//+ a00 a11 a22

  jacinv(0,0) = (-jacobian(1,2)*jacobian(2,1)+jacobian(1,1)*jacobian(2,2))/detjac;//-a12 a21 + a11 a22
  jacinv(0,1) = ( jacobian(0,2)*jacobian(2,1)-jacobian(0,1)*jacobian(2,2))/detjac;// a02 a21 - a01 a22
  jacinv(0,2) = (-jacobian(0,2)*jacobian(1,1)+jacobian(0,1)*jacobian(1,2))/detjac;//-a02 a11 + a01 a12
  jacinv(1,0) = ( jacobian(1,2)*jacobian(2,0)-jacobian(1,0)*jacobian(2,2))/detjac;// a12 a20 - a10 a22
  jacinv(1,1) = (-jacobian(0,2)*jacobian(2,0)+jacobian(0,0)*jacobian(2,2))/detjac;//-a02 a20 + a00 a22
  jacinv(1,2) = ( jacobian(0,2)*jacobian(1,0)-jacobian(0,0)*jacobian(1,2))/detjac;// a02 a10 - a00 a12
  jacinv(2,0) = (-jacobian(1,1)*jacobian(2,0)+jacobian(1,0)*jacobian(2,1))/detjac;//-a11 a20 + a10 a21
  jacinv(2,1) = ( jacobian(0,1)*jacobian(2,0)-jacobian(0,0)*jacobian(2,1))/detjac;// a01 a20 - a00 a21
  jacinv(2,2) = (-jacobian(0,1)*jacobian(1,0)+jacobian(0,0)*jacobian(1,1))/detjac;//-a01 a10 + a00 a11

  axes.Zero();
  axes(0,0) = 1.;
  axes(1,1) = 1.;
  axes(2,2) = 1.;
}

TPZGeoEl *TPZGeoCube::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc) {

  if(side<0 || side>26) {
	cout <<  "TPZGeoCube::CreateBCGeoEl unexpected side = " << side << endl;
    return 0;
  }
  if(side==26) {
    cout << "TPZGeoCube::CreateBCGeoEl with side = 26 not implemented\n";
    return 0;
  }

  if(side<8) {
    TPZManVector<int> nodeindexes(1);
    TPZGeoElPoint *gel;
    nodeindexes[0] = orig->NodeIndex(side);
    gel = new TPZGeoElPoint(nodeindexes,bc,*orig->Mesh());
	TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
    return gel;
  } 
  else 
    if(side > 7 && side < 20) {//side = 8 a 19 : arestas
      TPZManVector<int> nodes(2);
      nodes[0] = orig->SideNodeIndex(side,0);
      nodes[1] = orig->SideNodeIndex(side,1);
      TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeCube::SideConnectLocId(side,0)));
      TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeCube::SideConnectLocId(side,1)));
      TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
      return gel;
    } 
    else 
      if(side > 19) {//side = 20 a 25 : faces
		TPZManVector<int> nodes(4);
		int in;
		for (in=0;in<4;in++){
			nodes[in] = orig->SideNodeIndex(side,in);
		}
		TPZGeoElQ2d *gel = new TPZGeoElQ2d(nodes,bc,*orig->Mesh());
		for (in=0;in<8;in++){
			TPZGeoElSide(gel,in).SetConnectivity(TPZGeoElSide(orig,TPZShapeCube::SideConnectLocId(side,in)));
		}
		TPZGeoElSide(gel,8).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
      } 
      else 
		PZError << "TPZGeoCube::CreateBCGeoEl. Side = " << side << endl;
  return 0;
}

TPZIntPoints *TPZGeoCube::CreateSideIntegrationRule(int side, int order){

  if(side<0 || side>26) {
    PZError << "TPZGeoCube::CreateSideIntegrationRule. bad side number.\n";
    return 0;
  }
  //SideOrder corrige sides de 8 a 26 para 0 a 18
  if(side<8)   return new TPZInt1Point();//cantos 0 a 7
  if(side<20)  return new TPZInt1d(order);//lados 8 a 19
  if(side<26)  {//faces : 20 a 25
    return new TPZIntQuad(order,order);
  }
  if(side==26) {//integração do elemento
    return new TPZIntCube3D(order,order,order);
  }
  return 0;

}

