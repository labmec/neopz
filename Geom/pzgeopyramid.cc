// TPZGeoPyramid.c: implementation of the TPZGeoPyramid class.
//
//////////////////////////////////////////////////////////////////////

#include "pzgeopyramid.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "pzelgpoint.h"
#include "pzelg1d.h"
#include "pzelgt2d.h"
#include "pzelgq2d.h"
#include "pzelgpi3d.h"
#include "pzelgt3d.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"


void TPZGeoPyramid::Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {
	REAL T0xz = .5*(1.-pt[2]-pt[0]) / (1.-pt[2]);
	REAL T0yz = .5*(1.-pt[2]-pt[1]) / (1.-pt[2]);
	REAL T1xz = .5*(1.-pt[2]+pt[0]) / (1.-pt[2]);
	REAL T1yz = .5*(1.-pt[2]+pt[1]) / (1.-pt[2]);
	REAL lmez = (1.-pt[2]);
	phi(0,0)  = T0xz*T0yz*lmez;
	phi(1,0)  = T1xz*T0yz*lmez;
	phi(2,0)  = T1xz*T1yz*lmez;
	phi(3,0)  = T0xz*T1yz*lmez;
	phi(4,0)  = pt[2];
	REAL lmexmez = 1.-pt[0]-pt[2];
	REAL lmeymez = 1.-pt[1]-pt[2];
	REAL lmaxmez = 1.+pt[0]-pt[2];
	REAL lmaymez = 1.+pt[1]-pt[2];
	dphi(0,0) = -.25*lmeymez / lmez;
	dphi(1,0) = -.25*lmexmez / lmez;
	dphi(2,0) = -.25*(lmeymez+lmexmez-lmexmez*lmeymez/lmez) / lmez;

	dphi(0,1) =  .25*lmeymez / lmez;
	dphi(1,1) = -.25*lmaxmez / lmez;
	dphi(2,1) = -.25*(lmeymez+lmaxmez-lmaxmez*lmeymez/lmez) / lmez;

	dphi(0,2) =  .25*lmaymez / lmez;
	dphi(1,2) =  .25*lmaxmez / lmez;
	dphi(2,2) = -.25*(lmaymez+lmaxmez-lmaxmez*lmaymez/lmez) / lmez;

	dphi(0,3) = -.25*lmaymez / lmez;
	dphi(1,3) =  .25*lmexmez / lmez;
	dphi(2,3) = -.25*(lmaymez+lmexmez-lmexmez*lmaymez/lmez) / lmez;

	dphi(0,4) =  0.0;
	dphi(1,4) =  0.0;
	dphi(2,4) =  1.0;
}

void TPZGeoPyramid::Jacobian(TPZFMatrix & coord, TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){

	int nnodes = NNodes;

#ifdef DEBUG
	if (nnodes != 5) {
		PZError <<	"TPZGeoElPi3d.jacobian only implemented for"
					" 5 nodes, NumberOfNodes = " << nnodes << "\n";
	}
	if(param.NElements() != 3 || param[0] < -1. || param[0] > 1. ||
		param[1] < -1. || param[1] > 1. || param[2] < 0. || param[2] > 1.) {
		PZError << "TPZGeoElPi3d.jacobian. param out of range : "
					" param.NElements() = " << param.NElements() <<
					"\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2] << "\n";
		return;
	}
#endif

	REAL spacephi[5];
	TPZFMatrix phi(5,1,spacephi,5);
	REAL spacedphi[15];
	TPZFMatrix dphi(3,5,spacedphi,15);
	Shape(param,phi,dphi);
	jacobian.Zero();
	
	int i,j;
	for(i=0;i<5;i++) {
		for(j=0;j<3;j++) {
			jacobian(j,0) += coord(j,i)*dphi(0,i);
			jacobian(j,1) += coord(j,i)*dphi(1,i);
			jacobian(j,2) += coord(j,i)*dphi(2,i);
		}
	}

	detjac = -jacobian(0,2)*jacobian(1,1)*jacobian(2,0);//-a02 a11 a20
	detjac += jacobian(0,1)*jacobian(1,2)*jacobian(2,0);//+ a01 a12 a20
	detjac += jacobian(0,2)*jacobian(1,0)*jacobian(2,1);//+ a02 a10 a21
	detjac -= jacobian(0,0)*jacobian(1,2)*jacobian(2,1);//- a00 a12 a21
	detjac -= jacobian(0,1)*jacobian(1,0)*jacobian(2,2);//- a01 a10 a22
	detjac += jacobian(0,0)*jacobian(1,1)*jacobian(2,2);//+ a00 a11 a22

	jacinv(0,0) = (-jacobian(1,2)*jacobian(2,1)+jacobian(1,1)*jacobian(2,2))/detjac;//-a12 a21 + a11 a22
	jacinv(0,1) = ( jacobian(0,2)*jacobian(2,1)-jacobian(0,1)*jacobian(2,2))/detjac;//a02 a21 - a01 a22
	jacinv(0,2) = (-jacobian(0,2)*jacobian(1,1)+jacobian(0,1)*jacobian(1,2))/detjac;//-a02 a11 + a01 a12
	jacinv(1,0) = ( jacobian(1,2)*jacobian(2,0)-jacobian(1,0)*jacobian(2,2))/detjac;//a12 a20 - a10 a22
	jacinv(1,1) = (-jacobian(0,2)*jacobian(2,0)+jacobian(0,0)*jacobian(2,2))/detjac;//-a02 a20 + a00 a22
	jacinv(1,2) = ( jacobian(0,2)*jacobian(1,0)-jacobian(0,0)*jacobian(1,2))/detjac;//a02 a10 - a00 a12
	jacinv(2,0) = (-jacobian(1,1)*jacobian(2,0)+jacobian(1,0)*jacobian(2,1))/detjac;//-a11 a20 + a10 a21
	jacinv(2,1) = ( jacobian(0,1)*jacobian(2,0)-jacobian(0,0)*jacobian(2,1))/detjac;//a01 a20 - a00 a21
	jacinv(2,2) = (-jacobian(0,1)*jacobian(1,0)+jacobian(0,0)*jacobian(1,1))/detjac;//-a01 a10 + a00 a11

	axes.Zero();
	axes(0,0) = 1.;
	axes(1,1) = 1.;
	axes(2,2) = 1.;
}

void TPZGeoPyramid::X(TPZFMatrix & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result){
	REAL spacephi[10],spacedphi[20];
	int i,j;
	TPZFMatrix phi(5,1,spacephi,10);
	TPZFMatrix dphi(3,5,spacedphi,20);
	Shape(loc,phi,dphi);
	for(j=0;j<3;j++) {
		result[j] = 0.0;
		for(i=0;i<5;i++) result[j] += coord(j,i)*phi(i,0);
	}
}

TPZGeoEl *TPZGeoPyramid::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
	if(side<0 || side>18) {
		cout << "TPZGeoPyramid::CreateBCGeoEl Bad parameter side = " 
			 << side << "not implemented\n";
		return 0;
	}

   if(side==18) {
		cout << "TPZGeoElPi3d::CreateBCCompEl with side = 18 not implemented\n";
		return 0;
   }
   
	if(side<5) {
		TPZManVector<int> nodeindexes(1);
		TPZGeoElPoint *gel;
		nodeindexes[0] = orig->NodeIndex(side);
		gel = new TPZGeoElPoint(nodeindexes,bc,*orig->Mesh());
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	} 
	else if (side > 4 && side < 13) {//side =5 a 12 : lados
		TPZManVector<int> nodes(2);
		nodes[0] = orig->SideNodeIndex(side,0);
		nodes[1] = orig->SideNodeIndex(side,1);
		TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapePiram::SideConnectLocId(side,0)));
		TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapePiram::SideConnectLocId(side,1)));
		TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	}
	else if (side > 12) {//side = 13 a 17 : faces
		TPZManVector<int> nodes(4);//4o = -1 para face triangular
		int iside;
		for (iside=0;iside<4;iside++){
			nodes[iside] = orig->SideNodeIndex(side,iside);
		}
		TPZGeoElT2d *gelt;
		TPZGeoElQ2d *gelq;
		if(side==13) {
      		gelq = new TPZGeoElQ2d(nodes,bc,*orig->Mesh());
			for (iside=0; iside<8; iside++){
				TPZGeoElSide(gelq,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapePiram::SideConnectLocId(side,iside)));
			}
			TPZGeoElSide(gelq,8).SetConnectivity(TPZGeoElSide(orig,side));
			return gelq;
		} 
		else {
			nodes.Resize(3);
			gelt = new TPZGeoElT2d(nodes,bc,*orig->Mesh());
			for (iside=0; iside<6; iside++){
				TPZGeoElSide(gelt,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapePiram::SideConnectLocId(side,iside)));
			}
			TPZGeoElSide(gelt,6).SetConnectivity(TPZGeoElSide(orig,side));
			return gelt;
		}
	} 
	else 
		PZError << "TPZGeoPyramid::CreateBCGeoEl. Side = " << side << endl;
	return 0;
}


TPZIntPoints * TPZGeoPyramid::CreateSideIntegrationRule(int side, int order){
	if(side<0 || side>18) {
		PZError << "TPZGeoPyramid::CreateSideIntegrationRule. bad side number.\n";
		return 0;
	}
	//SideOrder corrige sides de 5 a 18 para 0 a 13
	if(side<5)   return new TPZInt1Point();//cantos 0 a 3
	if(side<13)  return new TPZInt1d(order);//lados 5 a 12
	if(side==13) return new TPZIntQuad(order,order);
	if(side<18)  {//faces : 14 a 17
		return new TPZIntTriang(order);
	}
	if(side==18) {//integração do elemento
		return new TPZIntPyram3D(order);
	}
	return 0;
}
