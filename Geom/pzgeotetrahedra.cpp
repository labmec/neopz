// TPZGeoTetrahedra.cpp: implementation of the TPZGeoTetrahedra class.
//
//////////////////////////////////////////////////////////////////////

#include "pzgeotetrahedra.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
//#include "pzelgpoint.h"
//#include "pzelg1d.h"
//#include "pzelgt2d.h"
//#include "pzelgq2d.h"
//#include "pzelgpi3d.h"
//#include "pzelgt3d.h"
#include "pzshapetetra.h"

using namespace pzshape;

using namespace pzshape;

namespace pzgeom {

void TPZGeoTetrahedra::Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {
  phi(0,0)  = 1-pt[0]-pt[1]-pt[2];
  phi(1,0)  = pt[0];
  phi(2,0)  = pt[1];
  phi(3,0)  = pt[2];

  dphi(0,0) = -1.0;
  dphi(1,0) = -1.0;
  dphi(2,0) = -1.0;
  dphi(0,1) =  1.0;
  dphi(1,1) =  0.0;
  dphi(2,1) =  0.0;
  dphi(0,2) =  0.0;
  dphi(1,2) =  1.0;
  dphi(2,2) =  0.0;
  dphi(0,3) =  0.0;
  dphi(1,3) =  0.0;
  dphi(2,3) =  1.0;
}

void TPZGeoTetrahedra::Jacobian(TPZFMatrix & coord, TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){

	int nnodes = NNodes;
#ifdef DEBUG
	if (nnodes != 4) {
		PZError << "TPZGeoTetrahedra.jacobian only implemented for"
		" 4 nodes, NumberOfNodes = " << nnodes << "\n";
	}
	if(param.NElements() != 3 || param[0] < 0. || param[0] > 1. ||
		param[1] < 0. || param[1] > 1. || param[2] < 0. || param[2] > 1.) {
		PZError << "TPZGeoTetrahedra.jacobian. param out of range : "
			" param.NElements() = " << param.NElements() <<
			"\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2] << "\n";
		return;
	}
#endif
	REAL spacephi[10];
	TPZFMatrix phi(4,1,spacephi,10);
	REAL spacedphi[20];
	TPZFMatrix dphi(3,4,spacedphi,20);
	Shape(param,phi,dphi);
	jacobian.Zero();
  
	int i,j;
	for(i=0;i<4;i++) {
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

void TPZGeoTetrahedra::X(TPZFMatrix & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result){
	REAL spacephi[10],spacedphi[20];
	int i,j;
	TPZFMatrix phi(4,1,spacephi,10);
	TPZFMatrix dphi(3,4,spacedphi,20);
	Shape(loc,phi,dphi);
	for(j=0;j<3;j++) {
		result[j] = 0.0;
		for(i=0;i<4;i++) 
			result[j] += coord(j,i)*phi(i,0);
	}
}

TPZGeoEl *TPZGeoTetrahedra::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
	if(side<0 || side>14){
		cout << "TPZGeoTetrahedra::CreateBCCompEl with bad side = " << side << "not implemented\n";	
		return 0;
	}

	if(side==14) {
		cout << "TPZGeoTetrahedra::CreateBCCompEl with side = 14 not implemented\n";
		return 0;
	}
	if(side<4) {
	  TPZManVector<int> nodeindexes(1);
	  //		TPZGeoElPoint *gel;
		nodeindexes[0] = orig->NodeIndex(side);
		int index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
		//		gel = new TPZGeoElPoint(nodeindexes,bc,*orig->Mesh());
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
   } else if (side > 3 && side < 10) {//side =4 a 9 : lados
 		TPZManVector<int> nodes(2);
		nodes[0] = orig->SideNodeIndex(side,0);
		nodes[1] = orig->SideNodeIndex(side,1);
		int index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
		//		TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeTetra::SideConnectLocId(side,0)));
		TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeTetra::SideConnectLocId(side,1)));
		TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	} else if (side > 9) {//side = 10 a 13 : faces
 		TPZManVector<int> nodes(3);
		int in;
		for (in=0;in<3;in++){
			nodes[in] = orig->SideNodeIndex(side,in);
		}
		int index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(ETriangle,nodes,bc,index);
		//		TPZGeoElT2d *gel = new TPZGeoElT2d(nodes,bc,*orig->Mesh());
		for (in=0;in<6;in++){
			TPZGeoElSide(gel,in).SetConnectivity(TPZGeoElSide(orig,TPZShapeTetra::SideConnectLocId(side,in)));
		}
		TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	} 
	else PZError << "TPZGeoTetrahedra::CreateBCGeoEl. Side = " << side << endl;
	return 0;
}


TPZIntPoints * TPZGeoTetrahedra::CreateSideIntegrationRule(int side, int order){
	if(side<0 || side>15) {
		PZError << "TPZGeoTetrahedra::CreateSideIntegrationRule. bad side number.\n";
   		return 0;
	}
	//SideOrder corrige sides de 4 a 14 para 0 a 10
	if(side<4)   return new TPZInt1Point();//cantos 0 a 3 : cria regra com um ponto
	if(side<10)  return new TPZInt1d(order);//lados 4 a 9
	if(side<14)  {//faces : 10 a 13
   		return new TPZIntTriang(order);
	}
	if(side==14) {//integração do elemento
   		return new TPZIntTetra3D(order);
	}
	return 0;
}

};
