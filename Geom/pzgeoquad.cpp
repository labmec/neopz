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
//#include "pzshapequad.h"

//using namespace pzshape;
using namespace std;

namespace pzgeom {

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


#ifdef DEBUG
//  const int nnodes = NNodes;
//  if (nnodes != 4) {
//    PZError << "TPZGeoQuad.jacobian only implemented for"
//      " 4 nodes, NumberOfNodes = " << nnodes << "\n";
//  }
  if( param[0] < -1.001 || param[0] > 1.001 || param[1] < -1.001 || param[1] > 1.001) {
	PZError << "TPZGeoQuad.jacobian. param out of range : "
	  " param.NElements() = " << param.NElements() <<
	  "\nparam[0] = " << param[0] << " param[1] = " << param[1] << "\n";
	//return;
  }
#endif
  jacobian.Resize(2,2); axes.Resize(2,3); jacinv.Resize(2,2);
  TPZFNMatrix<4> phi(4,1);
  TPZFNMatrix<8> dphi(2,4);
  TPZFNMatrix<6> axest(3,2);
  Shape(param,phi,dphi);
  jacobian.Zero();

  int spacedim = coord.Rows();
  TPZFMatrix VecMatrix(3,2,0.);
  for(int i = 0; i < 4; i++) {
	for(int j = 0; j < spacedim; j++) {
		VecMatrix(j,0) += coord(j,i)*dphi(0,i);
        VecMatrix(j,1) += coord(j,i)*dphi(1,i);
    }
  }
  VecMatrix.GramSchmidt(axest,jacobian);
  axest.Transpose(&axes);
  detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
  if(detjac)
  {
    jacinv(0,0) =  jacobian(1,1)/detjac;
    jacinv(1,1) =  jacobian(0,0)/detjac;
    jacinv(0,1) = -jacobian(0,1)/detjac;
    jacinv(1,0) = -jacobian(1,0)/detjac;
  }
  else
  {
    jacinv.Zero();
  }
}

void TPZGeoQuad::X(TPZFMatrix & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result){
  REAL spacephi[4],spacedphi[8];
  int i,j;
  TPZFMatrix phi(4,1,spacephi,4);
  TPZFMatrix dphi(2,4,spacedphi,8);
  int space = coord.Rows();
  Shape(loc,phi,dphi);
  result.Fill(0.);
  for(i=0;i<space;i++) {
    for(j=0;j<4;j++)
      //result[i] += phi(j,0)*NodePtr(j)->Coord(i);
	  result[i] += phi(j,0)*coord(i,j);
  }
}

bool TPZGeoQuad::MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix &JacToSide) {
     REAL qsi = InternalPar[0]; REAL eta = InternalPar[1];
     if( (fabs(qsi) - 1.) > 1e-5 || (fabs(eta) - 1.) > 1e-5 )
     {
         cout << "Point (qsi,eta) = (" << qsi << "," << eta << ") is out of TPZGeoQuad Master Element Range!\n";
         cout << "See TPZGeoQuad::MapToSide() method!\n";
		 DebugStop();
         exit(-1);
     }
	bool regularmap = true;
     TPZTransform Transf = pztopology::TPZQuadrilateral::SideToSideTransform(TPZGeoQuad::NSides - 1, side);
     SidePar.Resize(SideDimension(side));
     Transf.Apply(InternalPar,SidePar);

     int R = Transf.Mult().Rows();
     int C = Transf.Mult().Cols();

     JacToSide.Resize(R,C);
     for(int i = 0; i < R; i++)
     {
          for(int j = 0; j < C; j++) JacToSide(i,j) = Transf.Mult()(i,j);
     }
	return regularmap;
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
			TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,pztopology::TPZQuadrilateral::SideConnectLocId(side,iside)));
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
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,pztopology::TPZQuadrilateral::SideConnectLocId(side,0)));
		TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,pztopology::TPZQuadrilateral::SideConnectLocId(side,1)));
    TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
    return gel;
  }
  else PZError << "TPZGeoQuad::CreateBCCompEl has no bc.\n";
  return 0;
}


};
