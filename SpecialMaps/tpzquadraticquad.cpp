#include "tpzquadraticquad.h"
#include "pzshapequad.h"
#include "tpzgeoblend.h"
#include "tpzgeoelmapped.h"

#include "pzgeoelrefless.h.h"
#include "tpzgeoelrefpattern.h.h"
#include "pznoderep.h.h"

using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;

void TPZQuadraticQuad::Shape(TPZVec<REAL> &param,TPZFMatrix &phi,TPZFMatrix &dphi) {
   REAL qsi = param[0], eta = param[1];

   phi(0,0)  = -0.25*(-1. + eta)*(-1. + qsi)*(1. + eta + qsi);
   phi(1,0)  =  0.25*(-1. + eta)*(1. + eta - 1.*qsi)*(1. + qsi);
   phi(2,0)  =  0.25*( 1. + qsi)*(1. + eta)*(qsi + eta - 1.);
   phi(3,0)  = -0.25*( 1. + eta)*(-1. + eta - 1.*qsi)*(-1. + qsi);

   phi(4,0)  =  0.5*(-1. + eta)*(-1. + qsi)*( 1. + qsi);
   phi(5,0)  = -0.5*(-1. + eta)*( 1. + eta)*( 1. + qsi);
   phi(6,0)  = -0.5*( 1. + eta)*(-1. + qsi)*( 1. + qsi);
   phi(7,0)  =  0.5*(-1. + eta)*( 1. + eta)*(-1. + qsi);

   dphi(0,0) = -0.25*(-1. + eta)*(eta + 2.*qsi);
   dphi(1,0) = -0.25*(-1. + qsi)*(2.*eta + qsi);
   dphi(0,1) =  0.25*(-1. + eta)*(eta - 2.*qsi);
   dphi(1,1) =  0.25*(2.*eta - 1.*qsi)*(1. + qsi);
   dphi(0,2) =  0.25*(1. + eta)*(eta + 2.*qsi);
   dphi(1,2) =  0.25*(1. + qsi)*(2.*eta + qsi);
   dphi(0,3) = -0.25*(1. + eta)*(eta - 2.*qsi);
   dphi(1,3) = -0.25*(2.*eta - 1.*qsi)*(-1. + qsi);

   dphi(0,4) =  (-1. + eta)*qsi;
   dphi(1,4) =  0.5*(-1. + qsi)*(1. + qsi);
   dphi(0,5) = -0.5*(-1. + eta)*(1. + eta);
   dphi(1,5) = -1.*eta*(1. + qsi);
   dphi(0,6) = -1.*(1. + eta)*qsi;
   dphi(1,6) = -0.5*(-1. + qsi)*(1. + qsi);
   dphi(0,7) =  0.5*(-1. + eta)*(1. + eta);
   dphi(1,7) =  eta*(-1. + qsi);
}

void TPZQuadraticQuad::X(TPZFMatrix & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result) {

    TPZFNMatrix<9> phi(8,1);
    TPZFNMatrix<16> dphi(2,8);
    Shape(loc,phi,dphi);

    for(int i = 0; i < 3; i++) {
        result[i] = 0.0;
        for(int j = 0; j < 8; j++) result[i] += phi(j,0)*coord(i,j);
    }
}

void TPZQuadraticQuad::Jacobian(TPZFMatrix & coord, TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) {
  #ifdef DEBUG
  if (NNodes != 8) {
    PZError << "TPZGeoQuad.jacobian only implemented for 8, NumberOfNodes = " << NNodes << "\n";
  }

  if( param[0] < -1.001 || param[0] > 1.001 || param[1] < -1.001 || param[1] > 1.001) {
    PZError << "TPZGeoQuad.jacobian. param out of range : "
      " param.NElements() = " << param.NElements() <<
      "\nparam[0] = " << param[0] << " param[1] = " << param[1] << "\n";
  }
  #endif

  jacobian.Resize(2,2); axes.Resize(2,3); jacinv.Resize(2,2);

  REAL spacephi[8]; REAL spacedphi[16];
  TPZFMatrix phi(8,1,spacephi,8);
  TPZFMatrix dphi(2,8,spacedphi,16);
  Shape(param,phi,dphi);
  jacobian.Zero();

  TPZFMatrix VecMatrix(3,2,0.);
  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 3; j++) {
        VecMatrix(j,0) += coord(j,i)*dphi(0,i);
        VecMatrix(j,1) += coord(j,i)*dphi(1,i);
    }
  }

  TPZFMatrix axest;
  VecMatrix.GramSchmidt(axest,jacobian);
  axest.Transpose(&axes);

  detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
  jacinv(0,0) =  jacobian(1,1)/detjac;
  jacinv(1,1) =  jacobian(0,0)/detjac;
  jacinv(0,1) = -jacobian(0,1)/detjac;
  jacinv(1,0) = -jacobian(1,0)/detjac;
}

TPZGeoEl *TPZQuadraticQuad::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
	
	int ns = orig->NSideNodes(side);
	TPZManVector<int> nodeindices(ns);
	int in;
	for(in=0; in<ns; in++)
	{
		nodeindices[in] = orig->SideNodeIndex(side,in);
	}
	int index;
	
	TPZGeoMesh *mesh = orig->Mesh();
	MElementType type = orig->Type(side);
	
	TPZGeoEl *newel = mesh->CreateGeoBlendElement(type, nodeindices, bc, index);
	TPZGeoElSide me(orig,side);
	TPZGeoElSide newelside(newel,newel->NSides()-1);
	
	newelside.InsertConnectivity(me);
	newel->Initialize();
	
	return newel;
	
	
	
/*
	if(side==8) {//8
    TPZManVector<int> nodes(4);
    int i; 
    for (i=0;i<4;i++) {
      nodes[i] = orig->SideNodeIndex(side,i);
    }
    int index;
    TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EQuadrilateral,nodes,bc,index);
    int iside;
    for (iside = 0; iside <8; iside++){
      TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapeQuad::SideConnectLocId(side,iside)));
    }
    TPZGeoElSide(gel,8).SetConnectivity(TPZGeoElSide(orig,side));
    return gel;
  }
  else if(side>-1 && side<4) {//side = 0,1,2,3
    TPZManVector<int> nodeindexes(1);
    nodeindexes[0] = orig->SideNodeIndex(side,0); int index;
    TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
    TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
    return gel;
  }
  else if(side>3 && side<8) {
    TPZManVector<int> nodes(2);
    nodes[0] = orig->SideNodeIndex(side,0);
    nodes[1] = orig->SideNodeIndex(side,1); int index;
    TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
    TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeQuad::SideConnectLocId(side,0)));
    TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeQuad::SideConnectLocId(side,1)));
    TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
    return gel;
  }
  else PZError << "TPZGeoQuad::CreateBCCompEl has no bc.\n";
  return 0;
 */
}


/**
 * Creates a geometric element according to the type of the father element
 */

TPZGeoEl *TPZQuadraticQuad::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
									 TPZVec<int>& nodeindexes,
									 int matid,
									 int& index)
{
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}



///CreateGeoElement -> TPZQuadraticQuad

#define TPZGEOELEMENTQUADRATICQUADID 352
            template<>
                int TPZGeoElRefPattern<TPZQuadraticQuad>::ClassId() const {
              return TPZGEOELEMENTQUADRATICQUADID;
                }
                template class 
                    TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticQuad>, TPZGEOELEMENTQUADRATICQUADID>;

               
template class pzgeom::TPZNodeRep<8,TPZQuadraticQuad>;
template class TPZGeoElRefLess<TPZQuadraticQuad>;
