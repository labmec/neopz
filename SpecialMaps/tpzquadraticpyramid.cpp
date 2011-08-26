/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticPyramid methods. 
 */
#include "tpzquadraticpyramid.h"
#include "tpzgeoblend.h"
#include "tpzgeoelmapped.h"

#include "pzgeoelrefless.h.h"
#include "tpzgeoelrefpattern.h.h"
#include "pznoderep.h.h"
#include "pzshapepiram.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.specialmaps.quadraticpyramid"));
#endif

using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;

void TPZQuadraticPyramid::Shape(TPZVec<REAL> &param,TPZFMatrix &phi,TPZFMatrix &dphi) {
	
	REAL qsi = param[0], eta = param[1], zeta = param[2];
	
	phi(0,0)   =  ((-1. + eta)*(-1. + qsi)*(-1. + zeta)*(1. + eta + qsi + (2. + eta + qsi + 2.*eta*qsi)*zeta))/4.;
	phi(1,0)   =  ((-1. + eta)*(1. + qsi)*(-1. + zeta)*(-1. - eta + qsi + (-2. + qsi + eta*(-1. + 2.*qsi))*zeta))/4.;
	phi(2,0)   =  ((1. + eta)*(1. + qsi)*(-1. + zeta)*(1. + 2.*zeta - qsi*(1. + zeta) + eta*(-1. + (-1. + 2.*qsi)*zeta)))/4.;
	phi(3,0)   =  ((1. + eta)*(-1. + qsi)*(-1. + zeta)*(-1. + eta - qsi + (-2. + eta - qsi + 2.*eta*qsi)*zeta))/4.;
	phi(4,0)   =  zeta*(-1. + 2.*zeta);
	phi(5,0)   = -((-1. + eta)*(-1. + qsi*qsi)*(-1. + zeta)*(1. + eta*zeta))/2.;
	phi(6,0)   = -((-1. + eta*eta)*(1. + qsi)*(-1. + zeta)*(-1. + qsi*zeta))/2.;
	phi(7,0)   = -((1. + eta)*(-1. + qsi*qsi)*(-1. + zeta)*(-1. + eta*zeta))/2.;
	phi(8,0)   = -((-1. + eta*eta)*(-1. + qsi)*(-1. + zeta)*(1. + qsi*zeta))/2.;
	phi(9,0)   = -((-1. + eta)*(-1. + qsi)*(-1. + zeta)*zeta);
	phi(10,0)  =  (-1. + eta)*(1. + qsi)*(-1. + zeta)*zeta;
	phi(11,0)  = -((1. + eta)*(1. + qsi)*(-1. + zeta)*zeta);
	phi(12,0)  =  (1. + eta)*(-1. + qsi)*(-1. + zeta)*zeta;
	
	
	dphi(0,0)  =  ((-1. + eta)*(-1. + zeta)*(eta + zeta + eta*(-1. + 4*qsi)*zeta + 2.*qsi*(1. + zeta)))/4.;
	dphi(1,0)  =  ((-1. + qsi)*(-1. + zeta)*(qsi + zeta - qsi*zeta + 2.*eta*(1. + zeta + 2.*qsi*zeta)))/4.;
	dphi(2,0)  =  ((-1. + eta)*(-1. + qsi)*(-1. - 2.*eta*qsi + 2.*(2. + eta + qsi + 2.*eta*qsi)*zeta))/4.;
	
	dphi(0,1)  =  ((-1. + eta)*(-1. + zeta)*(-zeta + 2.*qsi*(1. + zeta) + eta*(-1. + zeta + 4*qsi*zeta)))/4.;
	dphi(1,1)  =  ((1. + qsi)*(-1. + zeta)*(qsi - (1. + qsi)*zeta + eta*(-2. + (-2. + 4*qsi)*zeta)))/4.;
	dphi(2,1)  =  ((-1. + eta)*(1. + qsi)*(1. - 2.*eta*qsi + 2.*(-2. + qsi + eta*(-1. + 2.*qsi))*zeta))/4.;
	
	dphi(0,2)  =  ((1. + eta)*(-1. + zeta)*(zeta - 2.*qsi*(1. + zeta) + eta*(-1. + zeta + 4*qsi*zeta)))/4.;
	dphi(1,2)  =  ((1. + qsi)*(-1. + zeta)*(qsi*(-1. + zeta) + zeta + eta*(-2. + (-2. + 4*qsi)*zeta)))/4.;
	dphi(2,2)  =  ((1. + eta)*(1. + qsi)*(-1. - 2.*eta*qsi + 2.*(2. - qsi + eta*(-1. + 2.*qsi))*zeta))/4.;
	
	dphi(0,3)  =  ((1. + eta)*(-1. + zeta)*(eta - zeta + eta*(-1. + 4*qsi)*zeta - 2.*qsi*(1. + zeta)))/4.;
	dphi(1,3)  =  ((-1. + qsi)*(-1. + zeta)*(qsi*(-1. + zeta) - zeta + 2.*eta*(1. + zeta + 2.*qsi*zeta)))/4.;
	dphi(2,3)  =  ((1. + eta)*(-1. + qsi)*(1. - 2.*(2. + qsi)*zeta + 2.*eta*(zeta + qsi*(-1. + 2.*zeta))))/4.;
	
	dphi(0,4)  =  0.;
	dphi(1,4)  =  0.;
	dphi(2,4)  = -1. + 4.*zeta;
	
	dphi(0,5)  = -((-1. + eta)*qsi*(-1. + zeta)*(1. + eta*zeta));
	dphi(1,5)  = -((-1. + qsi*qsi)*(-1. + zeta)*(1. + (-1. + 2.*eta)*zeta))/2.;
	dphi(2,5)  = -((-1. + eta)*(-1. + qsi*qsi)*(1. + eta*(-1. + 2.*zeta)))/2.;
	
	dphi(0,6)  = -((-1. + eta*eta)*(-1. + zeta)*(-1. + zeta + 2.*qsi*zeta))/2.;
	dphi(1,6)  = -(eta*(1. + qsi)*(-1. + zeta)*(-1. + qsi*zeta));
	dphi(2,6)  = -((-1. + eta*eta)*(1. + qsi)*(-1. + qsi*(-1. + 2.*zeta)))/2.;
	
	dphi(0,7)  = -((1. + eta)*qsi*(-1. + zeta)*(-1. + eta*zeta));
	dphi(1,7)  = -((-1. + qsi*qsi)*(-1. + zeta)*(-1. + zeta + 2.*eta*zeta))/2.;
	dphi(2,7)  = -((1. + eta)*(-1. + qsi*qsi)*(-1. + eta*(-1. + 2.*zeta)))/2.;
	
	dphi(0,8)  = -((-1. + eta*eta)*(-1. + zeta)*(1. + (-1. + 2.*qsi)*zeta))/2.;
	dphi(1,8)  = -(eta*(-1. + qsi)*(-1. + zeta)*(1. + qsi*zeta));
	dphi(2,8)  = -((-1. + eta*eta)*(-1. + qsi)*(1. + qsi*(-1. + 2.*zeta)))/2.;
	
	dphi(0,9)  = -((-1. + eta)*(-1. + zeta)*zeta);
	dphi(1,9)  = -((-1. + qsi)*(-1. + zeta)*zeta);
	dphi(2,9)  =  (-1. + eta)*(-1. + qsi)*(1. - 2.*zeta);
	
	dphi(0,10) =  (-1. + eta)*(-1. + zeta)*zeta;
	dphi(1,10) =  (1. + qsi)*(-1. + zeta)*zeta;
	dphi(2,10) =  (-1. + eta)*(1. + qsi)*(-1. + 2.*zeta);
	
	dphi(0,11) = -((1. + eta)*(-1. + zeta)*zeta);
	dphi(1,11) = -((1. + qsi)*(-1. + zeta)*zeta);
	dphi(2,11) =  (1. + eta)*(1. + qsi)*(1. - 2.*zeta);
	
	dphi(0,12) =  (1. + eta)*(-1. + zeta)*zeta;
	dphi(1,12) =  (-1. + qsi)*(-1. + zeta)*zeta;
	dphi(2,12) =  (1. + eta)*(-1. + qsi)*(-1. + 2.*zeta);
}

void TPZQuadraticPyramid::X(TPZFMatrix & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result) {
	
	TPZFNMatrix<13> phi(13,1);
	TPZFNMatrix<39> dphi(3,13);
	Shape(loc,phi,dphi);
	
	for(int i=0; i<3; i++)
	{
		result[i] = 0.0;
		for(int j=0; j<13; j++) 
		{
			result[i] += phi(j,0)*coord(i,j); 
		}
	}
}

void TPZQuadraticPyramid::Jacobian(TPZFMatrix & coord, TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) {
#ifdef DEBUG
	if (NNodes != 20) {
		PZError << "TPZGeoQuad.jacobian only implemented for 8, NumberOfNodes = " << NNodes << "\n";
	}
	
	if( param[0] < -1.001 || param[0] > 1.001 || param[1] < -1.001 || param[1] > 1.001) {
		PZError << "TPZGeoQuad.jacobian. param out of range : "
		" param.NElements() = " << param.NElements() <<
		"\nparam[0] = " << param[0] << " param[1] = " << param[1] << "\n";
	}
#endif
	
	jacobian.Resize(3,3); axes.Resize(3,3); jacinv.Resize(3,3);
	
	REAL spacephi[13]; REAL spacedphi[39];
	TPZFMatrix phi(13,1,spacephi,13);
	TPZFMatrix dphi(3,13,spacedphi,39);
	Shape(param,phi,dphi);
	jacobian.Zero();
	
	TPZFMatrix VecMatrix(3,3,0.);
	for(int i = 0; i < 13; i++) {
		for(int j = 0; j < 3; j++) {
			VecMatrix(j,0) += coord(j,i)*dphi(0,i);
			VecMatrix(j,1) += coord(j,i)*dphi(1,i);
			VecMatrix(j,2) += coord(j,i)*dphi(2,i);
		}
	}
	
	TPZFMatrix axest;
	VecMatrix.GramSchmidt(axest,jacobian);
	axest.Transpose(&axes);
	
	detjac = -jacobian(0,2)*jacobian(1,1)*jacobian(2,0)
	+ jacobian(0,1)*jacobian(1,2)*jacobian(2,0)
	+ jacobian(0,2)*jacobian(1,0)*jacobian(2,1)
	- jacobian(0,0)*jacobian(1,2)*jacobian(2,1)
	- jacobian(0,1)*jacobian(1,0)*jacobian(2,2)
	+ jacobian(0,0)*jacobian(1,1)*jacobian(2,2);
	
	if(IsZero(detjac))
	{
		std::stringstream sout;
		sout << "Singular Jacobian " << detjac;
		
#ifdef LOG4CXX
		LOGPZ_ERROR(logger , sout.str())
#endif
		
		detjac = ZeroTolerance();
	}
	jacinv(0,0) = (-jacobian(1,2)*jacobian(2,1)+jacobian(1,1)*jacobian(2,2))/detjac;//-a12 a21 + a11 a22
	jacinv(0,1) = ( jacobian(0,2)*jacobian(2,1)-jacobian(0,1)*jacobian(2,2))/detjac;// a02 a21 - a01 a22
	jacinv(0,2) = (-jacobian(0,2)*jacobian(1,1)+jacobian(0,1)*jacobian(1,2))/detjac;//-a02 a11 + a01 a12
	jacinv(1,0) = ( jacobian(1,2)*jacobian(2,0)-jacobian(1,0)*jacobian(2,2))/detjac;// a12 a20 - a10 a22
	jacinv(1,1) = (-jacobian(0,2)*jacobian(2,0)+jacobian(0,0)*jacobian(2,2))/detjac;//-a02 a20 + a00 a22
	jacinv(1,2) = ( jacobian(0,2)*jacobian(1,0)-jacobian(0,0)*jacobian(1,2))/detjac;// a02 a10 - a00 a12
	jacinv(2,0) = (-jacobian(1,1)*jacobian(2,0)+jacobian(1,0)*jacobian(2,1))/detjac;//-a11 a20 + a10 a21
	jacinv(2,1) = ( jacobian(0,1)*jacobian(2,0)-jacobian(0,0)*jacobian(2,1))/detjac;// a01 a20 - a00 a21
	jacinv(2,2) = (-jacobian(0,1)*jacobian(1,0)+jacobian(0,0)*jacobian(1,1))/detjac;//-a01 a10 + a00 a11
}


/**
 * Creates a geometric element according to the type of the father element
 */

TPZGeoEl *TPZQuadraticPyramid::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                             TPZVec<int>& nodeindexes,
                                             int matid,
                                             int& index)
{
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

TPZGeoEl *TPZQuadraticPyramid::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) 
{
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
}



///CreateGeoElement -> TPZQuadraticPyramid

#define TPZGEOELEMENTQUADRATICPYRAMIDID 3522 //verify the validity and range validity of this ID.
template<>
int TPZGeoElRefPattern<TPZQuadraticPyramid>::ClassId() const {
	return TPZGEOELEMENTQUADRATICPYRAMIDID;
}
template class 
TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticPyramid>, TPZGEOELEMENTQUADRATICPYRAMIDID>;


template class pzgeom::TPZNodeRep<13,TPZQuadraticPyramid>;
template class TPZGeoElRefLess<TPZQuadraticPyramid>;
