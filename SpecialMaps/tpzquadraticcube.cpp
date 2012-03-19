/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticCube methods. 
 */
#include "tpzquadraticcube.h"
#include "pzshapecube.h"
#include "tpzgeoblend.h"
#include "tpzgeoelmapped.h"

#include "pzgeoelrefless.h.h"
#include "tpzgeoelrefpattern.h.h"
#include "pznoderep.h.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.specialmaps.quadraticcube"));
#endif

using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;

void TPZQuadraticCube::Shape(TPZVec<REAL> &param,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
	
	REAL qsi = param[0], eta = param[1], zeta = param[2];
	
	phi(0,0)   =  1./8.*((-1. + eta)*(-1. + qsi)*(-1. + zeta)*(2. + eta + qsi + zeta));
	phi(1,0)   = -1./8.*((-1. + eta)*(1. + qsi)*(-1. + zeta)*(2. + eta - qsi + zeta));
	phi(2,0)   = -1./8.*((1. + eta)*(1. + qsi)*(-2. + eta + qsi - zeta)*(-1. + zeta));
	phi(3,0)   =  1./8.*((1 + eta)*(-1 + qsi)*(-2 + eta - qsi - zeta)*(-1 + zeta));
	
	phi(4,0)   = -1./8.*((-1. + eta)*(-1. + qsi)*(2. + eta + qsi - zeta)*(1. + zeta));
	phi(5,0)   =  1./8.*((-1. + eta)*(1. + qsi)*(2. + eta - qsi - zeta)*(1. + zeta));
	phi(6,0)   =  1./8.*((1. + eta)*(1. + qsi)*(1. + zeta)*(-2. + eta + qsi + zeta));
	phi(7,0)   = -1./8.*((1. + eta)*(-1. + qsi)*(1. + zeta)*(-2. + eta - qsi + zeta));
	
	phi(8,0)   = -1./4.*((-1. + eta)*(-1. + qsi*qsi)*(-1. + zeta));
	phi(9,0)   =  1./4.*((-1. + eta*eta)*(1. + qsi)*(-1. + zeta));
	phi(10,0)  =  1./4.*((1. + eta)*(-1. + qsi*qsi)*(-1. + zeta));
	phi(11,0)  = -1./4.*((-1. + eta*eta)*(-1. + qsi)*(-1. + zeta));
	
	phi(12,0)  = -1./4.*((-1. + eta)*(-1. + qsi)*(-1. + zeta*zeta));
	phi(13,0)  =  1./4.*((-1. + eta)*(1. + qsi)*(-1. + zeta*zeta));
	phi(14,0)  = -1./4.*((1. + eta)*(1. + qsi)*(-1. + zeta*zeta));
	phi(15,0)  =  1./4.*((1. + eta)*(-1. + qsi)*(-1. + zeta*zeta));
	
	phi(16,0)  =  1./4.*((-1. + eta)*(-1. + qsi*qsi)*(1. + zeta));
	phi(17,0)  = -1./4.*((-1. + eta*eta)*(1. + qsi)*(1. + zeta));
	phi(18,0)  = -1./4.*((1. + eta)*(-1. + qsi*qsi)*(1. + zeta));
	phi(19,0)  =  1./4.*((-1. + eta*eta)*(-1. + qsi)*(1. + zeta));
	
	
	dphi(0,0)  =  1./8.*((-1. + eta)*(-1. + zeta)*(1. + eta + 2.*qsi + zeta));
	dphi(1,0)  =  1./8.*((-1. + qsi)*(-1. + zeta)*(1. + 2.*eta + qsi + zeta));
	dphi(2,0)  =  1./8.*((-1. + eta)*(-1. + qsi)*(1. + eta + qsi + 2.*zeta));
	
	dphi(0,1)  = -1./8.*((-1. + eta)*(-1. + zeta)*(1 + eta - 2.*qsi + zeta));
	dphi(1,1)  =  1./8.*((1. + qsi)*(-1. - 2.*eta + qsi - zeta)*(-1. + zeta));
	dphi(2,1)  = -1./8.*((-1. + eta)*(1. + qsi)*(1. + eta - qsi + 2.*zeta));
	
	dphi(0,2)  =  1./8.*((1. + eta)*(-1. + zeta)*(1. - eta - 2.*qsi + zeta));
	dphi(1,2)  =  1./8.*((1. + qsi)*(-1. + zeta)*(1. - 2.*eta - qsi + zeta));
	dphi(2,2)  = -1./8.*((1. + eta)*(1. + qsi)*(-1. + eta + qsi - 2.*zeta));
	
	dphi(0,3)  =  1./8.*((1. + eta)*(-1. + eta - 2.*qsi - zeta)*(-1. + zeta));
	dphi(1,3)  = -1./8.*((-1. + qsi)*(-1. + zeta)*(1. - 2.*eta + qsi + zeta));
	dphi(2,3)  =  1./8.*((1. + eta)*(-1. + qsi)*(-1. + eta - qsi - 2.*zeta));
	
	dphi(0,4)  = -1./8.*((-1. + eta)*(1. + eta + 2.*qsi - zeta)*(1. + zeta));
	dphi(1,4)  = -1./8.*((-1. + qsi)*(1. + 2.*eta + qsi - zeta)*(1. + zeta));
	dphi(2,4)  = -1./8.*((-1. + eta)*(-1. + qsi)*(1. + eta + qsi - 2.*zeta));
	
	dphi(0,5)  =  1./8.*((-1. + eta)*(1. + eta - 2.*qsi - zeta)*(1. + zeta));
	dphi(1,5)  = -1./8.*((1. + qsi)*(1. + zeta)*(-1. - 2.*eta + qsi + zeta));
	dphi(2,5)  =  1./8.*((-1. + eta)*(1. + qsi)*(1. + eta - qsi - 2.*zeta));
	
	dphi(0,6)  =  1./8.*((1. + eta)*(1. + zeta)*(-1. + eta + 2.*qsi + zeta));
	dphi(1,6)  =  1./8.*((1. + qsi)*(1. + zeta)*(-1. + 2.*eta + qsi + zeta));
	dphi(2,6)  =  1./8.*((1. + eta)*(1. + qsi)*(-1. + eta + qsi + 2.*zeta));
	
	dphi(0,7)  = -1./8.*((1. + eta)*(1. + zeta)*(-1. + eta - 2.*qsi + zeta));
	dphi(1,7)  =  1./8.*((-1. + qsi)*(1. - 2.*eta + qsi - zeta)*(1. + zeta));
	dphi(2,7)  =  1./8.*((1. + eta)*(-1. + qsi)*(1. - eta + qsi - 2.*zeta));
	
	dphi(0,8)  = -1./2.*((-1. + eta)*qsi*(-1. + zeta));
	dphi(1,8)  = -1./4.*((-1. + qsi*qsi)*(-1. + zeta));
	dphi(2,8)  = -1./4.*((-1. + eta)*(-1. + qsi*qsi));
	
	dphi(0,9)  =  1./4.*((-1. + eta*eta)*(-1. + zeta));
	dphi(1,9)  =  1./2.*(eta*(1. + qsi)*(-1. + zeta));
	dphi(2,9)  =  1./4.*((-1. + eta*eta)*(1. + qsi));
	
	dphi(0,10) =  1./2.*((1. + eta)*qsi*(-1. + zeta));
	dphi(1,10) =  1./4.*((-1. + qsi*qsi)*(-1. + zeta));
	dphi(2,10) =  1./4.*((1. + eta)*(-1. + qsi*qsi));
	
	dphi(0,11) = -1./4.*((-1. + eta*eta)*(-1. + zeta));
	dphi(1,11) = -1./2.*(eta*(-1. + qsi)*(-1. + zeta));
	dphi(2,11) = -1./4.*((-1. + eta*eta)*(-1. + qsi));
	
	dphi(0,12) = -1./4.*((-1. + eta)*(-1. + zeta*zeta));
	dphi(1,12) = -1./4.*((-1. + qsi)*(-1. + zeta*zeta));
	dphi(2,12) = -1./2.*((-1. + eta)*(-1. + qsi)*zeta);
	
	dphi(0,13) =  1./4.*((-1. + eta)*(-1. + zeta*zeta));
	dphi(1,13) =  1./4.*((1. + qsi)*(-1. + zeta*zeta));
	dphi(2,13) =  1./2.*((-1. + eta)*(1. + qsi)*zeta);
	
	dphi(0,14) = -1./4.*((1. + eta)*(-1. + zeta*zeta));
	dphi(1,14) = -1./4.*((1. + qsi)*(-1. + zeta*zeta));
	dphi(2,14) = -1./2.*((1. + eta)*(1. + qsi)*zeta);
	
	dphi(0,15) =  1./4.*((1. + eta)*(-1. + zeta*zeta));
	dphi(1,15) =  1./4.*((-1. + qsi)*(-1. + zeta*zeta));
	dphi(2,15) =  1./2.*((1. + eta)*(-1. + qsi)*zeta);
	
	dphi(0,16) =  1./2.*((-1. + eta)*qsi*(1. + zeta));
	dphi(1,16) =  1./4.*((-1. + qsi*qsi)*(1. + zeta));
	dphi(2,16) =  1./4.*((-1. + eta)*(-1. + qsi*qsi));
	
	dphi(0,17) = -1./4.*((-1. + eta*eta)*(1. + zeta));
	dphi(1,17) = -1./2.*(eta*(1. + qsi)*(1. + zeta));
	dphi(2,17) = -1./4.*((-1. + eta*eta)*(1. + qsi));
	
	dphi(0,18) = -1./2.*((1. + eta)*qsi*(1. + zeta));
	dphi(1,18) = -1./4.*((-1. + qsi*qsi)*(1. + zeta));
	dphi(2,18) = -1./4.*((1. + eta)*(-1. + qsi*qsi));
	
	dphi(0,19) =  1./4.*((-1. + eta*eta)*(1. + zeta));
	dphi(1,19) =  1./2.*(eta*(-1. + qsi)*(1. + zeta));
	dphi(2,19) =  1./4.*((-1. + eta*eta)*(-1. + qsi));
}

void TPZQuadraticCube::X(TPZFMatrix<REAL> & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result) {
	
	TPZFNMatrix<20> phi(20,1);
	TPZFNMatrix<60> dphi(3,20);
	Shape(loc,phi,dphi);
	
	for(int i=0; i<3; i++)
	{
		result[i] = 0.0;
		for(int j=0; j<20; j++) 
		{
			result[i] += phi(j,0)*coord(i,j); 
		}
	}
}

void TPZQuadraticCube::Jacobian(TPZFMatrix<REAL> & coord, TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) {
#ifdef DEBUG
	if (NNodes != 20) {
		PZError << "TPZQuadraticCube.jacobian only implemented for 20, NumberOfNodes = " << NNodes << "\n";
	}
	
	if( fabs(param[0]) > 1.001 || fabs(param[1]) > 1.001 || fabs(param[2]) > 1.001) {
		PZError << "TPZQuadraticCube.jacobian. param out of range : "
		" param.NElements() = " << param.NElements() <<
		"\nparam[0] = " << param[0] << " param[1] = " << param[1] << "\n";
	}
#endif
	
	jacobian.Resize(3,3); axes.Resize(3,3); jacinv.Resize(3,3);
    jacobian.Zero(); axes.Zero();
    for(int d = 0; d < 3; d++) axes(d,d) = 1.;
	
	REAL spacephi[20]; REAL spacedphi[60];
	TPZFMatrix<REAL> phi(20,1,spacephi,20);
	TPZFMatrix<REAL> dphi(3,20,spacedphi,60);
	Shape(param,phi,dphi);
	for(int i = 0; i < 20; i++) {
		for(int j = 0; j < 3; j++) {
			jacobian(j,0) += coord(j,i)*dphi(0,i);
			jacobian(j,1) += coord(j,i)*dphi(1,i);
			jacobian(j,2) += coord(j,i)*dphi(2,i);
		}
	}
	
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

TPZGeoEl *TPZQuadraticCube::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                             TPZVec<int>& nodeindexes,
                                             int matid,
                                             int& index)
{
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

TPZGeoEl *TPZQuadraticCube::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) 
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



///CreateGeoElement -> TPZQuadraticCube

#define TPZGEOELEMENTQUADRATICCUBEID 3521 //verify the validity and range validity of this ID.
template<>
int TPZGeoElRefPattern<TPZQuadraticCube>::ClassId() const {
	return TPZGEOELEMENTQUADRATICCUBEID;
}
template class 
TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticCube>, TPZGEOELEMENTQUADRATICCUBEID>;


template class pzgeom::TPZNodeRep<20,TPZQuadraticCube>;
template class TPZGeoElRefLess<TPZQuadraticCube>;
