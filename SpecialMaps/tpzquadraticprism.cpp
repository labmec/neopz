/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticPrism methods. 
 */
#include "tpzquadraticprism.h"
#include "tpzgeoblend.h"
#include "tpzgeoelmapped.h"

#include "pzgeoelrefless.h.h"
#include "tpzgeoelrefpattern.h.h"
#include "pznoderep.h.h"
#include "pzshapepiram.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.specialmaps.quadraticprism"));
#endif

using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;

void TPZQuadraticPrism::Shape(TPZVec<REAL> &param,TPZFMatrix &phi,TPZFMatrix &dphi) {
	
	REAL qsi = param[0], eta = param[1], zeta = param[2];
	
	phi(0,0)   = -0.5*(-1. + eta + qsi)*(-1. + zeta)*(2.*(eta + qsi) + zeta);
	phi(1,0)   =  0.5*qsi*(-1. + zeta)*(2. - 2.*qsi + zeta);
	phi(2,0)   =  0.5*eta*(-1. + zeta)*(2. - 2.*eta + zeta);
	phi(3,0)   = -0.5*(-1. + eta + qsi)*(1. + zeta)*(-2.*(eta + qsi) + zeta);
	phi(4,0)   =  0.5*qsi*(1 + zeta)*(-2. + 2.*qsi + zeta);
	phi(5,0)   =  0.5*eta*(1. + zeta)*(-2. + 2.*eta + zeta);
	phi(6,0)   =  2.*qsi*(-1. + eta + qsi)*(-1. + zeta);
	phi(7,0)   = -2.*eta*qsi*(-1. + zeta);
	phi(8,0)   =  2.*eta*(-1. + eta + qsi)*(-1. + zeta);
	phi(9,0)   =  (-1. + eta + qsi)*(-1. + zeta*zeta);
	phi(10,0)  =  qsi - qsi*zeta*zeta;
	phi(11,0)  =  eta - eta*zeta*zeta;
	phi(12,0)  = -2.*qsi*(-1. + eta + qsi)*(1. + zeta);
	phi(13,0)  =  2.*eta*qsi*(1. + zeta);
    phi(14,0)  = -2.*eta*(-1. + eta + qsi)*(1. + zeta);
	//--------------------------------------
	dphi(0,0)  =  (1. - 2.*eta - 2.*qsi - 0.5*zeta)*(-1. + zeta);
	dphi(1,0)  =  (1. - 2.*eta - 2.*qsi - 0.5*zeta)*(-1. + zeta);
	dphi(2,0)  =  (-1. + eta + qsi)*(0.5 - eta - qsi - zeta);
	
	dphi(0,1)  = -1. + qsi*(2. - 2.*zeta) + 0.5*zeta + 0.5*zeta*zeta;
	dphi(1,1)  =  0.;
	dphi(2,1)  =  qsi*(0.5 - 1.*qsi + 1.*zeta);
	
	dphi(0,2)  =  0.;
	dphi(1,2)  = -1. + eta*(2. - 2.*zeta) + 0.5*zeta + 0.5*zeta*zeta;
	dphi(2,2)  =  eta*(0.5 - eta + zeta);
	
	dphi(0,3)  =  (-1. + 2.*eta + 2.*qsi - 0.5*zeta)*(1. + zeta);
	dphi(1,3)  =  (-1. + 2.*eta + 2.*qsi - 0.5*zeta)*(1. + zeta);
	dphi(2,3)  =  (-1. + eta + qsi)*(-0.5 + eta + qsi - zeta);
	
	dphi(0,4)  =  (-1. + 2.*qsi + 0.5*zeta)*(1. + zeta);
	dphi(1,4)  =  0.;
	dphi(2,4)  =  qsi*(-0.5 + qsi + zeta);
	
	dphi(0,5)  =  0.;
	dphi(1,5)  = -1. - 0.5*zeta + 0.5*zeta*zeta + eta*(2. + 2.*zeta);
	dphi(2,5)  =  eta*(-0.5 + eta + zeta);
	
	dphi(0,6)  =  (-2. + 2.*eta + 4.*qsi)*(-1. + zeta);
	dphi(1,6)  =  2.*qsi*(-1. + zeta);
	dphi(2,6)  =  2.*qsi*(-1. + eta + qsi);
	
	dphi(0,7)  = -2.*eta*(-1. + zeta);
	dphi(1,7)  = -2.*qsi*(-1. + zeta);
	dphi(2,7)  = -2.*eta*qsi;
	
	dphi(0,8)  =  2.*eta*(-1. + zeta);
	dphi(1,8)  =  (-2. + 4.*eta + 2.*qsi)*(-1. + zeta);
	dphi(2,8)  =  2.*eta*(-1. + eta + qsi);
	
	dphi(0,9)  = -1. + zeta*zeta;
	dphi(1,9)  = -1. + zeta*zeta;
	dphi(2,9)  =  2.*(-1. + eta + qsi)*zeta;
	
	dphi(0,10) =  1. - zeta*zeta;
	dphi(1,10) =  0.;
	dphi(2,10) = -2.*qsi*zeta;
	
	dphi(0,11) =  0.;
	dphi(1,11) =  1. - zeta*zeta;
	dphi(2,11) = -2.*eta*zeta;
	
	dphi(0,12) =  (2. - 2.*eta - 4.*qsi)*(1. + zeta);
	dphi(1,12) = -2.*qsi*(1. + zeta);
	dphi(2,12) = -2.*qsi*(-1. + eta + qsi);
    
    dphi(0,13) =  2.*eta*(1. + zeta);
	dphi(1,13) =  2.*qsi*(1. + zeta);
	dphi(2,13) =  2.*eta*qsi;
    
    dphi(0,14) = -2.*eta*(1. + zeta);
	dphi(1,14) =  (2. - 4.*eta - 2.*qsi)*(1. + zeta);
	dphi(2,14) = -2.*eta*(-1. + eta + qsi);
}

void TPZQuadraticPrism::X(TPZFMatrix & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result) {
	
	TPZFNMatrix<15> phi(15,1);
	TPZFNMatrix<45> dphi(3,15);
	Shape(loc,phi,dphi);
	
	for(int i=0; i<3; i++)
	{
		result[i] = 0.0;
		for(int j=0; j<15; j++) 
		{
			result[i] += phi(j,0)*coord(i,j); 
		}
	}
}

void TPZQuadraticPrism::Jacobian(TPZFMatrix & coord, TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) {
#ifdef DEBUG
	if (NNodes != 15) {
		PZError << "TPZQuadraticPrism.jacobian only implemented for 15, NumberOfNodes = " << NNodes << "\n";
	}
	
	if( param[0] < 0. || param[0] > 1.001 || param[1] < 0. || param[1] > 1.001-param[0] || fabs(param[2]) > 1.001) {
		PZError << "TPZQuadraticPrism.jacobian. param out of range : "
		" param.NElements() = " << param.NElements() <<
		"\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2]<< "\n";
	}
#endif
	
	jacobian.Resize(3,3); axes.Resize(3,3); jacinv.Resize(3,3);
    jacobian.Zero(); axes.Zero();
    for(int d = 0; d < 3; d++) axes(d,d) = 1.;
	
	REAL spacephi[15]; REAL spacedphi[45];
	TPZFMatrix phi(15,1,spacephi,15);
	TPZFMatrix dphi(3,15,spacedphi,45);
	Shape(param,phi,dphi);	
	for(int i = 0; i < 15; i++) {
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

TPZGeoEl *TPZQuadraticPrism::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                             TPZVec<int>& nodeindexes,
                                             int matid,
                                             int& index)
{
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

TPZGeoEl *TPZQuadraticPrism::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) 
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



///CreateGeoElement -> TPZQuadraticPrism

#define TPZGEOELEMENTQUADRATICPRISM 3535 //verify the validity and range validity of this ID.
template<>
int TPZGeoElRefPattern<TPZQuadraticPrism>::ClassId() const {
	return TPZGEOELEMENTQUADRATICPRISM;
}
template class 
TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticPrism>, TPZGEOELEMENTQUADRATICPRISM>;


template class pzgeom::TPZNodeRep<15,TPZQuadraticPrism>;
template class TPZGeoElRefLess<TPZQuadraticPrism>;
