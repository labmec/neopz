//
//  TPZQuadSphere.cpp
//  PZ
//
//  Created by Philippe Devloo on 3/21/14.
//
//

#include "TPZQuadSphere.h"
#include "tpzgeomid.h"
#include "tpzgeoelmapped.h"
#include "tpzgeoelrefpattern.h"

#include "pzlog.h"

#ifdef LOG4CXX
static log4cxx::LoggerPtr logger(Logger::getLogger("pz.geom.pzgeoquad"));
#endif

TPZFMatrix<REAL> TensorProd(TPZFMatrix<REAL> &mat1, TPZFMatrix<REAL> &mat2);

namespace pzgeom {
	
	TPZGeoEl *TPZQuadSphere::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
	{
    
		int ns = orig->NSideNodes(side);
		TPZManVector<long> nodeindices(ns);
		int in;
		for(in=0; in<ns; in++)
		{
			nodeindices[in] = orig->SideNodeIndex(side,in);
		}
		long index;
		
		TPZGeoMesh *mesh = orig->Mesh();
		MElementType type = orig->Type(side);
		
		TPZGeoEl *newel = mesh->CreateGeoBlendElement(type, nodeindices, bc, index);
		TPZGeoElSide me(orig,side);
		TPZGeoElSide newelside(newel,newel->NSides()-1);
		
		newelside.InsertConnectivity(me);
		newel->Initialize();
		
		return newel;
	}
	
	void TPZQuadSphere::X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
	{
		TPZManVector<REAL,3> xqsi(3,0.); // will store (x,y,z) from (qsi,eta)
		TPZManVector<REAL,3> xqsiLxc(3,0.); // will store (x,y,z)-xc
		TPZGeoQuad::X(nodes,loc,xqsi); // gives the map (qsi,eta) to (x,y,z)

		REAL norm = 0.;
		for (int i = 0; i < 3; i++) { // Does xqsi-xc and calculates its norm
			xqsiLxc[i] = xqsi[i] - fxc[i];
			norm += xqsiLxc[i] * xqsiLxc[i];
		}
		norm = sqrt(norm);
		
		for (int i = 0; i < 3; i++) {
			result[i] = fxc[i] + xqsiLxc[i] * fR / norm;
		}

	}

	void TPZQuadSphere::Jacobian(const TPZFMatrix<REAL> & coord, TPZVec<REAL>& par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
	{
		
		// will first do dxdqsi and then d(phi*v)/dx, finaly d(phi*v)/dx * dxdqsi, where phi = 1/norm(xqsi - xc) and v = (xqsi - xc) * fR
		
		TPZManVector<REAL,3> xqsi(3,0.); // will store (x,y,z) from (qsi,eta)
		TPZManVector<REAL,3> xqsiLxc(3,0.); // will store (x,y,z)-xc
		TPZGeoQuad::X(coord,par,xqsi); // gives the map (qsi,eta) to (x,y,z)
		REAL norm = 0.;
		for (int i = 0; i < 3; i++) { // Does xqsi-xc and calculates its norm
			xqsiLxc[i] = xqsi[i] - fxc[i];
			norm += xqsiLxc[i] * xqsiLxc[i];
		}
		norm = sqrt(norm);
		
		TPZFNMatrix<6,REAL> dxdqsi(3,2,0.); // But it is a (3,2) matrix. It is set (3,3) because of the later products
		TPZGeoQuad::Jacobian(coord, par, jacobian, axes, detjac, jacinv); // first calculate the derivative dxdqsi (note a lot of dummies in the parameters)
		TPZFMatrix<REAL> axest;
		axes.Transpose(&axest);
		axest.Multiply(jacobian, dxdqsi);
		jacobian.Zero();
		

		TPZFNMatrix<3,REAL> gradphi(3,1,0.), v(3,1,0.); // here phi = 1/norm(xqsi - xc) and v = (xqsi - xc) * fR
		TPZFNMatrix<9,REAL> gradv(3,3,0.); // here v = (xqsi - xc) * fR
		REAL phi = 1./norm;
		
		for (int i = 0; i < 3; i++) {
			v(i,0) = xqsiLxc[i] * fR;
			gradv(i,i) = fR;
			gradphi(i,0) = - (1. / (norm*norm*norm) ) * xqsiLxc[i];
		}
		
		TPZFNMatrix <9,REAL> DphivDx(3,3,0.); // will store d(phi*v)/dx
		DphivDx = TensorProd(v,gradphi) + phi*gradv;
		
		TPZFNMatrix <6,REAL> GradX(3,2,0.); // stores d(phi*v)/dx * dxdqsi
		DphivDx.Multiply(dxdqsi, GradX);
		
		/*
		// Amplifying
		TPZManVector<REAL,3> minx(3,0.),maxx(3,0.);
		
		int spacedim = coord.Rows();
		
		for (int j=0; j<spacedim; j++) {
			minx[j] = coord.GetVal(j,0);
			maxx[j] = coord.GetVal(j,0);
		}
		TPZFNMatrix<6,REAL> VecMatrix(3,2,0.);
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < spacedim; j++) {
				minx[j] = minx[j] < coord.GetVal(j,i) ? minx[j]:coord.GetVal(j,i);
				maxx[j] = maxx[j] > coord.GetVal(j,i) ? maxx[j]:coord.GetVal(j,i);
			}
		}
		REAL delx = 0.;
		for (int j=0; j<spacedim; j++) {
			delx = delx > (maxx[j]-minx[j]) ? delx : (maxx[j]-minx[j]);
		}
		
		GradX *= 1./delx;
		*/
		
		GradX.Print("GradX");
		
		GradX.GramSchmidt(axest,jacobian);
		axest.Transpose(&axes);
		detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
		
		if(IsZero(detjac))
		{
#ifdef DEBUG
			std::stringstream sout;
			sout << "Singular Jacobian " << detjac;
			LOGPZ_ERROR(logger, sout.str())
#endif
			detjac = ZeroTolerance();
		}
		
		jacinv(0,0) =  jacobian(1,1)/detjac;
		jacinv(1,1) =  jacobian(0,0)/detjac;
		jacinv(0,1) = -jacobian(0,1)/detjac;
		jacinv(1,0) = -jacobian(1,0)/detjac;
		
		/*
		// Left without the amplification for learning purposes
		jacobian *= delx;
		jacinv *= 1./delx;
		detjac *= (delx*delx);
		 */
	}
	
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	/** @brief Creates a geometric element according to the type of the father element */
	TPZGeoEl *TPZQuadSphere::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
																						TPZVec<long>& nodeindexes,
																						int matid,
																						long& index)
	
	{
		return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
	}
	
	
	
}

TPZFMatrix<REAL> TensorProd(TPZFMatrix<REAL> &vec1, TPZFMatrix<REAL> &vec2)
{
	TPZFNMatrix<9,REAL> res(3,3,0.);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			res(i,j) = vec1(i,0) * vec2(j,0);
		}
	}
	return res;
}

/**
 * @ingroup geometry
 * @brief Id for three dimensional arc element
 */

template<>
int TPZGeoElRefPattern<pzgeom::TPZQuadSphere>::ClassId() const {
	return TPZGEOELEMENTQUADSPHEREID;
}

template class TPZRestoreClass< TPZGeoElRefPattern<pzgeom::TPZQuadSphere>, TPZGEOELEMENTQUADSPHEREID>;


template class TPZGeoElRefLess<pzgeom::TPZQuadSphere>;

