/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticLine methods. 
 */
#include "tpzquadraticline.h"
#include "pzshapequad.h"
#include "tpzgeoblend.h"
#include "tpzgeoelmapped.h"

#include "pzgeoelrefless.h.h"
#include "tpzgeoelrefpattern.h.h"
#include "pznoderep.h.h"

using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;

void TPZQuadraticLine::Shape(TPZVec<REAL> &param,TPZFMatrix &phi,TPZFMatrix &dphi) {
    
	REAL qsi = param[0];
	
	phi(0,0)  = -qsi*(1.-qsi)/2.;
	phi(1,0)  = +qsi*(1.+qsi)/2.;
	phi(2,0)  = (1.-qsi)*(1.+qsi);
	
	dphi(0,0) = qsi-0.5;
	dphi(0,1) = qsi+0.5;
	dphi(0,2) = -2.*qsi;
}

void TPZQuadraticLine::X(TPZFMatrix & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result) {
	TPZFNMatrix<9> phi(NNodes,1);
	TPZFNMatrix<16> dphi(1,NNodes);
	Shape(loc,phi,dphi);
	
	for(int i = 0; i < 3; i++){
		result[i] = 0.0;
		for(int j = 0; j < NNodes; j++) result[i] += phi(j,0)*coord(i,j);
	}
}

void TPZQuadraticLine::Jacobian(TPZFMatrix & coord, TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) {
#ifdef DEBUG
	if( param[0] < -1.001 || param[0] > 1.001) {
		PZError << "TPZQuadraticLine.jacobian. param out of range : "
		" param.NElements() = " << param.NElements() <<
		"\nparam[0] = " << param[0] << "\n";
	}
#endif
	
	jacobian.Resize(1,1); axes.Resize(1,3); jacinv.Resize(1,1);
	
	TPZFNMatrix<3> phi(NNodes,1);
	TPZFNMatrix<3> dphi(1,NNodes);
	Shape(param,phi,dphi);
	jacobian.Zero();
	
	TPZFNMatrix<3> VecMatrix(3,1,0.);
	for(int i = 0; i < NNodes; i++) {
		for(int j = 0; j < 3; j++) {
			VecMatrix(j,0) += coord(j,i)*dphi(0,i);
		}
	}
	
	TPZFNMatrix<9> axest;
	VecMatrix.GramSchmidt(axest,jacobian);
	axest.Transpose(&axes);
	
	detjac = jacobian(0,0);
	jacinv(0,0) =  1./jacobian(0,0);
}

TPZGeoEl *TPZQuadraticLine::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
	
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

#include "tpzgeoelmapped.h"
/**
 * Creates a geometric element according to the type of the father element
 */

TPZGeoEl *TPZQuadraticLine::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
											 TPZVec<int>& nodeindexes,
											 int matid,
											 int& index)
{
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}



///CreateGeoElement -> TPZQuadraticLine

#define TPZGEOELEMENTQUADRATICLINEID 3111
template<>
int TPZGeoElRefPattern<TPZQuadraticLine>::ClassId() const {
	return TPZGEOELEMENTQUADRATICLINEID;
}
template class 
TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticLine>, TPZGEOELEMENTQUADRATICLINEID>;



template class pzgeom::TPZNodeRep<3,TPZQuadraticLine>;
template class TPZGeoElRefLess<TPZQuadraticLine>;
