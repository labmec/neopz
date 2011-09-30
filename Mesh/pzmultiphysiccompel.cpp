/*
 *  pzmultiphysiccompel.cpp
 *  PZ
 *
 *  Created by Agnaldo on 9/12/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "pzmultiphysiccompel.h"

#include "pzcompel.h"
#include "pzgeoel.h"
#include "pztrnsform.h"

#include "pzgeopoint.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "pzgeopyramid.h"


#include "pzlog.h"

#include <set.h>

using namespace pzgeom;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzmultiphysiccompEl"));
#endif

template <class TGeometry>
TPZMultiphysicCompEl<TGeometry>::TPZMultiphysicCompEl() : TPZCompEl(), fElementVec(0){
}

template <class TGeometry>
TPZMultiphysicCompEl<TGeometry>::TPZMultiphysicCompEl(TPZCompMesh &mesh, TPZGeoEl *ref, int &index) :TPZCompEl(mesh, ref, index), fElementVec(0) {
}

template<class TGeometry>
TPZMultiphysicCompEl<TGeometry>::~TPZMultiphysicCompEl(){	
}

template <class TGeometry>
void TPZMultiphysicCompEl<TGeometry>::SetElementVec(TPZManVector<TPZCompEl *> elemVec){
	
	if(!elemVec){
		LOGPZ_ERROR(logger, "TPZMultiphysicCompEl.SetElementVec called with zero pointer.");
	}
	
	int dim;
	dim=elemVec.size();
	fElementVec.Resize(dim);
			
	fElementVec = elemVec;
}

template <class TGeometry>
void TPZMultiphysicCompEl<TGeometry>::AffineTransform(TPZManVector<TPZTransform> &tr)
{
	int nel;
	nel=fElementVec.size();
	tr.Resize(nel);
	
	TPZGeoEl  *geoel;
	int side, dim;
	for (int i = 0; i<dim; i++) {
		geoel = fElementVec[i]->Reference();
		dim =  geoel->Dimension();
		side = geoel->NSides()-1;
		TPZTransform tt(dim);
		tr[i] = geoel->BuildTransform2(side, geoel, tr[i]);  
	}
}

template <class TGeometry>
void TPZMultiphysicCompEl<TGeometry>::GetReferenceIndexVec(TPZManVector<TPZCompMesh *> cmeshVec, std::set<int> &refIndexVec){
	
	if(cmeshVec.NElements() == 0) return;
	TPZCompMesh *cmesh = cmeshVec[0];
	TPZGeoMesh *gmesh = cmesh->Reference();
	gmesh->ResetReference();
	int isub;
	int ncm = cmeshVec.NElements();
	for (isub=0; isub<ncm; isub++) {
		cmeshVec[isub]->LoadReferences();
	}
	int ncel;
	TPZStack<TPZCompElSide> sidevec;
	for(int i = 0; i< ncm; i++){
		ncel = cmeshVec[i]->NElements();
		for (int j=0; j<ncel; j++) {
			TPZCompEl * cel = cmeshVec[i]->ElementVec()[j];
			if(cel){
				TPZGeoEl *geoel = cel->Reference();
				
				#ifdef DEBUG
				if (!geoel){
					PZError << "Error at " << __PRETTY_FUNCTION__ << " Geometry element null!\n";
					DebugStop();
				}
				#endif
				
				int ns = geoel->NSides();
				TPZGeoElSide *geoside = new TPZGeoElSide(geoel,ns-1);
				sidevec.Resize(0);
				geoside->HigherLevelCompElementList2(sidevec, 1,1);
				int nel = sidevec.NElements();
				if (nel==0){
					//std::cout << "Incluindo elemento " << geoel->Index() << std::endl;
					refIndexVec.insert(geoel->Index());
				}
			}
		}
	}
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Number of elements : " << refIndexVec.size() << std::endl;
		sout <<"Reference index of elements : "<< std::endl;
		set<int>::iterator it;
		for (it=refIndexVec.begin() ; it != refIndexVec.end(); it++ )
		sout << " " << *it;
		sout << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

template <class TGeometry>
TPZCompEl * TPZMultiphysicCompEl<TGeometry>::Clone(TPZCompMesh &mesh) const {

	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
    return 0;
}

template <class TGeometry>
TPZCompEl* TPZMultiphysicCompEl<TGeometry>::ClonePatchEl(TPZCompMesh &mesh,
						std::map<int,int> & gl2lcConMap,
						std::map<int,int> & gl2lcElMap) const {
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
    return 0;
}

template <class TGeometry>
int TPZMultiphysicCompEl<TGeometry>::NConnects() const {
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";

    return 0;
}

template <class TGeometry>
int TPZMultiphysicCompEl<TGeometry>::ConnectIndex(int i) const {
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
    return -1;
}

template <class TGeometry>
int TPZMultiphysicCompEl<TGeometry>::Dimension() const {
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
    return -1;
}

template <class TGeometry>
void TPZMultiphysicCompEl<TGeometry>::ComputeSolution(TPZVec<REAL> &qsi,
				TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix &axes){
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}

template <class TGeometry>
void TPZMultiphysicCompEl<TGeometry>::ComputeSolution(TPZVec<REAL> &qsi,
					 TPZVec<REAL> &normal,
					 TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
					 TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes){
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}

template <class TGeometry>
void TPZMultiphysicCompEl<TGeometry>::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
					 const TPZFMatrix &axes, TPZVec<REAL> &sol, TPZFMatrix &dsol){
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}

template <class TGeometry>
void TPZMultiphysicCompEl<TGeometry>::SetConnectIndex(int inode, int index){
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}
///---------------------------------------------------------------	
template class TPZMultiphysicCompEl<pzgeom::TPZGeoPoint>;
template class TPZMultiphysicCompEl<pzgeom::TPZGeoLinear>;
template class TPZMultiphysicCompEl<pzgeom::TPZGeoTriangle>;
template class TPZMultiphysicCompEl<pzgeom::TPZGeoQuad>;
template class TPZMultiphysicCompEl<pzgeom::TPZGeoCube>;
template class TPZMultiphysicCompEl<pzgeom::TPZGeoPrism>;
template class TPZMultiphysicCompEl<pzgeom::TPZGeoTetrahedra>;
template class TPZMultiphysicCompEl<pzgeom::TPZGeoPyramid>;


TPZCompEl * CreateMultiphysicsPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicCompEl<pzgeom::TPZGeoPoint>(mesh, gel, index); 
}


TPZCompEl * CreateMultiphysicsLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicCompEl<pzgeom::TPZGeoLinear>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicCompEl<pzgeom::TPZGeoTriangle >(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicCompEl<pzgeom::TPZGeoQuad>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicCompEl<pzgeom::TPZGeoCube >(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicCompEl<pzgeom::TPZGeoPrism>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicCompEl<pzgeom::TPZGeoTetrahedra>(mesh,gel,index);
}


TPZCompEl * CreateMultiphysicsPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicCompEl<pzgeom::TPZGeoPyramid >(mesh,gel,index);
}

