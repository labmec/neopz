/*
 *  pzmultiphysiccompel.cpp
 *  PZ
 *
 *  Created by Agnaldo on 9/12/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "pzmultiphysicscompel.h"

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
#include "pzmaterial.h"
#include "pzelmat.h"

#include "pzlog.h"

#include <set.h>

using namespace pzgeom;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzmultiphysiccompEl"));
#endif

template <class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::TPZMultiphysicsCompEl() : TPZMultiphysicsElement(), fElementVec(0){
}

template <class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::TPZMultiphysicsCompEl(TPZCompMesh &mesh, TPZGeoEl *ref, int &index) :TPZMultiphysicsElement(mesh, ref, index), fElementVec(0) {
}

template<class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::~TPZMultiphysicsCompEl(){	
}

//template <class TGeometry>
//void TPZMultiphysicsCompEl<TGeometry>::SetElementVec(TPZManVector<TPZCompEl *> elemVec){
//	
//	if(!elemVec){
//		LOGPZ_ERROR(logger, "TPZMultiphysicCompEl.SetElementVec called with zero pointer.");
//	}
//	
//	int dim;
//	dim=elemVec.size();
//	fElementVec.Resize(dim);
//			
//	fElementVec = elemVec;
//}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::AffineTransform(TPZManVector<TPZTransform> &trVec) const
{
	int nel;
	int side, dim;
	nel=fElementVec.size();
	trVec.Resize(nel);
	TPZGeoEl *gelmf = Reference();
	if(gelmf->Dimension() ==  1)   // whether the geometric element is TPZGeoPoint
		side = gelmf->NSides();
	else
		side = gelmf->NSides()-1;
	
	TPZGeoEl  *geoel;
	for (int i = 0; i<nel; i++) {
		geoel = fElementVec[i]->Reference();
		dim =  geoel->Dimension();
		TPZTransform tr(dim);
		trVec[i] = gelmf->BuildTransform2(side, geoel, tr);  
	}
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::GetReferenceIndexVec(TPZManVector<TPZCompMesh *> cmeshVec, std::set<int> &refIndexVec){
	
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
void TPZMultiphysicsCompEl<TGeometry>::Print(std::ostream & out) const {
	out << "Output for a computable element index: " << fIndex << std::endl;
	if(this->Reference())
	{
		out << "Center coordinate: ";
		TPZVec< REAL > centerMaster( this->Reference()->Dimension(),0. );
		TPZVec< REAL > centerEuclid( 3,0.);
		this->Reference()->CenterPoint(this->Reference()->NSides()-1,centerMaster);
		this->Reference()->X(centerMaster,centerEuclid);
		out << centerEuclid << std::endl;
	}
	if(this->Material())
	{
		out << "Material id " << this->Material()->Id() << "\n";
	}
	else {
		out << "No material\n";
	}
	
	out << "Number of connects = " << NConnects() << " Node indexes : ";
	int nod;
	for(nod=0; nod< NConnects(); nod++)
	{
		out << ConnectIndex(nod) <<  ' ' ;
	}
	
	TPZManVector<TPZTransform> tr;
	AffineTransform(tr);
	for(int ii=0;ii<tr.size();ii++)
	{
		out << std::endl; 
		out << "Transformacao para o elemento de index "<< fElementVec[ii]->Index() <<"  pertencente a malha computacional " << ii+1 << std::endl;
		out << tr[ii] << std::endl;
	}
	
	out << std::endl;
}


template <class TGeometry>
TPZCompEl * TPZMultiphysicsCompEl<TGeometry>::Clone(TPZCompMesh &mesh) const {

	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
    return 0;
}

template <class TGeometry>
TPZCompEl* TPZMultiphysicsCompEl<TGeometry>::ClonePatchEl(TPZCompMesh &mesh,
						std::map<int,int> & gl2lcConMap,
						std::map<int,int> & gl2lcElMap) const {
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
    return 0;
}

template <class TGeometry>
int TPZMultiphysicsCompEl<TGeometry>::NConnects() const {
    return fConnectIndexes.NElements();
}

template <class TGeometry>
int TPZMultiphysicsCompEl<TGeometry>::ConnectIndex(int i) const {
    return fConnectIndexes[i];
}

template <class TGeometry>
int TPZMultiphysicsCompEl<TGeometry>::Dimension() const {
	if(!fElementVec.size() || !fElementVec[0])
		return -1;
	return fElementVec[0]->Dimension();
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::ComputeSolution(TPZVec<REAL> &qsi,
				TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix &axes){
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::ComputeSolution(TPZVec<REAL> &qsi,
					 TPZVec<REAL> &normal,
					 TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
					 TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes){
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
					 const TPZFMatrix &axes, TPZVec<REAL> &sol, TPZFMatrix &dsol){
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::SetConnectIndex(int inode, int index){
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef)
{

}//void

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
	TPZAutoPointer<TPZMaterial> material = Material();
	if(!material){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
		ek.Reset();
		ef.Reset();
		return;
	}
	
	/*  {
	 std::stringstream sout;
	 sout << __PRETTY_FUNCTION__ << " material id " << material->Id();
	 LOGPZ_DEBUG(logger,sout.str());
	 }*/
	InitializeElementMatrix(ek,ef);
	
	if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
	
	//TPZMaterialData data;
//	this->InitMaterialData(data);
//	data.p = this->MaxOrder();
//	
//	int dim = Dimension();
//	TPZManVector<REAL,3> intpoint(dim,0.);
//	REAL weight = 0.;
//	
//	TPZAutoPointer<TPZIntPoints> intrule = GetIntegrationRule().Clone();
//    int order = material->IntegrationRuleOrder(data.p);
//    if(material->HasForcingFunction())
//    {
//        order = intrule->GetMaxOrder();
//    }
//    TPZManVector<int,3> intorder(dim,order);
//    intrule->SetOrder(intorder);
//	//    material->SetIntegrationRule(intrule, data.p, dim);
//	
//	int intrulepoints = intrule->NPoints();
//	for(int int_ind = 0; int_ind < intrulepoints; ++int_ind){
//		intrule->Point(int_ind,intpoint,weight);
//		this->ComputeShape(intpoint, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphix);
//		weight *= fabs(data.detjac);
//		data.intPtIndex = int_ind;
//		this->ComputeRequiredData(data, intpoint);
//		material->Contribute(data,weight,ek.fMat,ef.fMat);
//	}//loop over integratin points
	
}//CalcStiff

//---------------------------------------------------------------	
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoCube>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid>;


TPZCompEl * CreateMultiphysicsPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint>(mesh, gel, index); 
}


TPZCompEl * CreateMultiphysicsLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle >(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoCube >(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra>(mesh,gel,index);
}


TPZCompEl * CreateMultiphysicsPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid >(mesh,gel,index);
}

