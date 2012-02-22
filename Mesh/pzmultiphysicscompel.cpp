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
#include "pzmaterial.h"
#include "tpzautopointer.h"
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
#include "pzconnect.h"
#include "pzmaterialdata.h"
#include "pzinterpolationspace.h"
#include "pzlog.h"

#include <set>

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
//	if(gelmf->Dimension() ==  1){   // whether the geometric element is TPZGeoPoint
//		side = gelmf->NSides();
//	}
//	else{
//		side = gelmf->NSides()-1;
//	}
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
		std::set<int>::iterator it;
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
		out << "\n"<<std::endl; 
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


template<class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::Solution(TPZVec<REAL> &qsi, int var,TPZVec<REAL> &sol) {
	
	if(var >= 100) {
		TPZCompEl::Solution(qsi,var,sol);
		return;
	}
	int dim = this->Dimension();
	
	TPZAutoPointer<TPZMaterial> material = this->Material();
	if(!material){
		sol.Resize(0);
		return;
	}
	
	TPZManVector<TPZTransform> trvec;
	AffineTransform(trvec);
	
	//int neqsi= qsi.size();
	TPZVec<REAL> myqsi;
	myqsi.resize(qsi.size());
	
	const int numdof = material->NStateVariables();
	int nref = fElementVec.size();
	TPZVec<TPZMaterialData> datavec;
	datavec.resize(nref);
	
	for (int iref = 0; iref<nref; iref++)
	{		
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
		datavec[iref].p = msp->MaxOrder();
        datavec[iref].sol.resize(1);
		datavec[iref].sol[0].Resize(numdof);
		datavec[iref].sol.Fill(0.);
        datavec[iref].dsol.resize(1);
		datavec[iref].dsol[0].Redim(dim,numdof);
		datavec[iref].dsol[0].Zero();
		datavec[iref].axes.Redim(dim,3);
		datavec[iref].axes.Zero();
		
		trvec[iref].Apply(qsi, myqsi);
		msp->ComputeSolution(myqsi, datavec[iref].sol, datavec[iref].dsol, datavec[iref].axes);
		
		datavec[iref].x.Resize(3);
		msp->Reference()->X(myqsi, datavec[iref].x);
	}	
		int solSize = material->NSolutionVariables(var);
		sol.Resize(solSize);
		sol.Fill(0.);
		material->Solution(datavec, var, sol);
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::ComputeSolution(TPZVec<REAL> &qsi,
				TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix &axes){
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}//method

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::ComputeSolution(TPZVec<REAL> &qsi,
					 TPZVec<REAL> &normal,
					 TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix &leftaxes,
					 TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix &rightaxes){
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
					 const TPZFMatrix &axes, TPZSolVec &sol, TPZGradSolVec &dsol){
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
	const int ncon = this->NConnects();
	int numeq = 0;
	int ic;
		
	for(ic=0; ic<ncon; ic++)
	{
		numeq += Connect(ic).NDof(*Mesh());
	}
	
	int nref = this->fElementVec.size();
	int nstate;
	//nstate=1;
	for (int iref=0; iref<nref; iref++) {
		
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
		nstate += msp->Material()->NStateVariables();
	}
		
	const int numstate = nstate;
	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,1);
	ek.fBlock.SetNBlocks(ncon);
	ef.fBlock.SetNBlocks(ncon);
	ek.fNumStateVars = numstate;
	ef.fNumStateVars = numstate;

	int i;
	for(i=0; i<ncon; i++){
        int ndof = Connect(i).NDof(*Mesh());
#ifdef DEBUG
        TPZConnect &c = Connect(i);
        if (c.NShape()*c.NState() != ndof) {
            DebugStop();
        }
#endif
		ek.fBlock.Set(i,ndof);
		ef.fBlock.Set(i,ndof);
	}
	ek.fConnect.Resize(ncon);
	ef.fConnect.Resize(ncon);
	for(i=0; i<ncon; i++){
		(ek.fConnect)[i] = ConnectIndex(i);
		(ef.fConnect)[i] = ConnectIndex(i);
	}
	
}//void

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::InitMaterialData(TPZVec<TPZMaterialData > &dataVec)
{
	this->Material()->FillDataRequirements(dataVec);
	const int dim = this->Dimension();
	int nref = this->fElementVec.size();
	
#ifdef DEBUG
	if (nref != dataVec.size()) {
		PZError << "Error at " << __PRETTY_FUNCTION__ << " The number of materials can not be different from the size of the fElementVec !\n";
		DebugStop();
	}
#endif
	
	TPZVec<int> nshape(nref);
	for (int iref = 0; iref<nref; iref++) 
	{
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
		const int nstate = msp->Material()->NStateVariables();
		nshape[iref] =  msp->NShapeF();
		dataVec[iref].phi.Redim(nshape[iref],1);
		dataVec[iref].dphix.Redim(dim,nshape[iref]);
		dataVec[iref].axes.Redim(dim,3);
		dataVec[iref].jacobian.Redim(dim,dim);
		dataVec[iref].jacinv.Redim(dim,dim);
		dataVec[iref].x.Resize(3);
	
	
		if (dataVec[iref].fNeedsSol)
		{
            dataVec[iref].sol.resize(1);
            dataVec[iref].dsol.resize(1);
			dataVec[iref].sol[0].Resize(nstate);
			dataVec[iref].dsol[0].Redim(dim,nstate);
		}
	}
	
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
	
	InitializeElementMatrix(ek,ef);
	
	if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
	
	TPZVec<TPZMaterialData> datavec;
	const int nref = fElementVec.size(); 
	datavec.resize(nref);
	InitMaterialData(datavec);
	
	TPZManVector<TPZTransform> trvec;
	AffineTransform(trvec);
	
	int dim = Dimension();
	TPZAutoPointer<TPZIntPoints> intrule;
			
	TPZManVector<REAL,3> intpoint(dim,0.), intpointtemp(dim,0.);
	REAL weight = 0.;
	
	TPZVec<int> ordervec;
	ordervec.resize(nref);
	for (int iref=0;  iref<nref; iref++) 
	{
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
		datavec[iref].p = msp->MaxOrder();
		ordervec[iref] = datavec[iref].p; 
	}
	int order = material->IntegrationRuleOrder(ordervec);
		
	TPZGeoEl *ref = this->Reference();
	intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, order);
					
	TPZManVector<int,3> intorder(dim,order);
	intrule->SetOrder(intorder);	
	int intrulepoints = intrule->NPoints();
		
	TPZFMatrix jac, axe, jacInv;
	REAL detJac; 
	for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
	{		
		intrule->Point(int_ind,intpointtemp,weight);
		ref->Jacobian(intpointtemp, jac, axe, detJac , jacInv);
		weight *= fabs(detJac);
		for (int iref=0; iref<fElementVec.size(); iref++)
		{			
			TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
			trvec[iref].Apply(intpointtemp, intpoint);
			
			msp->ComputeShape(intpoint, datavec[iref].x, datavec[iref].jacobian, datavec[iref].axes, 
							  datavec[iref].detjac, datavec[iref].jacinv, datavec[iref].phi, datavec[iref].dphix);
			datavec[iref].intPtIndex = int_ind;
			msp->ComputeRequiredData(datavec[iref], intpoint);
		}
		material->Contribute(datavec,weight,ek.fMat,ef.fMat);
	}//loop over integratin points
	
	
}//CalcStiff


//#include "tpzpoint.h"
//#include "tpzline.h"
//#include "tpzquadrilateral.h"
//#include "tpztriangle.h"
//#include "tpzcube.h"
//#include "tpztetrahedron.h"
//#include "tpzprism.h"
//#include "tpzpyramid.h"
//#include "pzshapepoint.h"
//template<>
//void TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
//	if(dimension == 0) std::cout << "A point element has no graphical representation\n";
//}
//
//template<class TSHAPE>
//void TPZMultiphysicsCompEl<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
//	if(dimension == TSHAPE::Dimension && Material()->Id() > 0) {
//		new typename TSHAPE::GraphElType(this,&grafgrid);
//	}
//}

#include "pzgraphel.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "pzgraphel1d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pztrigraph.h"
#include "tpzgraphelt2dmapped.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel.h"
#include "pzmeshid.h"

template<class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension)
{
	TPZGeoEl *ref = Reference();
	TPZAutoPointer<TPZMaterial> material = Material();
	int mat = material->Id();
	int nsides = ref->NSides();
	
	if(dimension == 2 && mat > 0){
		if(nsides == 9){
			new TPZGraphElQ2dd(this,&grmesh);
			return;
		}
		if(nsides == 7){
			new TPZGraphElT2dMapped(this,&grmesh);
			return;
		}
	}//2d
	
	if(dimension == 3 && mat > 0){
		if(nsides == 27){
			new TPZGraphElQ3dd(this,&grmesh);
			return;
		}//cube
		if(nsides == 21){
			new TPZGraphElPrismMapped(this,&grmesh);
			return;
		}//prism
		if(nsides == 15){
			new TPZGraphElT3d(this,&grmesh);
			return;
		}//tetra
		if(nsides == 19){
			new TPZGraphElPyramidMapped(this,&grmesh);
			return;
		}//pyram
	}//3d
	
	if(dimension == 1 && mat > 0){
		new TPZGraphEl1dd(this,&grmesh);
	}//1d
}



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

