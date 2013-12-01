/**
 * @file
 * @brief Contains the implementation of the Multiphysics computational element methods.
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
//#include "pzcompelwithmem.h"

#include <set>

using namespace pzgeom;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzmultiphysiccompEl"));
#endif

template <class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::TPZMultiphysicsCompEl() : TPZMultiphysicsElement(), fElementVec(0){
}

template <class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::TPZMultiphysicsCompEl(TPZCompMesh &mesh, TPZGeoEl *ref, long &index) :TPZMultiphysicsElement(mesh, ref, index), fElementVec(0) {
}

template<class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::~TPZMultiphysicsCompEl(){	
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::AffineTransform(TPZManVector<TPZTransform> &trVec) const
{
	long nel;
	int side, dim;
	nel=fElementVec.size();
	trVec.Resize(nel);
	TPZGeoEl *gelmf = Reference();
	side = gelmf->NSides()-1;
	TPZGeoEl  *geoel;
	for (long i = 0; i<nel; i++) {
        if (!fElementVec[i]) {
            continue;
        }
		geoel = fElementVec[i]->Reference();
		dim =  geoel->Dimension();
		TPZTransform tr(dim);
		trVec[i] = gelmf->BuildTransform2(side, geoel, tr);  
	}
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::GetReferenceIndexVec(TPZManVector<TPZCompMesh *> cmeshVec, std::set<long> &refIndexVec){
	
	if(cmeshVec.NElements() == 0) return;
	TPZCompMesh *cmesh = cmeshVec[0];
	TPZGeoMesh *gmesh = cmesh->Reference();
	gmesh->ResetReference();
	long isub;
	long ncm = cmeshVec.NElements();
	for (isub=0; isub<ncm; isub++) {
		cmeshVec[isub]->LoadReferences();
	}
	long ncel;
	TPZStack<TPZCompElSide> sidevec;
	for(long i = 0; i< ncm; i++){
		ncel = cmeshVec[i]->NElements();
		for (long j=0; j<ncel; j++) {
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
				long nel = sidevec.NElements();
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
		std::set<long>::iterator it;
		for (it=refIndexVec.begin() ; it != refIndexVec.end(); it++ )
			sout << " " << *it;
		sout << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
}

/*{
    //ComputeNodElCon();
	out << "\n\t\tCOMPUTABLE GRID INFORMATIONS:\n\n";
	out << "TITLE-> " << fName << "\n\n";
	
	out << "number of connects            = " << NConnects() << std::endl;
	out << "number of elements            = " << NElements() << std::endl;
	out << "number of materials           = " << NMaterials() << std::endl;
	//  out << "number of nodal bound cond    = " << NBCConnects() << endl;
	
	out << "\n\t Connect Information:\n\n";
	int i, nelem = NConnects();
	for(i=0; i<nelem; i++) {
		if(fConnectVec[i].SequenceNumber() == -1) {
			if(fConnectVec[i].HasDependency()) {
				cout << "TPZCompMesh::Print inconsistency of connect\n";
				cout << "Index " << i << ' ';
				fConnectVec[i].Print(*this,std::cout);
			}
			continue;
		}
		out << "Index " << i << ' ';
		fConnectVec[i].Print(*this,out);
	}
	out << "\n\t Computable Element Information:\n\n";
	nelem = NElements();
	for(i=0; i<nelem; i++) {
		if(!fElementVec[i]) continue;
		TPZCompEl *el = fElementVec[i];
		out << "Index " << i << ' ';
		el->Print(out);
		if(!el->Reference()) continue;
		out << "\tReference Index = " << el->Reference()->Index() << std::endl << std::endl;
	}
	out << "\n\tMaterial Information:\n\n";
	std::map<int, TPZMaterial * >::const_iterator mit;
	nelem = NMaterials();
	for(mit=fMaterialVec.begin(); mit!= fMaterialVec.end(); mit++) {
        TPZMaterial *mat = mit->second;
        if (!mat) {
            DebugStop();
        }
		mat->Print(out);
	}

}*/

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::Print(std::ostream & out) const {
    
 	out <<"\nOutput for a computable element index: " << fIndex;
	if(this->Reference())
	{
		out << "\nCenter coordinate: ";
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
	
	out << "Number of connects = " << NConnects();
    out << "\nNode indexes : ";
	int nod;
	for(nod=0; nod< NConnects(); nod++)
	{
		out << ConnectIndex(nod) <<  ' ' ;
	}
    if(this->Reference())
        out << "\nReference Index = " << this->Reference()->Index() << std::endl;
    
    out << "\nComputational elements this multi-physics element: \n";
    long nmesh  = fElementVec.size();
    for(long iel = 0; iel< nmesh; iel++ ){
        
        out << "\nComputational element belonging to the mesh " << iel+1 <<":";
        TPZCompEl *cel = fElementVec[iel];
        if(!cel){
            out << "\n There is not element to computational mesh " << iel+1 <<"\n";
            continue;
        }
        cel->Print(out);
        if(!cel->Reference()) continue;
		out << "\tReference Index = " << cel->Reference()->Index();
        
        TPZManVector<TPZTransform> tr;
        AffineTransform(tr);
        out << "\n\tAffine transformation of the multiphysics element for this computational element:"<<"\n";
        out << std::endl;
		out <<"\t" << tr[iel];
    }
}


template <class TGeometry>
TPZCompEl * TPZMultiphysicsCompEl<TGeometry>::Clone(TPZCompMesh &mesh) const {
	
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
    return 0;
}

template <class TGeometry>
TPZCompEl* TPZMultiphysicsCompEl<TGeometry>::ClonePatchEl(TPZCompMesh &mesh,
														  std::map<long,long> & gl2lcConMap,
														  std::map<long,long> & gl2lcElMap) const {
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
    return 0;
}

template <class TGeometry>
int TPZMultiphysicsCompEl<TGeometry>::NConnects() const {
    return fConnectIndexes.NElements();
}

template <class TGeometry>
long TPZMultiphysicsCompEl<TGeometry>::ConnectIndex(int i) const {
    return fConnectIndexes[i];
}

template <class TGeometry>
int TPZMultiphysicsCompEl<TGeometry>::Dimension() const {
    for(int el = 0; el < fElementVec.NElements(); el++)
    {
        if(fElementVec[el])
        {
            return fElementVec[el]->Dimension();
        }
    }
	return -1;
}



template<class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::Integrate(int variable, TPZVec<REAL> & value){
	TPZMaterial * material = this->Material();
	if(!material){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " : no material for this element\n";
		return;
	}
	if (!this->Reference()){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " : no reference element\n";
		return;
	}

	long nref = fElementVec.size();
	TPZVec<TPZMaterialData> datavec;
	datavec.resize(nref);	
	
#ifdef DEBUG
	if (nref != datavec.size()) {
		PZError << "Error at " << __PRETTY_FUNCTION__ << " The number of materials can not be different from the size of the fElementVec !\n";
		DebugStop();
	}
#endif
	
	TPZVec<int> nshape(nref);
	for (long iref = 0; iref<nref; iref++) 
	{
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
        if(!msp) continue;
        msp->InitMaterialData(datavec[iref]);
	}	
	
	
//	TPZManVector<REAL, 3> intpoint(dim,0.);
//	const int varsize = material->NSolutionVariables(variable);
//	value.Resize(varsize);
	value.Fill(0.);
//	
//	const TPZIntPoints &intrule = this->GetIntegrationRule();
//	int npoints = intrule.NPoints(), ip, iv;
//	TPZManVector<REAL> sol(varsize);
//	for(ip=0;ip<npoints;ip++){
//		intrule.Point(ip,intpoint,weight);
//		sol.Fill(0.);
//		this->Solution(intpoint, variable, sol);
//		//Tiago: Next call is performed only for computing detcaj. The previous method (Solution) has already computed jacobian.
//		//       It means that the next call would not be necessary if I wrote the whole code here.
//		this->Reference()->Jacobian(intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
//		weight *= fabs(data.detjac);
//		for(iv = 0; iv < varsize; iv++) {
//#if !BUILD_COMPLEX_PROJECTS	
//			DebugStop();
//#else
//			value[iv] += sol[iv]*weight;
//#endif
//		}//for iv
//	}//for ip
}//method


template<class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::Solution(TPZVec<REAL> &qsi, int var,TPZVec<STATE> &sol) 
{
	
	if(var >= 100) {
		TPZCompEl::Solution(qsi,var,sol);
		return;
	}
	
	TPZMaterial * material = this->Material();
	if(!material){
		sol.Resize(0);
		return;
	}
	
	TPZManVector<TPZTransform> trvec;
	AffineTransform(trvec);
	
	TPZVec<REAL> myqsi;
	myqsi.resize(qsi.size());
	
	long nref = fElementVec.size();
	TPZVec<TPZMaterialData> datavec;
	datavec.resize(nref);
    
	for (long iref = 0; iref<nref; iref++)
	{		
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
        if(!msp) continue;
        msp->InitMaterialData(datavec[iref]);
        trvec[iref].Apply(qsi, myqsi);
        datavec[iref].p = msp->MaxOrder();
        
//        TPZMaterialData::MShapeFunctionType shapetype = datavec[iref].fShapeType;
//        if(shapetype==datavec[iref].EVecShape){
//            msp->ComputeRequiredData(datavec[iref], myqsi);
//            continue;
//        }
        
        msp->ComputeShape(myqsi,datavec[iref]);
        msp->ComputeSolution(myqsi, datavec[iref]);

		datavec[iref].x.Resize(3);
		msp->Reference()->X(myqsi, datavec[iref].x);
	}
	
	material->Solution(datavec, var, sol);
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::ComputeSolution(TPZVec<REAL> &qsi,
													   TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes){
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}//method

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::ComputeSolution(TPZVec<REAL> &qsi,
													   TPZVec<REAL> &normal,
													   TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix<REAL> &leftaxes,
													   TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes){
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
													   const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol){
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::SetConnectIndex(int inode, long index){
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
        long neqThisConn = Connect(ic).NDof(*Mesh());
		numeq += neqThisConn;
	}
	
	long nref = this->fElementVec.size();
	int nstate = 0;
	//nstate=1;
    int numloadcases = 1;
	for (long iref=0; iref<nref; iref++) {
		
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
        if (! msp) {
            continue;
        }
        TPZMaterial *mat = msp->Material();
		nstate += mat->NStateVariables();
        numloadcases = mat->NumLoadCases();
        
	}
	
	const int numstate = nstate;
	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,numloadcases);
	ek.fBlock.SetNBlocks(ncon);
	ef.fBlock.SetNBlocks(ncon);
	ek.fNumStateVars = numstate;
	ef.fNumStateVars = numstate;
	
	int i;
	for(i=0; i<ncon; i++){
        unsigned int ndof = Connect(i).NDof(*Mesh());
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
	long nref = this->fElementVec.size();
	
#ifdef DEBUG
	if (nref != dataVec.size()) {
		PZError << "Error at " << __PRETTY_FUNCTION__ << " The number of materials can not be different from the size of the fElementVec !\n";
		DebugStop();
	}
#endif
	
	TPZVec<int> nshape(nref);
	for (long iref = 0; iref < nref; iref++)
	{
        if(fElementVec[iref])
        {
            dataVec[iref].gelElId = fElementVec[iref]->Reference()->Id();
        }
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
        if (!msp) {
            continue;
        }
        msp->InitMaterialData(dataVec[iref]);
	}
    
    this->Material()->FillDataRequirements(dataVec);
	
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
	TPZMaterial * material = Material();
	if(!material){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
		ek.Reset();
		ef.Reset();
		return;
	}
	
	InitializeElementMatrix(ek,ef);
	
	if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
	
	TPZManVector<TPZMaterialData,3> datavec;
	const long nref = fElementVec.size(); 
	datavec.resize(nref);
	InitMaterialData(datavec);
	
	TPZManVector<TPZTransform> trvec;
	AffineTransform(trvec);
	
	int dim = Dimension();
	TPZAutoPointer<TPZIntPoints> intrule;
	
	TPZManVector<REAL,3> intpoint(dim,0.), intpointtemp(dim,0.);
	REAL weight = 0.;
	
	TPZManVector<int> ordervec;
	//ordervec.resize(nref);
	for (long iref=0;  iref<nref; iref++) 
	{
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
        int svec;
        if(msp)
        {
            ordervec.Resize(ordervec.size()+1);
            svec = ordervec.size();
        }
        else
        {
            continue;
        }
		datavec[iref].p = msp->MaxOrder();
		ordervec[svec-1] = datavec[iref].p;
	}
	int order = material->IntegrationRuleOrder(ordervec);
	
	TPZGeoEl *ref = this->Reference();
	intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, order);
	
	TPZManVector<int,3> intorder(dim,order);
	intrule->SetOrder(intorder);	
	int intrulepoints = intrule->NPoints();
    if(intrulepoints > 1000) {DebugStop();}
	
	TPZFMatrix<REAL> jac, axe, jacInv;
	REAL detJac; 
	for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
	{		
		intrule->Point(int_ind,intpointtemp,weight);
		ref->Jacobian(intpointtemp, jac, axe, detJac , jacInv);
		weight *= fabs(detJac);

        long ElemVecSize = fElementVec.size();
		for (long iref = 0; iref < ElemVecSize; iref++)
		{
			TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
            if (!msp) {
                continue;
            }
			trvec[iref].Apply(intpointtemp, intpoint);
			
			//msp->ComputeShape(intpoint, datavec[iref].x, datavec[iref].jacobian, datavec[iref].axes, 
            // datavec[iref].detjac, datavec[iref].jacinv, datavec[iref].phi, datavec[iref].dphix);
			datavec[iref].intLocPtIndex = int_ind;
            datavec[iref].NintPts = intrulepoints;

			msp->ComputeRequiredData(datavec[iref], intpoint);
		}
		material->Contribute(datavec,weight,ek.fMat,ef.fMat);
	}//loop over integratin points
	
	
}//CalcStiff

/** Returns the maximum interpolation order of all connected elements */
template <class TGeometry>
int TPZMultiphysicsCompEl<TGeometry>::IntegrationOrder()
{
	const long nref = fElementVec.size(); 
    TPZVec<int> ordervec;
	ordervec.resize(nref);
	for (long iref=0;  iref<nref; iref++) 
	{
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
        if(!msp) continue;
		ordervec[iref] =  msp->MaxOrder();
	}
	TPZMaterial * material = Material();
	int order = material->IntegrationRuleOrder(ordervec);
    return order;
}



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
#include "pzbndcond.h"

template<class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension)
{
	
	
	TPZGeoEl *ref = Reference();
	if (ref->Dimension() != dimension) {
		return;
	}
	TPZMaterial * material = Material();
	
	TPZBndCond * BDC = dynamic_cast<TPZBndCond * > (material);
	
	if (BDC) {
		return;
	}
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


TPZCompEl * CreateMultiphysicsPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint>(mesh, gel, index); 
}

TPZCompEl * CreateMultiphysicsLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle >(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoCube >(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid >(mesh,gel,index);
}
