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
#include "pzcompelwithmem.h"

#include "pzbndcond.h"

#include <set>

using namespace pzgeom;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzmultiphysiccompEl"));
#endif

template <class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::TPZMultiphysicsCompEl() : TPZMultiphysicsElement(), fElementVec(0){
}

template <class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::TPZMultiphysicsCompEl(TPZCompMesh &mesh, const TPZMultiphysicsCompEl<TGeometry> &copy) : TPZMultiphysicsElement(mesh,copy),
fElementVec(copy.fElementVec), fConnectIndexes(copy.fConnectIndexes)
{
  DebugStop(); // only implemented to use withmem. Hope it is not called
}

template <class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::TPZMultiphysicsCompEl(TPZCompMesh &mesh,
                                                        const TPZMultiphysicsCompEl<TGeometry> &copy,
                                                        std::map<long,long> & gl2lcConMap,
                                                        std::map<long,long> & gl2lcElMap){
    
  DebugStop(); // if this is called, withmem should not work
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
	int side, dim, dimmf;
	nel=fElementVec.size();
	trVec.Resize(nel);
	TPZGeoEl *gelmf = Reference();
    dimmf = gelmf->Dimension();
	side = gelmf->NSides()-1;
	TPZGeoEl  *geoel;
	for (long i = 0; i<nel; i++) {
        if (!fElementVec[i]) {
            continue;
        }
		geoel = fElementVec[i].Element()->Reference();
		dim =  geoel->Dimension();
        if (dim == dimmf) {
            TPZTransform tr(dim);
            trVec[i] = gelmf->BuildTransform2(side, geoel, tr);
        }
        else
        {
            TPZTransform LocalTransf(dimmf);
            TPZGeoElSide thisgeoside(gelmf,gelmf->NSides()-1);
            TPZGeoElSide neighgeoside = fElementVec[i].Reference();
            thisgeoside.SideTransform3(neighgeoside, LocalTransf);
            TPZGeoElSide highdim(neighgeoside.Element(), neighgeoside.Element()->NSides()-1);
            trVec[i] = neighgeoside.SideToSideTransform(highdim).Multiply(LocalTransf);

        }
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
    out << "\nConnect indexes : ";
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
        
        out << "\nComputational element belonging to the mesh " << iel+1 <<":\n";
        TPZCompEl *cel = fElementVec[iel].Element();
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
            return fElementVec[el].Element()->Dimension();
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
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
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
	TPZManVector<TPZMaterialData,2> datavec;
	datavec.resize(nref);
    
	for (long iref = 0; iref<nref; iref++)
	{
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        if(!msp) continue;
        msp->InitMaterialData(datavec[iref]);
        trvec[iref].Apply(qsi, myqsi);
        datavec[iref].p = msp->MaxOrder();
        
        TPZMaterialData::MShapeFunctionType shapetype = datavec[iref].fShapeType;
        if(shapetype==datavec[iref].EVecShape){
            msp->ComputeRequiredData(datavec[iref], myqsi);
        
        }
        else
        {
            msp->ComputeShape(myqsi,datavec[iref]);
        }
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
	fConnectIndexes[inode] = index;
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
		
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
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
void TPZMultiphysicsCompEl<TGeometry>::InitializeElementMatrix(TPZElementMatrix &ef)
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
        
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        if (! msp) {
            continue;
        }
        TPZMaterial *mat = msp->Material();
        nstate += mat->NStateVariables();
        numloadcases = mat->NumLoadCases();
        
    }
    
    const int numstate = nstate;
    ef.fMat.Redim(numeq,numloadcases);
    ef.fBlock.SetNBlocks(ncon);
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
        ef.fBlock.Set(i,ndof);
    }
    ef.fConnect.Resize(ncon);
    for(i=0; i<ncon; i++){
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
            dataVec[iref].gelElId = fElementVec[iref].Element()->Reference()->Id();
        }
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        if (!msp) {
            continue;
        }
        // precisa comentar essa parte se for calcular os vetores no pontos de integracao.
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
	
	TPZManVector<REAL,3> intpointtemp(3,0.);
	REAL weight = 0.;
	
	TPZManVector<int,3> ordervec;
	//ordervec.resize(nref);
	for (long iref=0;  iref<nref; iref++)
	{
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
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
    if(intrulepoints > 1000) {
        DebugStop();
    }
	
	TPZFMatrix<REAL> jac, axe, jacInv;
	REAL detJac;
	for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
	{
		intrule->Point(int_ind,intpointtemp,weight);
		ref->Jacobian(intpointtemp, jac, axe, detJac , jacInv);
		weight *= fabs(detJac);
        for (int i = 0; i < fElementVec.size(); i++) {
            TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[i].Element());
            if (!msp) {
                continue;
            }
            datavec[i].intLocPtIndex = int_ind;
        }
        
		this->ComputeRequiredData(datavec,intpointtemp,trvec);
        
		material->Contribute(datavec,weight,ek.fMat,ef.fMat);
	}//loop over integratin points
}//CalcStiff

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::CalcResidual(TPZElementMatrix &ef)
{
    TPZMaterial * material = Material();
    if(!material){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        ef.Reset();
        return;
    }
    
    InitializeElementMatrix(ef);
    
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
    
    TPZManVector<int,3> ordervec;
    //ordervec.resize(nref);
    for (long iref=0;  iref<nref; iref++)
    {
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
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
    if(intrulepoints > 1000) {
        DebugStop();
    }
    
    TPZFMatrix<REAL> jac, axe, jacInv;
    REAL detJac;
    for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
    {
        intrule->Point(int_ind,intpointtemp,weight);
        ref->Jacobian(intpointtemp, jac, axe, detJac , jacInv);
        weight *= fabs(detJac);
        datavec[0].intLocPtIndex = int_ind;
        datavec[1].intLocPtIndex = int_ind;
        
        this->ComputeRequiredData(datavec,intpointtemp,trvec);
        
        material->Contribute(datavec,weight,ef.fMat);
    }//loop over integratin points
}//CalcStiff

/**
 * @brief Compute the integral of a variable
 */
template <class TGeometry>
TPZVec<STATE> TPZMultiphysicsCompEl<TGeometry>::IntegrateSolution(int var) const
{
    TPZManVector<STATE> result;
    TPZMaterial * material = Material();
    if(!material){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        return result;
    }
    
    if (this->NConnects() == 0) return result;//boundary discontinuous elements have this characteristic
    
    TPZMultiphysicsCompEl<TGeometry> *thisnonconst = (TPZMultiphysicsCompEl<TGeometry> *) this;
    
    TPZManVector<TPZMaterialData,3> datavec;
    const long nref = fElementVec.size();
    datavec.resize(nref);
    thisnonconst->InitMaterialData(datavec);
    
    TPZManVector<TPZTransform> trvec;
    AffineTransform(trvec);
    
    int dim = Dimension();
    TPZAutoPointer<TPZIntPoints> intrule;
    
    TPZManVector<REAL,3> intpoint(dim,0.), intpointtemp(dim,0.);
    REAL weight = 0.;
    
    TPZManVector<int,3> ordervec;
    //ordervec.resize(nref);
    for (long iref=0;  iref<nref; iref++)
    {
        datavec[iref].fNeedsSol = true;
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
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
    if(intrulepoints > 1000) {
        DebugStop();
    }
    int nvar = material->NSolutionVariables(var);
    TPZManVector<STATE> solout(var);
    result.Resize(nvar, 0.);
    
    TPZFMatrix<REAL> jac, axe, jacInv;
    REAL detJac;
    for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
    {
        intrule->Point(int_ind,intpointtemp,weight);
        ref->Jacobian(intpointtemp, jac, axe, detJac , jacInv);
        weight *= fabs(detJac);
        datavec[0].intLocPtIndex = int_ind;
        datavec[1].intLocPtIndex = int_ind;
        
        
        thisnonconst->ComputeRequiredData(datavec,intpointtemp,trvec);
        
        material->Solution(datavec, var, solout);
        
        for (int iv=0; iv<nvar; iv++) {
            result[iv] += weight*solout[iv];
        }

    }//loop over integratin points
    
    return result;
}



template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::ComputeRequiredData(TPZVec<TPZMaterialData> &datavec,
                                                        TPZVec<REAL> &intpointtemp, TPZManVector<TPZTransform> &trvec)
{
	long ElemVecSize = fElementVec.size();
	for (long iref = 0; iref < ElemVecSize; iref++)
	{
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
		if (!msp) {
			continue;
		}
        TPZManVector<REAL,3> intpoint(msp->Reference()->Dimension(),0.);
		trvec[iref].Apply(intpointtemp, intpoint);
		
		msp->ComputeRequiredData(datavec[iref], intpoint);
	}
}//ComputeRequiredData


template <class TGeometry>
const TPZIntPoints & TPZMultiphysicsCompEl<TGeometry>::GetIntegrationRule()
{
	if (fElementVec.NElements() == 0) {
		DebugStop(); // This case should be treated before
	}
	
	TPZMaterial * material = Material();
	if(!material){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
		DebugStop();
	}
	
	const long nref = fElementVec.size(); 
	TPZManVector<TPZMaterialData,3> datavec;
	datavec.resize(nref);
	this->InitMaterialData(datavec);
	int dim = Dimension();
	TPZManVector<REAL,3> intpoint(dim,0.), intpointtemp(dim,0.);
	TPZManVector<int> ordervec;
	//ordervec.resize(nref);
	for (long iref=0;  iref<nref; iref++) 
	{
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
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
	const TPZIntPoints * intrulePtr = ref->CreateSideIntegrationRule(ref->NSides()-1, order);
	const TPZIntPoints &intrule = *intrulePtr;
	return intrule;
}


template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::EvaluateError(  void (*fp)(const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv),
                           TPZVec<REAL> &errors,TPZBlock<REAL> * flux )
{
	int NErrors = this->Material()->NEvalErrors();
	errors.Resize(NErrors);
	errors.Fill(0.);
	TPZMaterial * material = this->Material();
	//TPZMaterial * matptr = material.operator->();
	if(!material){
		PZError << "TPZInterpolatedElement::EvaluateError : no material for this element\n";
		Print(PZError);
		return;
	}
	if(dynamic_cast<TPZBndCond *>(material)) {
		LOGPZ_INFO(logger,"Exiting EvaluateError - null error - boundary condition material.");
		return;
	}
	int problemdimension = Mesh()->Dimension();
	if(Reference()->Dimension() < problemdimension) return;
	
	// Adjust the order of the integration rule
	//Cesar 2007-06-27 ==>> Begin
	//this->MaxOrder is usefull to evaluate polynomial function of the aproximation space.
	//fp can be any function and max order of the integration rule could produce best results
	int dim = Dimension();
	TPZAutoPointer<TPZIntPoints> intrule = this->GetIntegrationRule().Clone();
	int maxIntOrder = intrule->GetMaxOrder();
	TPZManVector<int,3> prevorder(dim), maxorder(dim, maxIntOrder);
	//end
	intrule->GetOrder(prevorder);
	
	intrule->SetOrder(maxorder);
	
	int ndof = material->NStateVariables();
	int nflux = material->NFluxes();
	TPZManVector<STATE,10> u_exact(ndof);
	TPZFNMatrix<90,STATE> du_exact(dim+1,ndof);
	TPZManVector<REAL,10> intpoint(3), values(NErrors);
	values.Fill(0.0);
	REAL weight;
	TPZManVector<STATE,9> flux_el(nflux,0.);
	
	if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
	
	TPZManVector<TPZMaterialData,3> datavec;
	const long nref = fElementVec.size();
	datavec.resize(nref);
	InitMaterialData(datavec);
	
	TPZManVector<TPZTransform> trvec;
	AffineTransform(trvec);
	 
	int nintpoints = intrule->NPoints();
  TPZFNMatrix<9,REAL> jac, axe, jacInv;
	REAL detJac;
  TPZGeoEl *ref = this->Reference();
	
	for(int nint = 0; nint < nintpoints; nint++) {
		
		intrule->Point(nint,intpoint,weight);
    
        ref->Jacobian(intpoint, jac, axe, detJac , jacInv);
		this->ComputeRequiredData(datavec,intpoint,trvec);
    
		weight *= fabs(datavec[0].detjac);
		// this->ComputeSolution(intpoint, data.phi, data.dphix, data.axes, data.sol, data.dsol);
		//this->ComputeSolution(intpoint, data);
		//contribuicoes dos erros
		if(fp) {
			fp(datavec[0].x,u_exact,du_exact);
      material->Errors(datavec[0],u_exact,du_exact,values);
      
			for(int ier = 0; ier < NErrors; ier++)
				errors[ier] += values[ier]*weight;
		}
		
	}//fim for : integration rule
	//Norma sobre o elemento
	for(int ier = 0; ier < NErrors; ier++){
		errors[ier] = sqrt(errors[ier]);
	}//for ier
	
	intrule->SetOrder(prevorder);
}

/** Returns the maximum interpolation order of all connected elements */
template <class TGeometry>
int TPZMultiphysicsCompEl<TGeometry>::IntegrationOrder()
{
	const long nref = fElementVec.size();
    TPZVec<int> ordervec;
	ordervec.resize(nref);
	for (long iref=0;  iref<nref; iref++)
	{
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
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


/*
template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint> >;
template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear> >;
template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle> >;
template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad> >;
template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoCube> >;
template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism> >;
template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra> >;
template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid> >;
*/


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

//--------------------- WITH MEMORY ----------------------

TPZCompEl *CreateMultiphysicsPointElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
//	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint> >(mesh,gel,index) ;
//	index = -1;
//	return NULL;
}
TPZCompEl *CreateMultiphysicsLinearElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
//	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear> >(mesh,gel,index);
//	index = -1;
//	return NULL;
}
TPZCompEl *CreateMultiphysicsQuadElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
//	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad> >(mesh,gel,index);
//	index = -1;
//	return NULL;
}
TPZCompEl *CreateMultiphysicsTriangleElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
//	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle > >(mesh,gel,index);
//	index = -1;
//	return NULL;
}
TPZCompEl *CreateMultiphysicsCubeElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
//	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoCube > >(mesh,gel,index);
//	index = -1;
//	return NULL;
}
TPZCompEl *CreateMultiphysicsPrismElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
//	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism> >(mesh,gel,index);
//	index = -1;
//	return NULL;
}
TPZCompEl *CreateMultiphysicsPyramElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
//	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid > >(mesh,gel,index);
//	index = -1;
//	return NULL;
}
TPZCompEl *CreateMultiphysicsTetraElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
//	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra> >(mesh,gel,index);
//	index = -1;
//	return NULL;
}

