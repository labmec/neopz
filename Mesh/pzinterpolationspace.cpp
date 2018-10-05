/**
 * @file
 * @brief Contains the implementation of the TPZInterpolationSpace methods.
 */

#include "pzinterpolationspace.h"
#include "pzmaterialdata.h"
#include "pzbndcond.h"
#include "pzelmat.h"
#include "pzquad.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"
#include "pztransfer.h"
#include "tpzchangeel.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZInterpolationSpace"));
#endif

TPZInterpolationSpace::TPZInterpolationSpace()
: TPZCompEl()
{
	fPreferredOrder = -1;
}

TPZInterpolationSpace::TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy)
: TPZCompEl(mesh, copy)
{
	fPreferredOrder = copy.fPreferredOrder;
}

TPZInterpolationSpace::TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy, std::map<int64_t,int64_t> &gl2lcElMap)
: TPZCompEl(mesh, copy, gl2lcElMap)
{
	fPreferredOrder = copy.fPreferredOrder;
}

TPZInterpolationSpace::TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy, int64_t &index)
: TPZCompEl(mesh, copy, index)
{
	fPreferredOrder = copy.fPreferredOrder;
}

TPZInterpolationSpace::TPZInterpolationSpace(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index)
: TPZCompEl(mesh,gel,index)
{
	fPreferredOrder = mesh.GetDefaultOrder();
}

TPZInterpolationSpace::~TPZInterpolationSpace(){
}

int TPZInterpolationSpace::MaxOrder(){
	const int n = this->NConnects();
	int result = -1;
	int side;
	for(int i = 0; i < n; i++){
        // skip unidentified connects
        if (ConnectIndex(i) < 0) {
            continue;
        }
		side = this->Connect(i).Order();
		if (side > result) result = side;
	}//i
	return result;
}

/** @brief Adjust the integration rule according to the polynomial order of shape functions. */
void TPZInterpolationSpace::AdjustIntegrationRule()
{
    int order = MaxOrder();
    int integrationruleorder = 0;
    TPZMaterial * mat = this->Material();
    if (mat) {
        integrationruleorder = mat->IntegrationRuleOrder(order);
    }else
    {
        integrationruleorder = order + order;
    }
    SetIntegrationRule(integrationruleorder);
}

int TPZInterpolationSpace::ComputeIntegrationOrder() const {
    DebugStop();
	return 0;
}

void TPZInterpolationSpace::Print(std::ostream &out) const {
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZCompEl::Print(out);
    out << "PreferredSideOrder " << fPreferredOrder << std::endl;
}

void TPZInterpolationSpace::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                                         TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                                         REAL &detjac, TPZFMatrix<REAL> &jacinv,
                                         TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphidx){
	TPZGeoEl * ref = this->Reference();
	if (!ref){
		PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
		return;
	}//if

	ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
	this->Shape(intpoint,phi,dphi);
    this->Convert2Axes(dphi, jacinv, dphidx);
    
}

void TPZInterpolationSpace::ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data){
    
	
    this->ComputeShape(intpoint,data.x,data.jacobian,data.axes,data.detjac,data.jacinv,data.phi,data.dphi,data.dphix);
    
}

REAL TPZInterpolationSpace::InnerRadius(){
	if (!this->Reference()){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Reference() == NULL\n";
		return 0.;
	}
	return this->Reference()->ElementRadius();
}

void TPZInterpolationSpace::InitMaterialData(TPZMaterialData &data){
  data.gelElId = this->Reference()->Id();
    TPZMaterial *mat = Material();
#ifdef PZDEBUG
    if(!mat)
    {
        mat= Material();
        DebugStop();
    }
#endif
	mat->FillDataRequirements(data);
	const int dim = this->Dimension();
	const int nshape = this->NShapeF();
	const int nstate = this->Material()->NStateVariables();
    data.fShapeType = TPZMaterialData::EScalarShape;
	data.phi.Redim(nshape,1);
	data.dphi.Redim(dim,nshape);
	data.dphix.Redim(dim,nshape);
	data.axes.Redim(dim,3);
	data.jacobian.Redim(dim,dim);
	data.jacinv.Redim(dim,dim);
	data.x.Resize(3);
	if (data.fNeedsSol){
        int64_t nsol = data.sol.size();
        for (int64_t is=0; is<nsol; is++) {
            data.sol[is].Resize(nstate);
            data.dsol[is].Redim(dim,nstate);            
        }
	}
}//void

void TPZInterpolationSpace::ComputeRequiredData(TPZMaterialData &data,
                                                TPZVec<REAL> &qsi){
    data.intGlobPtIndex = -1;
    this->ComputeShape(qsi, data);
    
    if (data.fNeedsSol){
        if (data.phi.Rows()){//if shape functions are available
            this->ComputeSolution(qsi, data);
        }
        else{//if shape functions are not available
            this->ComputeSolution(qsi, data.sol, data.dsol, data.axes);
        }
    }//fNeedsSol
	
    data.x.Resize(3, 0.0);
    Reference()->X(qsi, data.x);
    data.xParametric = qsi;
    TPZManVector<REAL,3> x_center(3,0.0);
    TPZVec<REAL> center_qsi(3,0.0);
    
    if (Reference()->Dimension() == 2){
        if (Reference()->Type() == EQuadrilateral) {
            center_qsi[0] = 0.0;
            center_qsi[1] = 0.0;
        } else if (Reference()->Type() == ETriangle) {
            center_qsi[0] = 0.25;
            center_qsi[1] = 0.25;
        }
    }

    Reference()->X(center_qsi, x_center);
    data.XCenter = x_center;
    
    if (data.fNeedsHSize){
        data.HSize = 2.*this->InnerRadius();
    }//fNeedHSize
    
    if (data.fNeedsNormal){
        this->ComputeNormal(data);
    }//fNeedsNormal
}//void


void TPZInterpolationSpace::ComputeNormal(TPZMaterialData & data)
{
	data.normal.Resize(3,0.);
	
	int thisFace, neighbourFace, i, dim;
	TPZGeoEl * thisGeoEl, * neighbourGeoEl;
	TPZManVector<REAL,3> thisCenter(3,0.), neighbourCenter(3,0.), thisXVol(3,0.), neighbourXVol(3,0.), vec(3), axes1(3), axes2(3);
	
	thisGeoEl = this->Reference();
	thisFace = thisGeoEl->NSides() - 1;
    TPZGeoElSide thisside(thisGeoEl,thisFace);
    TPZCompMesh *cmesh = this->Mesh();
    int thiseldim = thisGeoEl->Dimension();
	
	TPZGeoElSide neighbourGeoElSide = thisGeoEl->Neighbour(thisFace);
    int matid = neighbourGeoElSide.Element()->MaterialId();
    while (neighbourGeoElSide != thisside) {
        matid = neighbourGeoElSide.Element()->MaterialId();
        if(cmesh->FindMaterial(matid) && neighbourGeoElSide.Element()->Dimension() > thiseldim)
        {
            break;
        }
        neighbourGeoElSide = neighbourGeoElSide.Neighbour();
    }
	neighbourGeoEl = neighbourGeoElSide.Element();
	neighbourFace = neighbourGeoEl->NSides() - 1;
    
    
	if(neighbourGeoEl == thisGeoEl)
	{
		// normal evaluation makes no sense since the internal element side doesn't have a neighbour.
		return; // place a breakpoint here if this is an issue
	}
    
	thisGeoEl->CenterPoint(thisFace, thisCenter);
	neighbourGeoEl->CenterPoint(neighbourFace, neighbourCenter);
    
	
	thisGeoEl->X(thisCenter,thisXVol);
	neighbourGeoEl->X(neighbourCenter,neighbourXVol);
	
	for(i = 0; i < 3; i++)
		vec[i] = -neighbourXVol[i] + thisXVol[i];// vector towards the center of the neighbour element
	
	dim = thisGeoEl->Dimension();
	
	switch(dim)
	{
		case(0): // normal points towards the x-direction
			data.normal[0] = 1.;
			data.normal[1] = 0.;
			data.normal[2] = 0.;
			break;
		case(1):
			for(i = 0 ; i < 3; i ++) axes1[i] = data.axes(0,i); // rib direction
			this->VectorialProd(axes1, vec, axes2);
			this->VectorialProd(axes2, axes1, data.normal, true);
			break;
		case(2):
			for(i = 0; i < 3; i++)
			{
				axes1[i] = data.axes(0,i);
				axes2[i] = data.axes(1,i);
			}
			this->VectorialProd(axes1, axes2, data.normal, true);
			break;
		case(3):// in this case the normal becomes senseless. A null vector is passed instead
			break;
		default:
			PZError << "TPZInterpolationSpace::ComputeNormal - unhandled element dimension\n";
	}
	
	// ensuring the normal vector points towards the neighbour element
	
	REAL dot = 0.;
	for(i = 0; i < 3; i++) dot += data.normal[i] * vec[i];
	
	if(dot < 0.)
		for(i = 0; i < 3; i++) data.normal[i] *= -1.;
	
}

void TPZInterpolationSpace::VectorialProd(TPZVec<REAL> & ivec, TPZVec<REAL> & jvec, TPZVec<REAL> & kvec, bool unitary)
{
	kvec.Resize(3);
	kvec[0] =  ivec[1]*jvec[2] - ivec[2]*jvec[1];
	kvec[1] = -ivec[0]*jvec[2] + ivec[2]*jvec[0];
	kvec[2] =  ivec[0]*jvec[1] - ivec[1]*jvec[0];
	
	if(unitary)
	{
		REAL size = 0.;
		int i;
		for(i = 0; i < 3; i++)size += kvec[i] * kvec[i];
		size = sqrt(size);
		//if(size <= 1.e-9)PZError << "\nTPZInterpolationSpace::VectorialProd - null result\n";
		for(i = 0; i < 3; i++)kvec[i] /= size;
	}
}

void TPZInterpolationSpace::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){
    
    TPZMaterial * material = Material();
    if(!material){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        int matid = Reference()->MaterialId();
        PZError << "Material Id which is missing = " << matid << std::endl;
        ek.Reset();
        ef.Reset();
        return;
    }
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << " material id " << material->Id();
        LOGPZ_DEBUG(logger,sout.str());
    }
#endif
    this->InitializeElementMatrix(ek,ef);
    
    if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
    
    
    TPZMaterialData data;
    //data.p = this->MaxOrder();
    this->InitMaterialData(data);
    data.p = this->MaxOrder();
    
    int dim = Dimension();
    TPZManVector<REAL,3> intpoint(dim,0.);
    REAL weight = 0.;
    
    TPZAutoPointer<TPZIntPoints> intrule = GetIntegrationRule().Clone();
    
//#ifdef PZDEBUG
//    {

        TPZManVector<int,3> intorder(dim,this->MaxOrder()*2);
        TPZAutoPointer<TPZIntPoints> intrule_clone = intrule->Clone();
        intrule_clone->SetOrder(intorder);
    
        if(intrule_clone->NPoints() > intrule->NPoints()){
            std::cout << "Element " << fIndex << " Bad integration rule points needed = " << intrule_clone->NPoints() << "; points obtained  = " <<  intrule->NPoints() << std::endl;
            intrule = intrule_clone;
        }
    
//    }
//#endif
    
//    if(material->HasForcingFunction())
//    {
//        int maxorder = intrule->GetMaxOrder();
//        if (maxorder > order) {
//            order = maxorder;
//        }
//    }
//    TPZManVector<int,3> intorder(dim,order);
//    intrule->SetOrder(intorder);
    
    int intrulepoints = intrule->NPoints();
    for(int int_ind = 0; int_ind < intrulepoints; ++int_ind){
        intrule->Point(int_ind,intpoint,weight);
        data.intLocPtIndex = int_ind;
        this->ComputeRequiredData(data, intpoint);
        weight *= fabs(data.detjac);
        
        material->Contribute(data, weight, ek.fMat, ef.fMat);
    }//loop over integratin points
    
}//CalcStiff

void TPZInterpolationSpace::CalcResidual(TPZElementMatrix &ef){
	
	TPZMaterial * material = Material();
	if(!material){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
		ef.Reset();
		return;
	}
	
	this->InitializeElementMatrix(ef);
	
	if (this->NConnects() == 0) return; //boundary discontinuous elements have this characteristic
	
	TPZMaterialData data;
	this->InitMaterialData(data);
	data.p = this->MaxOrder();
	
	int dim = Dimension();
	TPZManVector<REAL,3> intpoint(dim,0.);
	REAL weight = 0.;
	
	TPZAutoPointer<TPZIntPoints> intrule = GetIntegrationRule().Clone();
//    if(material->HasForcingFunction())
//    {
//        int maxorder = intrule->GetMaxOrder();
//        if (maxorder > order) {
//            order = maxorder;
//        }
//    }
//    TPZManVector<int,3> intorder(dim,order);
//    intrule->SetOrder(intorder);
	//  material->SetIntegrationRule(intrule, data.p, dim);
	
	int intrulepoints = intrule->NPoints();
	for(int int_ind = 0; int_ind < intrulepoints; ++int_ind){
		intrule->Point(int_ind,intpoint,weight);
		
		//this->ComputeShape(intpoint, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphix);
		
		//weight *= fabs(data.detjac);
		
		data.intLocPtIndex = int_ind;
		
		this->ComputeRequiredData(data, intpoint);
        weight *= fabs(data.detjac);
		
		material->Contribute(data,weight,ef.fMat);
		
	}//loop over integratin points
	
}//CalcResidual

void TPZInterpolationSpace::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef){
    TPZMaterial *mat = this->Material();
	const int numdof = mat->NStateVariables();
	const int nshape = this->NShapeF();
	const int ncon = this->NConnects();
    const int numloadcases = mat->NumLoadCases();
    
    ek.fMesh = Mesh();
    ek.fType = TPZElementMatrix::EK;
    ek.fOneRestraints = GetShapeRestraints();
    ef.fOneRestraints = ek.fOneRestraints;

    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
    
	ek.fBlock.SetNBlocks(ncon);
	ef.fBlock.SetNBlocks(ncon);
	ek.fNumStateVars = numdof;
	ef.fNumStateVars = numdof;
	int i;
    int numeq=0;
	for(i=0; i<ncon; i++){
        TPZConnect &c = Connect(i);
        int nshape = c.NShape();
#ifdef PZDEBUG
        if (nshape != NConnectShapeF(i,c.Order())) {
            DebugStop();
        }
#endif
        int nstate = c.NState();
        
#ifdef PZDEBUG
        if(nstate != numdof)
        {
            DebugStop();
        }
#endif
		ek.fBlock.Set(i,nshape*nstate);
		ef.fBlock.Set(i,nshape*nstate);
        numeq += nshape*nstate;
	}
	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,numloadcases);
	ek.fConnect.Resize(ncon);
	ef.fConnect.Resize(ncon);
	for(i=0; i<ncon; i++){
		(ef.fConnect)[i] = ConnectIndex(i);
		(ek.fConnect)[i] = ConnectIndex(i);
	}
}//void

void TPZInterpolationSpace::InitializeElementMatrix(TPZElementMatrix &ef){
    TPZMaterial *mat = this->Material();
	const int numdof = mat->NStateVariables();
	const int ncon = this->NConnects();
	const int nshape = this->NShapeF();
	const int numeq = nshape*numdof;
    const int numloadcases = mat->NumLoadCases();
    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
	ef.fMat.Redim(numeq,numloadcases);
	ef.fBlock.SetNBlocks(ncon);
	ef.fNumStateVars = numdof;
    ef.fOneRestraints = GetShapeRestraints();
	int i;
	for(i=0; i<ncon; i++){
        TPZConnect &c = Connect(i);
        unsigned int nshapec = c.NShape();
#ifdef PZDEBUG
        if (NConnectShapeF(i, c.Order()) != nshapec || c.NState() != numdof) {
            DebugStop();
        }
#endif
		ef.fBlock.Set(i,nshapec*numdof);
	}
	ef.fConnect.Resize(ncon);
	for(i=0; i<ncon; i++){
		(ef.fConnect)[i] = ConnectIndex(i);
	}
}//void

void TPZInterpolationSpace::Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) {
	if(var >= 100) {
		TPZCompEl::Solution(qsi,var,sol);
		return;
	}
	if(var == 99) {
		sol[0] = GetPreferredOrder();
        //        if (sol[0] != 2) {
        //            std::cout << __PRETTY_FUNCTION__ << " preferred order " << sol[0] << std::endl;
        //        }
		return;
	}
	
	TPZMaterial * material = this->Material();
	if(!material) {
		sol.Resize(0);
		return;
	}
	
    TPZMaterialData data;
    this->InitMaterialData(data);
    data.p = this->MaxOrder();
    this->ComputeShape(qsi, data);
	this->ComputeSolution(qsi,data);
    
	data.x.Resize(3);
    for (int i=0; i<qsi.size(); i++) {
        qsi[i] *= 0.999;
    }
	this->Reference()->X(qsi, data.x);
	
	int solSize = material->NSolutionVariables(var);
	sol.Resize(solSize);
	sol.Fill(0.);
	material->Solution(data, var, sol);
}

void TPZInterpolationSpace::InterpolateSolution(TPZInterpolationSpace &coarsel){
	// accumulates the transfer coefficients between the current element and the
	// coarse element into the transfer matrix, using the transformation t
	TPZMaterial * material = Material();
	if(!material) {
		PZError << __PRETTY_FUNCTION__ << " this->Material() == NULL " << std::endl;
		return;
	}
	
	TPZTransform<> t(Dimension());
	TPZGeoEl *ref = Reference();
	
	//Cedric 16/03/99
	//  Reference()->BuildTransform(NConnects(),coarsel.Reference(),t);
	t = Reference()->BuildTransform2(ref->NSides()-1,coarsel.Reference(),t);
	
	int locnod = NConnects();
	//   int cornod = coarsel.NConnects();
	int locmatsize = NShapeF();
	int cormatsize = coarsel.NShapeF();
	int nvar = material->NStateVariables();
	int dimension = Dimension();
	if (!dimension) {
		PZError << "\nExiting " << __PRETTY_FUNCTION__ << " - trying to interpolate a node solution.\n";
		return ;
	}
	
	//TPZFMatrix<REAL> loclocmat(locmatsize,locmatsize,0.);
	TPZFMatrix<STATE> loclocmat(locmatsize,locmatsize,0.);
	//TPZFMatrix<REAL> projectmat(locmatsize,nvar,0.);
	TPZFMatrix<STATE> projectmat(locmatsize,nvar,0.);
	
	TPZManVector<int,3> prevorder(dimension);
	TPZAutoPointer<TPZIntPoints> intrule = GetIntegrationRule().Clone();
	intrule->GetOrder(prevorder);
	
	int thismaxorder = this->MaxOrder();
	int coarsemaxorder = coarsel.MaxOrder();
	int maxorder = (thismaxorder > coarsemaxorder) ? thismaxorder : coarsemaxorder;
	// Cesar 2003-11-25 -->> To avoid integration warnings...
	maxorder = (2*maxorder > intrule->GetMaxOrder() ) ? intrule->GetMaxOrder() : 2*maxorder;
	TPZManVector<int,3> order(dimension,maxorder);
	for(int dim = 0; dim < dimension; dim++) {
		order[dim] = maxorder*2;
	}
	intrule->SetOrder(order);
	
	TPZFNMatrix<220> locphi(locmatsize,1);
	TPZFNMatrix<660> locdphi(dimension,locmatsize);//derivative of the shape function in the master domain
	TPZFNMatrix<660> locdphidx(dimension,locmatsize);//derivative of the shape function in the deformed domain
	
	TPZFNMatrix<220> corphi(cormatsize,1);
	TPZFNMatrix<660> cordphi(dimension,cormatsize);//derivative of the shape function in the master domain
	
	TPZManVector<REAL,3> int_point(dimension),coarse_int_point(dimension);
	TPZFNMatrix<9> jacobian(dimension,dimension),jacinv(dimension,dimension);
	TPZFNMatrix<9> axes(3,3,0.), coarseaxes(3,3,0.);
	REAL zero = 0.;
	TPZManVector<REAL,3> x(3,zero);
	//TPZManVector<TPZManVector<REAL,10>, 10> u(1);
	TPZSolVec u(1);
    u[0].resize(nvar);
	//TPZManVector<TPZFNMatrix<30>, 10> du(1);
	TPZGradSolVec du(1);
    du[0].Redim(dimension,nvar);
	
	int numintpoints = intrule->NPoints();
	REAL weight;
	int lin,ljn,cjn;
	TPZConnect *df;
	//   TPZBlock &coarseblock = coarsel.Mesh()->Block();
	
	for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {
		intrule->Point(int_ind,int_point,weight);
		REAL jac_det = 1.;
		this->ComputeShape(int_point, x, jacobian, axes, jac_det, jacinv, locphi, locdphi, locdphidx);
		weight *= jac_det;
		t.Apply(int_point,coarse_int_point);
		coarsel.ComputeSolution(coarse_int_point, u, du, coarseaxes);
		
		for(lin=0; lin<locmatsize; lin++) {
			for(ljn=0; ljn<locmatsize; ljn++) {
				loclocmat(lin,ljn) += weight*locphi(lin,0)*locphi(ljn,0);
			}
			for(cjn=0; cjn<nvar; cjn++) {
				projectmat(lin,cjn) += (STATE)weight*(STATE)locphi(lin,0)*u[0][cjn];
			}
		}
		jacobian.Zero();
	}
    REAL large = 0.;
    for (lin = 0; lin < locmatsize; lin++) {
        for (cjn = 0; cjn < locmatsize; cjn++) {
            REAL val = fabs(loclocmat(lin,cjn));
            large = large<val ? val : large;
        }
    }
    loclocmat *= 1./large;
    projectmat *= 1./large;
	
	loclocmat.SolveDirect(projectmat,ELU);
	// identify the non-zero blocks for each row
	TPZBlock<STATE> &fineblock = Mesh()->Block();
	int iv=0,in;
	for(in=0; in<locnod; in++) {
		df = &Connect(in);
		int64_t dfseq = df->SequenceNumber();
		int dfvar = fineblock.Size(dfseq);
		for(ljn=0; ljn<dfvar; ljn++) {
			fineblock(dfseq,0,ljn,0) = projectmat(iv/nvar,iv%nvar);
			iv++;
		}
	}
	intrule->SetOrder(prevorder);
}//InterpolateSolution

void TPZInterpolationSpace::CreateInterfaces(bool BetweenContinuous) {
	//nao verifica-se caso o elemento de contorno
	//eh maior em tamanho que o interface associado
	//caso AdjustBoundaryElement nao for aplicado
	//a malha eh criada consistentemente
	TPZGeoEl *ref = Reference();
	int nsides = ref->NSides();
	int InterfaceDimension = this->Material()->Dimension() - 1;
	int side;
	nsides--;//last face
	for(side=nsides;side>=0;side--) {
		if(ref->SideDimension(side) != InterfaceDimension) continue;
		TPZCompElSide thisside(this,side);
		if(this->ExistsInterface(thisside.Reference())) {
			//      cout << "TPZCompElDisc::CreateInterface inconsistent: interface already exists\n";
			continue;
		}
		TPZStack<TPZCompElSide> highlist;
		thisside.HigherLevelElementList(highlist,0,1);
		//a interface se cria uma vez so quando existem ambos
		//elementos esquerdo e direito (compu tacionais)
		if(!highlist.NElements()) {
			this->CreateInterface(side, BetweenContinuous);//s�tem iguais ou grande => pode criar a interface
		} else {
			int64_t ns = highlist.NElements();
			int64_t is;
			for(is=0; is<ns; is++) {//existem pequenos ligados ao lado atual
				const int higheldim = highlist[is].Reference().Dimension();
				if(higheldim != InterfaceDimension) continue;
				// 	TPZCompElDisc *del = dynamic_cast<TPZCompElDisc *> (highlist[is].Element());
				// 	if(!del) continue;
				
				TPZCompEl *del = highlist[is].Element();
				if(!del) continue;
				
				TPZCompElSide delside( del, highlist[is].Side() );
				TPZInterpolationSpace * delSp = dynamic_cast<TPZInterpolationSpace*>(del);
				if (!delSp) {
					PZError << "\nERROR AT " << __PRETTY_FUNCTION__ <<  " - CASE NOT AVAILABLE\n";
					return;
				}
				if ( delSp->ExistsInterface(delside.Reference()) ) {
					//          cout << "TPZCompElDisc::CreateInterface inconsistent: interface already exists\n";
				}
				else {
					delSp->CreateInterface(highlist[is].Side(), BetweenContinuous);
				}
			}
		}
	}
}

TPZInterfaceElement * TPZInterpolationSpace::CreateInterface(int side, bool BetweenContinuous)
{
	//  LOGPZ_INFO(logger, "Entering CreateInterface");
	
	TPZGeoEl *ref = Reference();
	if(!ref) {
		LOGPZ_WARN(logger, "Exiting CreateInterface Null reference reached - NULL interface returned");
		return NULL;
	}
	
	TPZCompElSide thisside(this,side);
	TPZStack<TPZCompElSide> list;
	list.Resize(0);
	thisside.EqualLevelElementList(list,0,1);//retorna distinto ao atual ou nulo
	const int64_t size = list.NElements();
    
    if (size > 1) {
        DebugStop();
    }
	//espera-se ter os elementos computacionais esquerdo e direito
	//ja criados antes de criar o elemento interface
	if(size){
		//Interface has the same material of the neighbour with lesser dimension.
		//It makes the interface have the same material of boundary conditions (TPZCompElDisc with interface dimension)
		int64_t index;
		
		TPZCompEl *list0 = list[0].Element();
		int list0side = list[0].Side();
		TPZCompElDisc * thisdisc  = dynamic_cast<TPZCompElDisc*>(this);
		TPZCompElDisc * neighdisc = dynamic_cast<TPZCompElDisc*>(list0);
		int thisside = side;
		int neighside = list0side;
		
		if (BetweenContinuous == false){
			//It means at least one element must be discontinuous
			if (!thisdisc  && !neighdisc){
				return NULL;
			}
		}
        TPZInterfaceElement * newcreatedinterface = NULL;
		
		if (Dimension() == list0->Dimension()) 
        {
            const int matid = this->Mesh()->Reference()->InterfaceMaterial(this->Material()->Id(), list0->Material()->Id() );
            TPZGeoEl *gel = ref->CreateBCGeoEl(side,matid); //isto acertou as vizinhanas da interface geometrica com o atual
            if(!gel) 
            {
#ifdef LOG4CXX
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "CreateBCGeoEl devolveu zero!@@@@";
                    LOGPZ_DEBUG(logger,sout.str());
                }
#endif
                DebugStop();
            }
            TPZCompElSide thiscompelside(this, thisside);
            TPZCompElSide neighcompelside(list0, neighside);
            int64_t index;
            newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,thiscompelside,neighcompelside);
            
        } else 
        {
            
            if(Dimension() > list0->Dimension())
            {
                const int matid = list0->Material()->Id();
                TPZGeoEl *gel = ref->CreateBCGeoEl(side,matid); //isto acertou as vizinhanas da interface geometrica com o atual
                if (!gel) {
                    DebugStop();
                }
                
                //o de volume eh o direito caso um deles seja BC
                //a normal aponta para fora do contorno
                TPZCompElSide thiscompelside(this, thisside);
                TPZCompElSide neighcompelside(list0, neighside);
                newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,thiscompelside,neighcompelside);
            } else {
                const int matid = this->Material()->Id();
                TPZGeoEl *gel = ref->CreateBCGeoEl(side,matid); //isto acertou as vizinhanas da interface geometrica com o atual
                if(!gel) DebugStop();
                //caso contrario ou caso ambos sejam de volume
                TPZCompElSide thiscompelside(this, thisside);
                TPZCompElSide neighcompelside(list0, neighside);
                newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,neighcompelside,thiscompelside);
            }
		}
		
		
		/** GeoBlend verifications ***/
#ifdef PZDEBUG
		{
			TPZGeoEl * faceGel = newcreatedinterface->Reference();
			TPZGeoEl * leftGel = newcreatedinterface->LeftElement()->Reference();
			TPZGeoEl * rightGel = newcreatedinterface->RightElement()->Reference();
			bool leftIsLinear = leftGel->IsLinearMapping();
			bool rightIsLinear = rightGel->IsLinearMapping();
			if(!leftIsLinear && !rightIsLinear){
				if(faceGel->IsGeoBlendEl() == false){
					std::cout << "\nError at " << __PRETTY_FUNCTION__ << "\n";
#ifdef LOG4CXX
					{
						std::stringstream sout;
						sout << "\nError at " << __PRETTY_FUNCTION__ << "\n";
						sout << "left gel:\n";
						leftGel->Print(sout);
						sout << "right gel:";
						rightGel->Print(sout);
						sout << "face gel:";
						faceGel->Print(sout);
						LOGPZ_ERROR(logger,sout.str());
					}
#endif
				}
			}
		}
#endif
		/** ***/
		
		return newcreatedinterface;
	}
	
	//If there is no equal level element, we try the lower elements.
	//Higher elements will not be considered by this method. In that case the interface must be created by the neighbour.
	TPZCompElSide lower = thisside.LowerLevelElementList(0);
	if(lower.Exists()){
		//Interface has the same material of the neighbour with lesser dimension.
		//It makes the interface has the same material of boundary conditions (TPZCompElDisc with interface dimension)
		int matid;
		int thisdim = this->Dimension();
		int neighbourdim = lower.Element()->Dimension();
		
		if (thisdim == neighbourdim){
			//      matid = this->Material()->Id();
			matid = this->Mesh()->Reference()->InterfaceMaterial(this->Material()->Id(), lower.Element()->Material()->Id() );
		}
		else { //one element is a boundary condition
			if (thisdim < neighbourdim) matid = this->Material()->Id();
			else
            {
                TPZCompEl *cel = lower.Element();
                if (!cel) {
                    DebugStop();
                }
                TPZMaterial *mat = cel->Material();
                if (!mat ) {
                    DebugStop();
                }
                matid = lower.Element()->Material()->Id();
            }
		}
		
		
		TPZCompEl *lowcel = lower.Element();
		int lowside = lower.Side();
		TPZCompElDisc * thisdisc  = dynamic_cast<TPZCompElDisc*>(this);
		TPZCompElDisc * neighdisc = dynamic_cast<TPZCompElDisc*>(lowcel);
		int thisside = side;
		int neighside = lowside;
		
		if (BetweenContinuous == false){
			//It means at least one element must be discontinuous
			if (!thisdisc && !neighdisc){
				return NULL;
			}
		}
        TPZInterfaceElement * newcreatedinterface = NULL;
        
		int64_t index;
		
        if(Dimension() == lowcel->Dimension()){///faces internas
            
            const int matid = this->Mesh()->Reference()->InterfaceMaterial(lowcel->Material()->Id(), this->Material()->Id() );
            TPZGeoEl *gel = ref->CreateBCGeoEl(side,matid); //isto acertou as vizinhanas da interface geometrica com o atual
            if(!gel) DebugStop();
            
            TPZCompElSide lowcelcompelside(lowcel, neighside);
            TPZCompElSide thiscompelside(this, thisside);
            int64_t index;
            newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,lowcelcompelside,thiscompelside);
        }
        else{
            
            if(Dimension() > lowcel->Dimension()){
                const int matid = lowcel->Material()->Id();
                TPZGeoEl *gel = ref->CreateBCGeoEl(side,matid);
                if (!gel) {
                    DebugStop();
                }
                
                //para que o elemento esquerdo seja de volume
                TPZCompElSide thiscompelside(this, thisside);
                TPZCompElSide lowcelcompelside(lowcel, neighside);
                newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,thiscompelside,lowcelcompelside);
            } else {
                const int matid = this->Material()->Id();
                TPZGeoEl *gel = ref->CreateBCGeoEl(side,matid);
                if (!gel) {
                    DebugStop();
                }
                TPZCompElSide thiscompelside(this, thisside);
                TPZCompElSide lowcelcompelside(lowcel, neighside);
#ifdef LOG4CXX_KEEP
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << __PRETTY_FUNCTION__ << " left element";
                    sout << lowcelcompelside << thiscompelside;
                    sout << "Left Element ";
                    lowcelcompelside.Element()->Print(sout);
                    sout << "Right Element ";
                    thiscompelside.Element()->Print(sout);
                    LOGPZ_DEBUG(logger,sout.str())
                }
#endif
                newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,lowcelcompelside,thiscompelside);
            }
        }		
		/** GeoBlend verifications ***/
#ifdef PZDEBUG
		{
			TPZGeoEl * faceGel = newcreatedinterface->Reference();
			TPZGeoEl * leftGel = newcreatedinterface->LeftElement()->Reference();
			TPZGeoEl * rightGel = newcreatedinterface->RightElement()->Reference();
			bool leftIsLinear = leftGel->IsLinearMapping();
			bool rightIsLinear = rightGel->IsLinearMapping();
			if(!leftIsLinear && !rightIsLinear){
				if(faceGel->IsGeoBlendEl() == false){
					std::cout << "\nError at " << __PRETTY_FUNCTION__ << "\n";
#ifdef LOG4CXX
					{
						std::stringstream sout;
						sout << "\nError at " << __PRETTY_FUNCTION__ << "\n";
						sout << "left gel:\n";
						leftGel->Print(sout);
						sout << "right gel:";
						rightGel->Print(sout);
						sout << "face gel:";
						faceGel->Print(sout);
						LOGPZ_ERROR(logger,sout.str());
					}
#endif
				}
			}
		}
#endif
		/** ***/
		
		return newcreatedinterface;
	}
	return NULL;
}

int TPZInterpolationSpace::ExistsInterface(TPZGeoElSide geosd){
	
	TPZGeoElSide  neighside = geosd.Neighbour();
	while(neighside.Element() && neighside.Element() != geosd.Element()){
		TPZCompElSide neighcompside = neighside.Reference();
		neighside = neighside.Neighbour();
		if(!neighcompside.Element()) continue;
		if(neighcompside.Element()->Type() == EInterface)
			return 1;
	}
	return 0;
}

void TPZInterpolationSpace::RemoveInterfaces(){
	
	int nsides = Reference()->NSides();
	if (!this->Material()){
		std::stringstream mess;
		mess << __PRETTY_FUNCTION__ << " - this->Material() == NULL, I can't RemoveInterfaces()";
		PZError << mess.str() << std::endl;
		LOGPZ_ERROR(logger, mess.str());
		return;
	}
	int InterfaceDimension = this->Material()->Dimension() - 1;
	int is;
	TPZStack<TPZCompElSide> list,equal;
	for(is=0;is<nsides;is++){
		TPZCompElSide thisside(this,is);
		if(thisside.Reference().Dimension() != InterfaceDimension) continue;
		// procurar na lista de elementos iguais
		list.Resize(0);// o lado atual �uma face
		//thisside.EqualLevelElementList(list,0,0);// monta a lista de elementos iguais
		RemoveInterface(is);// chame remove interface do elemento atual (para o side atual)
		thisside.HigherLevelElementList(list,0,0);// procurar na lista de elementos menores (todos)
		int64_t size = list.NElements(), i;            // 'isto pode incluir elementos interfaces'
		//tirando os elementos de interface da lista
		for(i=0;i<size;i++){
			if(list[i].Element()->Type() == EInterface) {
				LOGPZ_DEBUG(logger, "Removing interface element from the list of higher level elements");
				//This need to be done because otherwise list could be invalidated when an interface is removed.
				list[i] = TPZCompElSide();//tirando interface
			}
		}
		for(i=0;i<size;i++){// percorre os elementos menores
			if(!list[i].Element()) continue;
			TPZGeoElSide geolist = list[i].Reference();//TESTE
			if(geolist.Dimension() != InterfaceDimension) continue;
			equal.Resize(0);// para cada elemento menor e' preciso verificar a dimensao,
			list[i].EqualLevelElementList(equal,0,0);//montar a lista de elementos iguais (todos)
			equal.Push(list[i]);//n� �incorporado no m�odo anterior
			int64_t neq = equal.NElements(),k=-1;
			while(++k < neq) if(equal[k].Element()->Type() != EInterface) break;//procurando elemento descont�uo cujo
			if(!neq || k == neq){                               //lado faz parte da parti� do lado side do this
				LOGPZ_FATAL(logger, " Inconsistency of data");
				DebugStop();//elemento descont�uo n� achado: ERRO
			}// chame removeinterface do elemento menor
			
			TPZInterpolationSpace * equalkSp = dynamic_cast<TPZInterpolationSpace*>(equal[k].Element());
			if (!equalkSp){
				PZError << "\nERROR AT " << __PRETTY_FUNCTION__ <<  " - CASE NOT AVAILABLE\n";
				return;
			}
			equalkSp->RemoveInterface(equal[k].Side());
		}
	}
	
}

void TPZInterpolationSpace::RemoveInterface(int side) {
	
	TPZStack<TPZCompElSide> list;
	list.Resize(0);
	TPZCompElSide thisside(this,side);
	thisside.EqualLevelElementList(list,0,0);// monta a lista de elementos iguais
	int64_t size = list.NElements(),i=-1;
	while(++i < size) if(list[i].Element()->Type() == EInterface) break;// procura aquele que e derivado de TPZInterfaceEl
	if(!size || i == size){
#ifdef LOG4CXX_keep
        if (logger->isDebugEnabled())
		{
			std::stringstream sout;
			sout << __PRETTY_FUNCTION__ << " no interface element found\n";
			Print(sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		return;// nada a ser feito
	}
	// aqui existe a interface
	TPZCompEl *cel = list[i].Element();
#ifdef LOG4CXX
	TPZGeoEl *gel = cel->Reference();
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " element index " << Index() << " side " << std::endl;
		sout << "geometric element reference count " << gel->NumInterfaces();
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	delete cel;
}

void TPZInterpolationSpace::EvaluateError(std::function<void(const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv)> fp,
										  TPZVec<REAL> &errors,bool store_error){
	TPZMaterial * material = Material();
	//TPZMaterial * matptr = material.operator->();
	if (!material) {
		PZError << "TPZInterpolatedElement::EvaluateError : no material for this element\n";
		Print(PZError);
		return;
	}
	if (dynamic_cast<TPZBndCond *>(material)) {
		LOGPZ_INFO(logger, "Exiting EvaluateError - null error - boundary condition material.");
		return;
	}
	int NErrors = material->NEvalErrors();
	errors.Resize(NErrors);
	errors.Fill(0.);
	int problemdimension = Mesh()->Dimension();
	if(Reference()->Dimension() < problemdimension) return;
	
	// Adjust the order of the integration rule
	int dim = Dimension();
	TPZAutoPointer<TPZIntPoints> intrule = this->GetIntegrationRule().Clone();
	int maxIntOrder = intrule->GetMaxOrder();
    TPZManVector<int,3> prevorder(dim), maxorder(dim, maxIntOrder);
    //end
    intrule->GetOrder(prevorder);
    const int order_limit = 8;
    if(maxIntOrder > order_limit)
    {
        if (prevorder[0] > order_limit) {
            maxIntOrder = prevorder[0];
        }
        else
        {
            maxIntOrder = order_limit;
        }
    }
	
	intrule->SetOrder(maxorder);
	
	int ndof = material->NStateVariables();
	int nflux = material->NFluxes();
	TPZManVector<STATE,10> u_exact(ndof);
	TPZFNMatrix<9,STATE> du_exact(dim+1,ndof);
	TPZManVector<REAL,10> intpoint(problemdimension), values(NErrors);
	REAL weight;
	TPZManVector<STATE,9> flux_el(nflux,0.);
	
	TPZMaterialData data;
	this->InitMaterialData(data);
	int nintpoints = intrule->NPoints();
	
	for(int nint = 0; nint < nintpoints; nint++) {
		
        values.Fill(0.0);
		intrule->Point(nint,intpoint,weight);
        
        //in the case of the hdiv functions
        TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
        if(shapetype==data.EVecShape){
            this->ComputeRequiredData(data, intpoint);
        }
        else
        {
            this->ComputeShape(intpoint, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphi, data.dphix);
        }
		weight *= fabs(data.detjac);
		// this->ComputeSolution(intpoint, data.phi, data.dphix, data.axes, data.sol, data.dsol);
		//this->ComputeSolution(intpoint, data);
		//contribuicoes dos erros
        TPZGeoEl * ref = this->Reference();
        ref->X(intpoint, data.x);
		if(fp) {
			fp(data.x,u_exact,du_exact);
            
			if(data.fVecShapeIndex.NElements())
			{
				this->ComputeSolution(intpoint, data);
                				
				material->ErrorsHdiv(data,u_exact,du_exact,values);
                
			}
			else{
				this->ComputeSolution(intpoint, data.phi, data.dphix, data.axes, data.sol, data.dsol);
				material->Errors(data.x,data.sol[0],data.dsol[0],data.axes,flux_el,u_exact,du_exact,values);
			}
        
			for(int ier = 0; ier < NErrors; ier++)
				errors[ier] += weight*values[ier];
		}
		
	}//fim for : integration rule

    //Norma sobre o elemento
	for(int ier = 0; ier < NErrors; ier++){
		errors[ier] = sqrt(errors[ier]);
	}//for ier
	if(store_error)
    {
        int64_t index = Index();
        TPZFMatrix<STATE> &elvals = Mesh()->ElementSolution();
        if (elvals.Cols() < NErrors) {
            std::cout << "The element solution of the mesh should be resized before EvaluateError\n";
            DebugStop();
        }
        for (int ier=0; ier <NErrors; ier++) {
            elvals(index,ier) = errors[ier];
        }
    }
	intrule->SetOrder(prevorder);
	
}//method

void TPZInterpolationSpace::ComputeError(int errorid,
                                         TPZVec<REAL> &error){
	
	TPZMaterial * material = Material();
	if(!material){
		std::cout << "TPZCompElDisc::ComputeError : no material for this element\n";
		return;
	}
	
	TPZMaterialData data;
	this->InitMaterialData(data);
	
	REAL weight;
	int dim = Dimension();
	TPZVec<REAL> intpoint(dim,0.);
	
	TPZAutoPointer<TPZIntPoints> intrule = this->GetIntegrationRule().Clone();
	
	TPZManVector<int,3> prevorder(dim), maxorder(dim, this->MaxOrder());
	intrule->GetOrder(prevorder);
	intrule->SetOrder(maxorder);
	
	data.p = this->MaxOrder();
	data.HSize = 2.*this->InnerRadius();
	error.Fill(0.);
	int npoints = intrule->NPoints(), ip;
	for(ip=0;ip<npoints;ip++){
		intrule->Point(ip,intpoint,weight);
		this->ComputeShape(intpoint, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.phi, data.dphix);
		weight *= fabs(data.detjac);
		this->ComputeSolution(intpoint, data.phi, data.dphix, data.axes, data.sol, data.dsol);
		material->ContributeErrors(data,weight,error,errorid);
	}
	intrule->SetOrder(prevorder);
}

TPZVec<STATE> TPZInterpolationSpace::IntegrateSolution(int variable) const {
	TPZMaterial * material = Material();
    TPZVec<STATE> result;
	if(!material){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " : no material for this element\n";
		return result;
	}
	if (!this->Reference()){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " : no reference element\n";
		return result;
	}
	const int dim = this->Dimension();
    int meshdim = Mesh()->Dimension();
	REAL weight;
	TPZMaterialData data;
    TPZInterpolationSpace *thisnonconst = (TPZInterpolationSpace *) this;
    
    TPZInterpolationSpace *effective = thisnonconst;
    TPZTransform<REAL> tr(dim);
    if (dim != Mesh()->Dimension()) {
        TPZGeoElSide gelside(thisnonconst->Reference(),this->Reference()->NSides()-1);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            if (neighbour.Element()->Reference() && neighbour.Element()->Dimension() == Mesh()->Dimension()) {
                break;
            }
            neighbour = neighbour.Neighbour();
        }
        if (neighbour == gelside) {
            DebugStop();
        }
        gelside.SideTransform3(neighbour, tr);
        TPZTransform<REAL> tr2 = neighbour.Element()->SideToSideTransform(neighbour.Side(), neighbour.Element()->NSides()-1);
        tr = tr2.Multiply(tr);
        effective = dynamic_cast<TPZInterpolationSpace *> (neighbour.Element()->Reference());
        material = effective->Material();
    }
    
    if(!effective) DebugStop();
    TPZMaterialData data2d;
    thisnonconst->InitMaterialData(data2d);
	effective->InitMaterialData(data);
    data.fNeedsSol = true;
	
	TPZManVector<REAL, 3> intpoint(dim,0.);
	const int varsize = material->NSolutionVariables(variable);
    TPZManVector<STATE,3> value(varsize,0.);
	
	const TPZIntPoints &intrule = this->GetIntegrationRule();
	int npoints = intrule.NPoints(), ip, iv;
	TPZManVector<STATE> sol(varsize);
	for(ip=0;ip<npoints;ip++) {
		intrule.Point(ip,intpoint,weight);
        
        TPZManVector<REAL,3> intpointTR(meshdim,0.);
        tr.Apply(intpoint, intpointTR);
        //Tiago: Next call is performed only for computing detcaj. The previous method (Solution) has already computed jacobian.
        //       It means that the next call would not be necessary if I wrote the whole code here.
        thisnonconst->Reference()->Jacobian(intpoint, data2d.jacobian, data2d.axes, data2d.detjac, data2d.jacinv);
        effective->Reference()->Jacobian(intpointTR, data.jacobian, data.axes, data.detjac, data.jacinv);
        data.intLocPtIndex = ip;

        effective->ComputeRequiredData(data, intpointTR);
        
		sol.Fill(0.);
        material->Solution(data, variable, sol);
		weight *= fabs(data2d.detjac);
		for(iv = 0; iv < varsize; iv++) {
#ifdef STATE_COMPLEX
			value[iv] += sol[iv].real()*weight;
#else
            
            value[iv] += sol[iv]*weight;
#endif
		}//for iv
	}//for ip
    return value;
}//method


void TPZInterpolationSpace::Integrate(int variable, TPZVec<STATE> & value)//AQUIFRAN
{
    value = IntegrateSolution(variable);
}

//void TPZInterpolationSpace::IntegrateSolution(TPZVec<STATE> & value){
//	TPZMaterial * material = Material();
//	if(!material){
//		PZError << "Error at " << __PRETTY_FUNCTION__ << " : no material for this element\n";
//		return;
//	}
//	if (!this->Reference()){
//		PZError << "Error at " << __PRETTY_FUNCTION__ << " : no reference element\n";
//		return;
//	}
//	const int dim = this->Dimension();
//	REAL weight;
//	TPZMaterialData data;
//	this->InitMaterialData(data);
//	
//	TPZManVector<REAL, 3> intpoint(dim,0.);
//	const int varsize = material->NStateVariables();
//	value.Resize(varsize);
//	value.Fill(0.);
//	
//	const TPZIntPoints &intrule = this->GetIntegrationRule();
//	int npoints = intrule.NPoints(), ip, iv;
//	
//	for(ip=0;ip<npoints;ip++){
//		intrule.Point(ip,intpoint,weight);
//		data.sol.Fill(0.);
//		this->ComputeSolution(intpoint, data.sol, data.dsol, data.axes);
//		//Tiago: Next call is performet only for computing detcaj. The previous method (Solution) has already computed jacobian.
//		//       It means that the next call would not be necessary if I write the whole code here.
//		this->Reference()->Jacobian(intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
//		weight *= fabs(data.detjac);
//		for(iv = 0; iv < varsize; iv++){
//			value[iv] += data.sol[0][iv]*(STATE)weight;
//		}//for iv
//	}//for ip
//}//method

void TPZInterpolationSpace::ProjectFlux(TPZElementMatrix &ek, TPZElementMatrix &ef) {
	
	TPZMaterial * material = Material();
	if(!material){
		std::stringstream sout;
		sout << "Exiting ProjectFlux: no material for this element\n";
		Print(sout);
		LOGPZ_ERROR(logger,sout.str());
		ek.Reset();
		ef.Reset();
		return;
	}
	
	int num_flux = material->NFluxes();
	int dim = Dimension();
	int nshape = NShapeF();
	int ncon = NConnects();
	const TPZIntPoints &intrule = GetIntegrationRule();
	
	int numeq = nshape;
	ek.fMat.Resize(numeq,numeq);
	ek.fBlock.SetNBlocks(ncon);
	ef.fMat.Resize(numeq,num_flux);
	ef.fBlock.SetNBlocks(ncon);
	
	for(int i=0; i<ncon; ++i){
		(ef.fConnect)[i] = ConnectIndex(i);
		(ek.fConnect)[i] = ConnectIndex(i);
	}
	
	TPZMaterialData data;
	this->InitMaterialData(data);
	
	//TPZManVector<REAL> flux(num_flux,1);
	TPZManVector<STATE> flux(num_flux,1);
	TPZManVector<REAL,3> intpoint(dim);
	REAL weight = 0.;
	for(int int_ind = 0; int_ind < intrule.NPoints(); ++int_ind){
		
		intrule.Point(int_ind,intpoint,weight);
		this->ComputeShape(intpoint, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphi, data.dphix);
		weight *= fabs(data.detjac);
		this->ComputeSolution(intpoint, data.phi, data.dphix, data.axes, data.sol, data.dsol);
		
		material->Flux(data.x,data.sol[0],data.dsol[0],data.axes,flux);
		for(int in=0; in<nshape; in++){
			for(int ifl=0; ifl<num_flux; ifl++){
				(ef.fMat)(in,ifl) += flux[ifl]*(STATE)data.phi(in,0)*(STATE)weight;
			}//for ifl
			for(int jn = 0; jn<nshape; jn++){
				(ek.fMat)(in,jn) += data.phi(in,0)*data.phi(jn,0)*weight;
			}//for jn
		}//for in
	}//for int_ind
}//method

/** Save the element data to a stream */
void TPZInterpolationSpace::Write(TPZStream &buf, int withclassid) const
{
	TPZCompEl::Write(buf,withclassid);
	buf.Write(&fPreferredOrder,1);
}


void TPZInterpolationSpace::MinMaxSolutionValues(TPZVec<STATE> &min, TPZVec<STATE> &max){
#ifndef STATE_COMPLEX
	
	const int dim = Dimension();
	TPZManVector<REAL,3> intpoint(dim,0.);
	
	TPZAutoPointer<TPZIntPoints> intrule = GetIntegrationRule().Clone();
	TPZManVector<int,3> prevorder(dim,0);
	intrule->GetOrder(prevorder);
	
	TPZManVector<int,3> maxorder(dim,intrule->GetMaxOrder());
	intrule->SetOrder(maxorder);
	
	TPZSolVec sol;
	TPZGradSolVec dsol;
	TPZFNMatrix<9> axes(3,3,0.);
	REAL weight;
	
	int intrulepoints = intrule->NPoints();
	intrule->Point(0,intpoint,weight);
	this->ComputeSolution(intpoint, sol, dsol, axes);
	min = sol[0];
	max = sol[0];
	const int nvars = sol.NElements();
	for(int int_ind = 1; int_ind < intrulepoints; int_ind++){
		intrule->Point(int_ind,intpoint,weight);
		this->ComputeSolution(intpoint, sol, dsol, axes);
		for(int iv = 0; iv < nvars; iv++){
			if (sol[0][iv] < min[iv]) min[iv] = sol[0][iv];
			if (sol[0][iv] > max[iv]) max[iv] = sol[0][iv];
		}//iv
	}//loop over integratin points
	
	intrule->SetOrder(prevorder);
#else
	DebugStop();
#endif
	
}//void


void TPZInterpolationSpace::BuildTransferMatrix(TPZInterpolationSpace &coarsel, TPZTransform<> &t, TPZTransfer<STATE> &transfer){
	// accumulates the transfer coefficients between the current element and the
	// coarse element into the transfer matrix, using the transformation t
	TPZGeoEl *ref = Reference();
	int locnod = NConnects();
	int cornod = coarsel.NConnects();
	int locmatsize = NShapeF();
	int cormatsize = coarsel.NShapeF();
	
	// compare interpolation orders
	// the minimum interpolation order of this needs to be larger than the maximum interpolation order of coarse
	
	int mymaxorder = MaxOrder();
	int ic;
	int coarsemaxorder = this->MaxOrder();
	if(coarsemaxorder > mymaxorder) {
		std::stringstream sout;
		sout << "Exiting BuildTransferMatrix - compute the transfer matrix coarse "
        << coarsemaxorder << " me " << mymaxorder << std::endl;
		LOGPZ_ERROR(logger,sout.str());
		return;
	}
	TPZStack<int64_t> connectlistcoarse;
	TPZStack<int> corblocksize;
	TPZStack<int> dependencyordercoarse;
	connectlistcoarse.Resize(0);
	dependencyordercoarse.Resize(0);
	corblocksize.Resize(0);
	for(ic=0; ic<coarsel.NConnects(); ic++) connectlistcoarse.Push(coarsel.ConnectIndex(ic));
	coarsel.BuildConnectList(connectlistcoarse);
	TPZConnect::BuildDependencyOrder(connectlistcoarse,dependencyordercoarse,*coarsel.Mesh());
	
	// cornod = number of connects associated with the coarse element
	cornod = connectlistcoarse.NElements();
	int nvar = coarsel.Material()->NStateVariables();
	
	// number of blocks is cornod
	TPZBlock<REAL> corblock(0,cornod);
	int in;
	
	cormatsize = 0;
	for(in=0;in<cornod; in++) {
		int c = connectlistcoarse[in];
		unsigned int blsize = coarsel.Mesh()->ConnectVec()[c].NDof(*(coarsel.Mesh()))/nvar;
#ifdef PZDEBUG
        TPZConnect &con = coarsel.Mesh()->ConnectVec()[c];
        if(con.NShape() != blsize)
        {
            DebugStop();
        }
#endif
		corblock.Set(in,blsize);
		corblocksize.Push(blsize);
		cormatsize += blsize;
	}
	corblock.Resequence();
	
	//  REAL loclocmatstore[500] = {0.};
	// loclocmat is the inner product of the shape functions of the local element
	// loccormat is the inner product of the shape functions with the shape functions
	//    of the coarse element, both dependent and independent
	TPZFNMatrix<500,STATE> loclocmat(locmatsize,locmatsize);
	TPZFNMatrix<500,STATE> loccormat(locmatsize,cormatsize);
	loclocmat.Zero();
	loccormat.Zero();
	
	TPZAutoPointer<TPZIntPoints> intrule = GetIntegrationRule().Clone();
	int dimension = Dimension();
	
	TPZManVector<int> prevorder(dimension),order(dimension);
	intrule->GetOrder(prevorder);
	
	
	// compute the interpolation order of the shapefunctions squared
	int dim;
	for(dim=0; dim<dimension; dim++) {
		order[dim] = mymaxorder*2;
	}
	intrule->SetOrder(order);
	
	TPZBlock<REAL> locblock(0,locnod);
	
	for(in = 0; in < locnod; in++) {
        TPZConnect &c = Connect(in);
        unsigned int nshape = c.NShape();
#ifdef PZDEBUG
        if(NConnectShapeF(in, c.Order()) != nshape)
        {
            DebugStop();
        }
#endif
		locblock.Set(in,nshape);
	}
	locblock.Resequence();
	
	REAL locphistore[50]={0.},locdphistore[150]={0.};
	TPZFMatrix<REAL> locphi(locmatsize,1,locphistore,50);
	TPZFMatrix<REAL> locdphi(dimension,locmatsize,locdphistore,150);
	locphi.Zero();
	locdphi.Zero();
	// derivative of the shape function
	// in the master domain
	
	TPZFMatrix<REAL> corphi(cormatsize,1);
	TPZFMatrix<REAL> cordphi(dimension,cormatsize);
	// derivative of the shape function
	// in the master domain
	
	REAL jacobianstore[9],
	axesstore[9];
	TPZManVector<REAL> int_point(dimension),
	coarse_int_point(dimension);
	TPZFMatrix<REAL> jacobian(dimension,dimension,jacobianstore,9),jacinv(dimension,dimension);
	TPZFMatrix<REAL> axes(dimension,3,axesstore,9);
	TPZManVector<REAL> x(3);
	int_point.Fill(0.,0);
	REAL jac_det = 1.;
	ref->Jacobian( int_point, jacobian , axes, jac_det, jacinv);
	REAL multiplier = 1./jac_det;
	
	int numintpoints = intrule->NPoints();
	REAL weight;
	int lin,ljn,cjn;
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled() && coarsel.HasDependency()) {
        std::stringstream sout;
        coarsel.Print(sout);
        int nc = coarsel.NConnects();
        for (int ic=0; ic<nc ; ic++) {
            coarsel.Connect(ic).Print(*(coarsel.Mesh()),sout);
        }
        sout << "Connect list coarse "  << connectlistcoarse << std::endl;
        sout << "dependency order " << dependencyordercoarse << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	
	for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {
		
		intrule->Point(int_ind,int_point,weight);
		ref->Jacobian( int_point, jacobian , axes, jac_det, jacinv);
		ref->X(int_point, x);
		Shape(int_point,locphi,locdphi);
		weight *= jac_det;
		t.Apply(int_point,coarse_int_point);
		corphi.Zero();
		cordphi.Zero();
		coarsel.Shape(coarse_int_point,corphi,cordphi);
		
#ifdef LOG4CXX
        if (logger->isDebugEnabled() && coarsel.HasDependency()) {
            std::stringstream sout;
            corphi.Print("Coarse shape functions before expandShapeFunctions",sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
		coarsel.ExpandShapeFunctions(connectlistcoarse,dependencyordercoarse,corblocksize,corphi,cordphi);
#ifdef LOG4CXX
        if (logger->isDebugEnabled() && coarsel.HasDependency()) {
            std::stringstream sout;
            corphi.Print("Coarse shape functions after expandShapeFunctions",sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
		
		for(lin=0; lin<locmatsize; lin++) {
			for(ljn=0; ljn<locmatsize; ljn++) {
				loclocmat(lin,ljn) += weight*locphi(lin,0)*locphi(ljn,0)*multiplier;
			}
			for(cjn=0; cjn<cormatsize; cjn++) {
				loccormat(lin,cjn) += weight*locphi(lin,0)*corphi(cjn,0)*multiplier;
			}
		}
		jacobian.Zero();
	}
	loclocmat.SolveDirect(loccormat,ELDLt);
	
#ifdef LOG4CXX
    {
        std::stringstream sout;
        loccormat.Print("Element transfer matrix",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	
	for(in=0; in<locnod; in++) {
		//    int cind = connectlistcoarse[in];
		if(Connect(in).HasDependency()) continue;
		int64_t locblocknumber = Connect(in).SequenceNumber();
		int locblocksize = locblock.Size(in);
		int64_t locblockpos = locblock.Position(in);
		TPZStack<int> locblockvec;
		TPZStack<int> globblockvec;
		int64_t numnonzero = 0, jn;
		//      if(transfer.HasRowDefinition(locblocknumber)) continue;
		
		for(jn = 0; jn<cornod; jn++) {
			int corblocksize = corblock.Size(jn);
			int64_t corblockpos = corblock.Position(jn);
			int64_t cind = connectlistcoarse[jn];
			TPZConnect &con = coarsel.Mesh()->ConnectVec()[cind];
			if(con.HasDependency()) continue;
			int64_t corblocknumber = con.SequenceNumber();
			if(locblocksize == 0 || corblocksize == 0) continue;
			TPZFMatrix<STATE> smalll(locblocksize,corblocksize,0.);
			loccormat.GetSub(locblockpos,corblockpos, locblocksize,corblocksize,smalll);
			REAL tol = Norm(smalll);
			if(tol >= 1.e-10) {
				locblockvec.Push(jn);
				globblockvec.Push(corblocknumber);
				numnonzero++;
			}
		}
		if(transfer.HasRowDefinition(locblocknumber)) continue;
		transfer.AddBlockNumbers(locblocknumber,globblockvec);
		int64_t jnn;
		for(jnn = 0; jnn<numnonzero; jnn++) {
			jn = locblockvec[jnn];
			int corblocksize = corblock.Size(jn);
			int64_t corblockpos = corblock.Position(jn);
			if(corblocksize == 0 || locblocksize == 0) continue;
			TPZFMatrix<STATE> smalll(locblocksize,corblocksize,0.);
			loccormat.GetSub(locblockpos,corblockpos,locblocksize,corblocksize,smalll);
			transfer.SetBlockMatrix(locblocknumber,globblockvec[jnn],smalll);
		}
	}
	intrule->SetOrder(prevorder);
}


void TPZInterpolationSpace::ExpandShapeFunctions(TPZVec<int64_t> &connectlist, TPZVec<int> &dependencyorder, TPZVec<int> &blocksizes, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix) {
	int64_t numblocks =  connectlist.NElements();
	TPZCompMesh &mesh = *Mesh();
	int nhandled=0;
	int current_order = 0;
	int current_block =0;
	while(nhandled < numblocks) {
		if(dependencyorder[current_block] == current_order) {
			nhandled++;
			int64_t cind = connectlist[current_block];
			TPZConnect &con = mesh.ConnectVec()[cind];
			con.ExpandShape(cind,connectlist,blocksizes,phi,dphix);
		}
		current_block++;
		if(current_block == numblocks) {
			current_block = 0;
			current_order++;
		}
	}
}

int TPZInterpolationSpace::GetSideOrient(int side){
    
    TPZGeoEl *gel = this->Reference();
    int nsides = gel->NSides();
    int nnos = gel->NCornerNodes();
    // para elemento quadratico, ie, quadrilatero de 8 nos
    if (gel->Type()==EQuadrilateral && nnos>side) {
        int sideorient = gel->NormalOrientation(side);
        return sideorient;
    }
    
    if(side < nnos || side ==nsides-1) DebugStop();
    
    int sideorient = gel->NormalOrientation(side);
    return sideorient;
}

/**
 * @brief It set the normal orientation of the element by the side.
 * Only side that has dimension equal to my dimension minus one.
 * @param side: side of the reference elemen
 */
void TPZInterpolationSpace::SetSideOrient(int side, int sideorient){
    DebugStop();
}


/** Read the element data from a stream */
void TPZInterpolationSpace::Read(TPZStream &buf, void *context)
{
	TPZCompEl::Read(buf,context);
	buf.Read(&fPreferredOrder,1);
}

/// convert a shapefunction derivative in xi-eta to a function derivative in x,y,z
void TPZInterpolationSpace::Convert2Axes(const TPZFMatrix<REAL> &dphi, const TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &dphidx)
{
	int nshape = dphi.Cols();
    int dim = dphi.Rows();
    dphidx.Resize(dim,nshape);
	int ieq;
	switch(dim){
		case 0:
        {
            
        }
			break;
		case 1:
        {
			dphidx = dphi;
			dphidx *= jacinv.GetVal(0,0);
        }
			break;
		case 2:
        {
			for(ieq = 0; ieq < nshape; ieq++) {
				dphidx(0,ieq) = jacinv.GetVal(0,0)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,0)*dphi.GetVal(1,ieq);
				dphidx(1,ieq) = jacinv.GetVal(0,1)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,1)*dphi.GetVal(1,ieq);
			}
        }
			break;
		case 3:
        {
			for(ieq = 0; ieq < nshape; ieq++) {
				dphidx(0,ieq) = jacinv.GetVal(0,0)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,0)*dphi.GetVal(1,ieq) + jacinv.GetVal(2,0)*dphi.GetVal(2,ieq);
				dphidx(1,ieq) = jacinv.GetVal(0,1)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,1)*dphi.GetVal(1,ieq) + jacinv.GetVal(2,1)*dphi.GetVal(2,ieq);
				dphidx(2,ieq) = jacinv.GetVal(0,2)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,2)*dphi.GetVal(1,ieq) + jacinv.GetVal(2,2)*dphi.GetVal(2,ieq);
			}
        }
			break;
		default:
        {
			PZError << "Error at " << __PRETTY_FUNCTION__ << " please implement the " << dim << "d Jacobian and inverse\n";
        }
	} //switch
    
}

int TPZInterpolationSpace::ClassId() const{
    return Hash("TPZInterpolationSpace") ^ TPZCompEl::ClassId() << 1;
}
