/**
 * @file
 * @brief Contains the implementation of the TPZInterfaceElement methods.
 */

#include "pzelmat.h"
#include "TPZInterfaceEl.h"
#include "TPZCompElDisc.h"
#include "pzgeoelside.h"
#include "pzquad.h"
#include "TPZMaterial.h"
#include "pzconslaw.h"
#include "pzbndcond.h"
#include "pzintel.h"
#include "pzlog.h"
#include "pzinterpolationspace.h"
#include "pzmaterialdata.h"
#include "pzvec_extras.h"
#include "pzsubcmesh.h"

using namespace std;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzinterfacelement"));
static LoggerPtr logdata(Logger::getLogger("pz.material.axisymetric.data"));
#endif

void TPZInterfaceElement::SetLeftRightElements(TPZCompElSide & left, TPZCompElSide & right) {
	
	if(fLeftElSide.Element() && fRightElSide.Element()) this->DecreaseElConnected();
	
	TPZCompEl * cel = left.Element();
	if(cel) {
		TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement *>(cel);
		TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc *>(cel);
		if (!intel && !disc){
			PZError << __PRETTY_FUNCTION__ << " - Left element is not a TPZInterpolatedElement or TPZCompElDisc.\n";
			DebugStop();
		}
		
		this->fLeftElSide.SetElement( left.Element() );
		this->fLeftElSide.SetSide( left.Side() );
	}
	else {
		PZError << __PRETTY_FUNCTION__ << " - Left element is null.\n";
		DebugStop();
	}
	
	cel = right.Element();
	if (cel) {
		
		TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement *>(cel);
		TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc *>(cel);
		if (!intel && !disc){
			PZError << __PRETTY_FUNCTION__ << " - Right element is not a TPZInterpolatedElement or TPZCompElDisc.\n";
			DebugStop();
		}
		
		this->fRightElSide.SetElement( right.Element() );
		this->fRightElSide.SetSide( right.Side() );
	}
	else{
		PZError << __PRETTY_FUNCTION__ << " - Right element is null.\n";
		DebugStop();
	}
    TPZGeoEl *gel = Reference();
    if (gel->Dimension() != left.Element()->Dimension() || gel->Dimension() != right.Element()->Dimension()) {
        this->ComputeCenterNormal(fCenterNormal);
    }
    else
    {
        fCenterNormal.Resize(3, 0.);
    }

	
	this->IncrementElConnected();
    
    this->InitializeIntegrationRule();
    this->PrepareIntPtIndices();
    
}//method

void TPZInterfaceElement::DecreaseElConnected(){
	const int ncon = this->NConnects();
	for(int i = 0; i < ncon; i++){
		int index = this->ConnectIndex(i);
		fMesh->ConnectVec()[index].DecrementElConnected();
	}
}

void TPZInterfaceElement::IncrementElConnected(){
	const int ncon = this->NConnects();
	for(int i = 0; i < ncon; i++){
		int index = this->ConnectIndex(i);
		fMesh->ConnectVec()[index].IncrementElConnected();
	}
}

TPZInterfaceElement::~TPZInterfaceElement() {
    DecreaseElConnected();
    TPZGeoEl *gel = this->Reference();
    if (gel) {
        gel->ResetReference();
    }
    if (gel && gel->NumInterfaces() > 0) {
        gel->DecrementNumInterfaces();
        if (gel->NumInterfaces() == 0) {
            gel->RemoveConnectivities(); // deleta o elemento das vizinhancas
            TPZGeoMesh *gmesh = gel->Mesh();
            int index = gmesh->ElementIndex(gel); // identifica o index do elemento
            gmesh->ElementVec()[index] = NULL;
            delete gel; // deleta o elemento
        }
    }
    if(fIntegrationRule) delete fIntegrationRule;
}

TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int64_t &index,
                                         TPZCompElSide& left, TPZCompElSide& right)
: TPZRegisterClassId(&TPZInterfaceElement::ClassId),
TPZCompEl(mesh,geo,index)
{
	
	geo->SetReference(this);
	geo->IncrementNumInterfaces();
	
    if (!left.Element() || !right.Element()) {
        PZError << "Error at " << __PRETTY_FUNCTION__ << " left or right null elements\n";
    }
	if (left.Side() == -1 || right.Side() == -1){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " Side should not be -1\n";
		DebugStop();
	}
	
	this->SetLeftRightElements(left, right);
	
	this->IncrementElConnected();
}

TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int64_t &index)
: TPZRegisterClassId(&TPZInterfaceElement::ClassId),
TPZCompEl(mesh,geo,index), fLeftElSide(), fRightElSide(){
	geo->SetReference(this);
	geo->IncrementNumInterfaces();
	this->IncrementElConnected();
}

TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh, const TPZInterfaceElement &copy)
: TPZRegisterClassId(&TPZInterfaceElement::ClassId),
TPZCompEl(mesh,copy) {
	
	this->fLeftElSide.SetElement( mesh.ElementVec()[copy.fLeftElSide.Element()->Index()] );
	this->fLeftElSide.SetSide( copy.fLeftElSide.Side() );
	
	this->fRightElSide.SetElement( mesh.ElementVec()[copy.fRightElSide.Element()->Index()] );
	this->fRightElSide.SetSide( copy.fRightElSide.Side() );
	
#ifdef PZDEBUG
	if( !fLeftElSide.Element() || ! fRightElSide.Element() ) {
		cout << "Something wrong with clone of interface element\n";
		DebugStop();
	}
	if(fLeftElSide.Element()->Mesh() != &mesh || fRightElSide.Element()->Mesh() != &mesh) {
		cout << "The discontinuous elements should be cloned before the interface elements\n";
		DebugStop();
	}
#endif
	
    if (copy.fIntegrationRule) {
        fIntegrationRule = copy.fIntegrationRule->Clone();
    }
	fCenterNormal = copy.fCenterNormal;
	
	//TPZMaterial * mat = copy.Material();
	
	this->IncrementElConnected();
	
	if (this->Reference()){
		this->Reference()->IncrementNumInterfaces();
	}
	else{
		PZError << "ERROR at " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " - this->Reference() is NULL\n";
		DebugStop();
	}
}


TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,
                                         const TPZInterfaceElement &copy,
                                         std::map<int64_t,int64_t> &gl2lcConIdx,
                                         std::map<int64_t,int64_t> &gl2lcElIdx) :
TPZRegisterClassId(&TPZInterfaceElement::ClassId),TPZCompEl(mesh,copy)
{
	
	int64_t cplftIdx = copy.fLeftElSide.Element()->Index();
	int64_t cprgtIdx = copy.fRightElSide.Element()->Index();
    if (copy.fIntegrationRule) {
        fIntegrationRule = copy.fIntegrationRule->Clone();
    }
    
	if (gl2lcElIdx.find(cplftIdx) == gl2lcElIdx.end() || gl2lcElIdx.find(cprgtIdx) == gl2lcElIdx.end())
	{
		std::stringstream sout;
		sout << "ERROR in " << __PRETTY_FUNCTION__
		<< " trying to generate an interface for discontinuous elements that are not cloned."
		<< " Right idx = " << cprgtIdx << " left index = " << cplftIdx;
		LOGPZ_ERROR (logger,sout.str().c_str());
		DebugStop();
	}
	
	this->fLeftElSide.SetElement( mesh.ElementVec()[gl2lcElIdx[cplftIdx]] );
	this->fLeftElSide.SetSide( copy.fLeftElSide.Side() );
	
	this->fRightElSide.SetElement( mesh.ElementVec()[gl2lcElIdx[cprgtIdx]] );
	this->fRightElSide.SetSide( copy.fRightElSide.Side() );
	
#ifdef PZDEBUG
	if( !fLeftElSide.Element() || ! fRightElSide.Element() ) {
		cout << "Something wrong with clone of interface element\n";
		DebugStop();
	}
	if(fLeftElSide.Element()->Mesh() != &mesh || fRightElSide.Element()->Mesh() != &mesh) {
		cout << "The discontinuous elements should be cloned before the interface elements\n";
		DebugStop();
	}
#endif
	
	fCenterNormal = copy.fCenterNormal;
	//TPZMaterial * mat = copy.Material();
	this->IncrementElConnected();
	
	if (this->Reference()){
		this->Reference()->IncrementNumInterfaces();
	}
	else{
		PZError << "ERROR at " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " - this->Reference() is NULL\n";
		DebugStop();
	}
}



TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,const TPZInterfaceElement &copy,int64_t &index)
: TPZRegisterClassId(&TPZInterfaceElement::ClassId),
TPZCompEl(mesh,copy,index) {
	
	//ambos elementos esquerdo e direito j�foram clonados e moram na malha aglomerada
	//o geometrico da malha fina aponta para o computacional da malha aglomerada
	fCenterNormal = copy.fCenterNormal;
    if (copy.fIntegrationRule) {
        fIntegrationRule = copy.fIntegrationRule->Clone();
    }
	
	this->fLeftElSide.SetElement( mesh.ElementVec()[copy.fLeftElSide.Element()->Index()] );
	this->fLeftElSide.SetSide( copy.fLeftElSide.Side() );
	
	this->fRightElSide.SetElement( mesh.ElementVec()[copy.fRightElSide.Element()->Index()] );
	this->fRightElSide.SetSide( copy.fRightElSide.Side() );
	
#ifdef PZDEBUG
	if( !fLeftElSide.Element() || ! fRightElSide.Element() ) {
		cout << "TPZInterfaceElement::TPZInterfaceElement Something wrong with clone of interface element\n";
		DebugStop();
	}
	if(fLeftElSide.Element()->Mesh() != &mesh || fRightElSide.Element()->Mesh() != &mesh) {
		cout << "TPZInterfaceElement::TPZInterfaceElement The discontinuous elements should be cloned "
		<< "before the interface elements\n";
		DebugStop();
	}
#endif
	
	this->IncrementElConnected();
	
	if (this->Reference()){
		this->Reference()->IncrementNumInterfaces();
	}
	else{
		PZError << "ERROR at " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " - this->Reference() is NULL\n";
		DebugStop();
	}
}

TPZInterfaceElement::TPZInterfaceElement() : TPZRegisterClassId(&TPZInterfaceElement::ClassId),
TPZCompEl(), fLeftElSide(), fRightElSide(),
fCenterNormal(3,0.)
{
	//NOTHING TO BE DONE HERE
}

TPZCompEl * TPZInterfaceElement::CloneInterface(TPZCompMesh &aggmesh,int64_t &index, /*TPZCompElDisc **/TPZCompElSide &left, /*TPZCompElDisc **/TPZCompElSide &right) const {
	return  new TPZInterfaceElement(aggmesh, this->Reference(), index, left, right);
}

void TPZInterfaceElement::CalcResidual(TPZElementMatrix &ef){
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(Material());
	
#ifdef PZDEBUG
	if(!mat || mat->Name() == "no_name"){
		PZError << "TPZInterfaceElement::CalcResidual interface material null, do nothing\n";
		ef.Reset();
		DebugStop();
		return;
	}
#endif
	
	TPZInterpolationSpace * left = dynamic_cast<TPZInterpolationSpace*>(this->LeftElement());
	TPZInterpolationSpace * right = dynamic_cast<TPZInterpolationSpace*>(this->RightElement());
	
#ifdef PZDEBUG
	if (!left || !right){
		PZError << "\nError at TPZInterfaceElement::CalcResidual null neighbour\n";
		ef.Reset();
		DebugStop();
		return;
	}
	if(!left->Material() || !right->Material()){
		PZError << "\n Error at TPZInterfaceElement::CalcResidual null material\n";
		ef.Reset();
		DebugStop();
		return;
	}
#endif
	
    // data holds geometric data for integrating over the interface element
    // datal holds the approximation space data of the left element
    // datar holds the approximation space data of the right element
	TPZMaterialData data,datal,datar;
	const int dim = this->Dimension();
	const int diml = left->Dimension();
	const int dimr = right->Dimension();
    this->InitMaterialData(data);
	this->InitMaterialData(datal,left);
    this->InitMaterialData(datar, right);
	this->InitializeElementMatrix(ef);
	
	//LOOKING FOR MAX INTERPOLATION ORDER
	datal.p = left->MaxOrder();
	datar.p = right->MaxOrder();
	//Max interpolation order
	const int p = (datal.p > datar.p) ? datal.p : datar.p;
	
	TPZGeoEl *ref = Reference();
    int intorder = mat->IntegrationRuleOrder(p);
    
	TPZAutoPointer<TPZIntPoints> intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, intorder);
    if(mat->HasForcingFunction())
    {
        intorder = intrule->GetMaxOrder();
        TPZManVector<int,3> order(ref->Dimension(),intorder);
        intrule->SetOrder(order);
    }
	//   mat->SetIntegrationRule(intrule, p, dim);
	const int npoints = intrule->NPoints();
	
	//integration points in left and right elements: making transformations to neighbour elements
	TPZTransform<> TransfLeft, TransfRight;
	this->ComputeSideTransform(this->LeftElementSide(), TransfLeft);
	this->ComputeSideTransform(this->RightElementSide(), TransfRight);
	
	TPZManVector<REAL,3> intpoint(dim), LeftIntPoint(diml), RightIntPoint(dimr);
	REAL weight;
	//LOOP OVER INTEGRATION POINTS
	for(int ip = 0; ip < npoints; ip++){
		
		intrule->Point(ip,intpoint,weight);
		ref->Jacobian( intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
		weight *= fabs(data.detjac);
		
		this->Normal(data.axes,data.normal);
		
		TransfLeft.Apply( intpoint, LeftIntPoint );
		TransfRight.Apply( intpoint, RightIntPoint );
		
#ifdef PZDEBUG
		this->CheckConsistencyOfMappedQsi(this->LeftElementSide(), intpoint, LeftIntPoint);
		this->CheckConsistencyOfMappedQsi(this->RightElementSide(), intpoint, RightIntPoint);
#endif
		
		//left->ComputeShape(LeftIntPoint, data.x, datal.jacobian, datal.axes, datal.detjac, datal.jacinv, datal.phi, datal.dphix);
		//right->ComputeShape(RightIntPoint, data.x, datar.jacobian, datar.axes, datar.detjac, datar.jacinv, datar.phi, datar.dphix);
		
		this->ComputeRequiredData(datal, left, LeftIntPoint);
		this->ComputeRequiredData(datar, right, RightIntPoint);
        data.intLocPtIndex = ip;
        this->ComputeRequiredData(data,intpoint);
		mat->ContributeInterface(data, datal, datar, weight, ef.fMat);
		
	}//loop over integration points
	
}

int TPZInterfaceElement::NConnects() const {
	return this->NLeftConnects() + this->NRightConnects();
}

int TPZInterfaceElement::NLeftConnects() const{
	TPZCompEl * LeftEl  = fLeftElSide.Element();
	if (!LeftEl) return 0;
	return LeftEl->NConnects();
}

int TPZInterfaceElement::NRightConnects() const{
	TPZCompEl * RightEl = fRightElSide.Element();
	if (!RightEl) return 0;
	return RightEl->NConnects();
}

int TPZInterfaceElement::ComputeIntegrationOrder() const {
    
    TPZInterpolationSpace * left = dynamic_cast<TPZInterpolationSpace*>(this->LeftElement());
    TPZInterpolationSpace * right = dynamic_cast<TPZInterpolationSpace*>(this->RightElement());
    if (!left || !right) return -1;
    
    int leftmaxp = left->MaxOrder();
    int rightmaxp = right->MaxOrder();
    int p = (leftmaxp > rightmaxp) ? leftmaxp : rightmaxp;
    return 2*p;
}

int64_t TPZInterfaceElement::ConnectIndex(int i) const {
	
	const int nleftcon = this->NLeftConnects();
	const int nrightcon = this->NRightConnects();
	const int ncon = nleftcon + nrightcon;
	
	if(i < 0 || i >= ncon){
		PZError << "TPZInterfaceElement::ConnectIndex wrong argument i, i = " << i << endl;
		DebugStop();
		return -1;
	}
	
	if(i < nleftcon){ //required connect is associated to left neighbour
        TPZCompMesh *leftmesh = fLeftElSide.Element()->Mesh();
        if (leftmesh == Mesh()) {
            return fLeftElSide.Element()->ConnectIndex(i);
        }
        else
        {
            int64_t leftindex = fLeftElSide.Element()->ConnectIndex(i);
            TPZCompMesh *comm = Mesh()->CommonMesh(leftmesh);
            int64_t superind = leftmesh->PutinSuperMesh(leftindex, comm);
            int64_t hereindex = Mesh()->GetFromSuperMesh(superind, comm);
            return hereindex;
        }
	}
	
	if(i < ncon){ //required connect is associated to right neighbour
        TPZCompMesh *rightmesh = fRightElSide.Element()->Mesh();
        if (rightmesh == Mesh()) {
            return fRightElSide.Element()->ConnectIndex(i-nleftcon);
        }
        else
        {
            int64_t rightindex = fRightElSide.Element()->ConnectIndex(i-nleftcon);
            TPZCompMesh *comm = Mesh()->CommonMesh(rightmesh);
            int64_t superind = rightmesh->PutinSuperMesh(rightindex, comm);
            int64_t hereindex = Mesh()->GetFromSuperMesh(superind, comm);
            return hereindex;
        }
	}
	return -1;
}

void TPZInterfaceElement::Print(std::ostream &out) const {
	
	TPZCompEl::Print(out);
	out << "\nInterface element : \n";
	if(!LeftElement() || !LeftElement()->Reference()) out << "\tNULL LeftElement - this is inconsistent\n";
	else{
		out << "\tLeft CompEl Index: " << LeftElement()->Index() << endl;
		out << "\tLeft Geometric Index: " << LeftElement()->Reference()->Index() << endl;
		out << "\tLeft Geometric Id: " << LeftElement()->Reference()->Id() << endl;
		out << "\tElement Dimension " << LeftElement()->Reference()->Dimension() << endl;
	}
	
	if(!RightElement() || !RightElement()->Reference()) out << "\tNULL RightElement - this is inconsistent";
	else{
		out << "\tRight CompEl Index: " << RightElement()->Index() << endl;
		out << "\tRight Geometric Index: " << RightElement()->Reference()->Index() << endl;
		out << "\tRight Geometric Id: " << RightElement()->Reference()->Id() << endl;
		out << "\tElement Dimension " << RightElement()->Reference()->Dimension() << endl;
	}
	
	out << "\tMaterial id : " << Reference()->MaterialId() << endl;
	
	out << "\tNormal vector (at center point): ";
	out << "(" << fCenterNormal[0] << "," << fCenterNormal[1] << "," << fCenterNormal[2] << ")\n";
}

void TPZInterfaceElement::SetConnectIndex(int node, int64_t index) {
	cout << "TPZInterfaceElement::SetConnectIndex should never be called\n";
	DebugStop();
}

int TPZInterfaceElement::main(TPZCompMesh &cmesh){
	// esta func� testa o correto desempenho do algoritmo que cria e
	// deleta elementos de interface numa malha sujeita a refinamento h
	
	// InterfaceDimension �a dimens� do elemento de interface
	// verifica-se para cada lado de dimens� InterfaceDimension do
	// elemento que existe um elemento interface e que este �nico
	
	int64_t iel,iside,nel = cmesh.NElements();
	
	int InterfaceDimension;
	
	for(iel=0;iel<nel;iel++){
		TPZCompEl *cel = cmesh.ElementVec()[iel];
		if(!cel) continue;
		TPZGeoEl *geo = cel->Reference();
		InterfaceDimension = cel->Material()->Dimension() -1;
		if(!geo){
			PZError << "TPZInterfaceElement::main computational element with null reference\n";
			DebugStop();
		}
		int nsides = geo->NSides();
		for(iside=0;iside<nsides;iside++){
			if(geo->SideDimension(iside) != InterfaceDimension) continue;
			TPZCompElSide compside(cel,iside);
			if(ExistInterfaces(compside)){
				continue;
			} else {
				PZError << "TPZInterfaceEl::main interface error\t->\t";
				int nint = ExistInterfaces(compside);
				PZError << "number of existing interfaces : " << nint << endl;
				return 0;
			}
		}
	}//fim for iel
	if(!FreeInterface(cmesh)) return 0;
	return 1;
}

int TPZInterfaceElement::ExistInterfaces(TPZCompElSide &comp){
	
	TPZStack<TPZCompElSide> list;
	list.Resize(0);
	
	if(!comp.Exists()){
		PZError << "TPZInterfaceElement::ExistInterfaces null argument, do nothing it verify\n";
		return 1;//sem problemas
	}
	comp.HigherLevelElementList(list,0,0);
	int64_t cap = list.NElements();
	
	if(cap) {
		//caso existem elementos pequenos n� deve existir
		//interface associada ao lado atual, o lado atual
		//deve apontar para elemento computacional nulo
		TPZGeoElSide geo = comp.Reference(),neigh;
		neigh = geo.Neighbour();
		while(neigh.Exists() && neigh != geo){
			if(neigh.Element()->Reference()){
				PZError << "TPZInterfaceElement::ExistInterfaces error of data structure\n";
				DebugStop();
			}
			neigh = neigh.Neighbour();
		}
		//caso o vizinho n� existe todo bem
		//caso existe n� pode ter refer�cia computacional
		return 1;//sem problemas
	}
	
	//neste estagio o lado atual enxerga um elemento vizinho ou
	//est�comtido no lado de um elemento maior, portanto deve
	//ter associado um elemento interface
	TPZGeoElSide geo = comp.Reference();
	if(!geo.Exists()){
		PZError << "TPZInterfaceElement::ExistInterfaces error of data structure\n";
		DebugStop();
	}
	TPZGeoElSide  neigh = geo.Neighbour();
	int exists = 0;
	if(comp.Element()->Type() == EInterface) exists++;//o pr�rio �interface
	
	while(neigh.Element() && neigh.Element() != geo.Element()){
		TPZCompElSide comp = neigh.Reference();
		neigh = neigh.Neighbour();
		if(!comp.Element()) continue;
		if(comp.Element()->Type() == EInterface) exists++;
	}
	if(exists != 1) return 0;
	return 1;//existe uma nica interface
}

int TPZInterfaceElement::FreeInterface(TPZCompMesh &cmesh){
	
	int64_t iel,nel = cmesh.NElements();
	for(iel=0;iel<nel;iel++){
		TPZCompEl *cel = cmesh.ElementVec()[iel];
		if(!cel) continue;
		if(cel->Type() != EInterface) continue;//interessa s�interfaces
		TPZGeoEl *gel = cel->Reference();
		if(!gel){
			PZError << "TPZInterfaceElement::FreeInterface computational element with null reference\n";
			DebugStop();
            continue;
		}
		int nsides = gel->NSides();
		TPZCompElSide compside(cel,nsides-1);//face ou aresta
		TPZGeoElSide geo = compside.Reference();
		TPZGeoElSide neigh = geo.Neighbour();
		int exists = 0;
		while(neigh.Element() && neigh.Element() != geo.Element()){
			TPZCompElSide comp = neigh.Reference();
			neigh = neigh.Neighbour();
			if(!comp.Element()) continue;
			if(comp.Element()->Type() != EInterface) exists++;
		}
		//s�pode haver 1 ou 2 elementos de volume associados a um el. interface
		if(exists < 1 || exists > 2) return 0;
	}
	return 1;
}

void TPZInterfaceElement::ComputeCenterNormal(TPZVec<REAL> &normal){
	TPZManVector<REAL> qsi(3);
	int side = this->Reference()->NSides() - 1;
	this->Reference()->CenterPoint(side,qsi);
	this->ComputeNormal(qsi,normal);
}

void TPZInterfaceElement::ComputeNormal(TPZVec<REAL>&qsi, TPZVec<REAL> &normal){
	int dim = Reference()->Dimension();
	TPZFNMatrix<9> jacobian(dim,dim),jacinv(dim,dim),axes(dim,3);
	REAL detjac;
	this->Reference()->Jacobian(qsi,jacobian,axes,detjac,jacinv);
	this->ComputeNormal(axes,normal);
}

void TPZInterfaceElement::ComputeNormal(TPZFMatrix<REAL> &axes, TPZVec<REAL> &normal){
	
	normal.Resize(3);
	
	TPZCompEl * LeftEl = this->LeftElement();
	TPZCompEl * RightEl = this->RightElement();
    
    TPZMaterial *LeftElMaterial = LeftEl->Material();
    TPZMaterial *RightElMaterial = RightEl->Material();
    
    if (!LeftElMaterial || !RightElMaterial) {
        std::cout << "Interface elements created between elements without material\n";
        std::cout << "Material Ids missing ";
        if(!LeftElMaterial) std::cout << LeftEl->Reference()->MaterialId() << " ";
        if(!RightElMaterial) std::cout << RightEl->Reference()->MaterialId() << " ";
        std::cout << std::endl;
        DebugStop();
    }
	
	//  int dim = Reference()->Dimension();
	// TPZGeoEl *ref = Reference();
	//  int face = ref->NSides()-1;
	//face: lado do elemento bidimensional ou aresta
	//do unidimensional ou canto do ponto
	normal.Resize(3,0.);
	normal.Fill(0.);
	int faceleft,faceright;
	
	TPZManVector<REAL, 3> centleft(3),centright(3),point(3,0.),result(3,0.),xint(3),xvolleft(3),xvolright(3),vec(3),rib(3);
	REAL normalize;
	int i;
	
	faceleft = LeftEl->Reference()->NSides()-1;//lado interior do elemento esquerdo
	faceright = RightEl->Reference()->NSides()-1; // lado interior do element direito
	LeftEl->Reference()->CenterPoint(faceleft,centleft);//ponto centro do elemento de volume
	RightEl->Reference()->CenterPoint(faceright,centright);
	LeftEl->Reference()->X(centleft,xvolleft);
	RightEl->Reference()->X(centright,xvolright);
	for(i=0;i<3;i++) vec[i] = xvolright[i]-xvolleft[i];//n� deve ser nulo
	
	int myinterfacedim = Reference()->Dimension();
	int InterfaceDimension =  LeftEl->Material()->Dimension() - 1;
	if (myinterfacedim != InterfaceDimension) {
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << "the dimension of the interface element " << myinterfacedim << " is not compatible with the dimension of the material " << InterfaceDimension <<
		" Expect trouble ";
#ifdef LOG4CXX
		LOGPZ_WARN(logger,sout.str())
#else
		//		std::cout << sout.str() << std::endl;
#endif
		InterfaceDimension = myinterfacedim;
	}
	
	
	REAL vecnorm = sdot(vec, vec);
	if(InterfaceDimension && vecnorm < 1.e-10)
	{
		int index = fabs(axes(0,0)) < fabs(axes(0,1)) ? 0 : 1;
		index = fabs(axes(0,index)) < fabs(axes(0,2)) ? index : 2;
		vec[index] = 1.;
		LOGPZ_WARN(logger,"Left and Right element centers coincide")
	}
	
	
	switch(InterfaceDimension){
		case 0:
			normal[0] = 1.0;// a normal sempre apontar�na dire� positiva do eixo
			normal[1] = 0.;
			normal[2] = 0.;
			break;
		case 1:
			for(i=0;i<3;i++) rib[i] = axes(0,i);//dire� da aresta
			this->VetorialProd(rib,vec,result);
			this->VetorialProd(result,rib,normal);
			//normalizando a normal
			normalize = 0.;
			for(i=0;i<3;i++) normalize += normal[i]*normal[i];
			if(normalize == 0.0)
			{
				PZError << __PRETTY_FUNCTION__ << " null normal vetor\n";
#ifdef LOG4CXX
				{
					std::stringstream sout;
					Print(sout);
					LOGPZ_DEBUG(logger,sout.str())
				}
#endif
				
				DebugStop();
			}
			normalize = sqrt(normalize);
			for(i=0;i<3;i++) normal[i] = normal[i]/normalize;
			break;
		case 2:{
			TPZManVector<REAL,3> axes1(3), axes2(3);
			for(int iax = 0; iax < 3; iax++){
				axes1[iax] = axes(0,iax);
				axes2[iax] = axes(1,iax);
			}
			this->VetorialProd(axes1,axes2,normal);
		}
			break;
		default:
			PZError << "TPZInterfaceElement::NormalToFace in case that not treated\n";
			normal.Resize(0);
			DebugStop();
			return;
	}
	
	//to guarantee the normal points from left to right neighbours:
	REAL dot = 0.;
	for(i=0; i<3; i++) dot += normal[i]*vec[i];
	if(dot < 0.) {
		for(i=0; i<3; i++) normal[i] = -normal[i];
	}
}

void TPZInterfaceElement::VetorialProd(TPZVec<REAL> &ivet,TPZVec<REAL> &jvet,TPZVec<REAL> &kvet){
	
	kvet.Resize(3);
	kvet[0] =  ivet[1]*jvet[2] - ivet[2]*jvet[1];
	kvet[1] = -ivet[0]*jvet[2] + ivet[2]*jvet[0];
	kvet[2] =  ivet[0]*jvet[1] - ivet[1]*jvet[0];
}

void TPZInterfaceElement::CenterNormal(TPZVec<REAL> &CenterNormal) const{
	const int n = fCenterNormal.NElements();
	CenterNormal.Resize(n);
	for(int i = 0; i < n; i++) CenterNormal[i] = fCenterNormal[i];
}

void TPZInterfaceElement::Normal(TPZFMatrix<REAL> &axes, TPZVec<REAL> &normal){
	TPZGeoEl * gel = this->Reference();
	if(gel->IsLinearMapping()) return this->CenterNormal(normal);
	return this->ComputeNormal(axes, normal);
}

void TPZInterfaceElement::Normal(TPZVec<REAL>&qsi, TPZVec<REAL> &normal){
	TPZGeoEl * gel = this->Reference();
	if(gel->IsLinearMapping()) return this->CenterNormal(normal);
	return this->ComputeNormal(qsi, normal);
}

/**
 * @brief Performs an error estimate on the elemen
 * @param fp function pointer which computes the exact solution
 * @param errors [out] the L2 norm of the error of the solution
 * @param flux [in] value of the interpolated flux values
 */
void TPZInterfaceElement::EvaluateError(std::function<void(const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv)> func,
                           TPZVec<REAL> &errors, bool store_error)
{
	errors.Fill(0.0);
}


/**
 * returns the unique identifier for reading/writing objects to streams
 */
int TPZInterfaceElement::ClassId() const{
    return Hash("TPZInterfaceElement") ^ TPZCompEl::ClassId() << 1;
}

#ifndef BORLAND
template class
TPZRestoreClass<TPZInterfaceElement>;
#endif

/**
 Save the element data to a stream
 */
void TPZInterfaceElement::Write(TPZStream &buf, int withclassid) const
{
	TPZCompEl::Write(buf,withclassid);
	int64_t leftelindex = fLeftElSide.Element()->Index();
	int64_t rightelindex = fRightElSide.Element()->Index();
	if ( (this->Index() < leftelindex) || (this->Index() < rightelindex) ){
		PZError << __PRETTY_FUNCTION__ << endl
		<< "Indices of neighbours are less than interface index:" << endl
		<< "Left: " << leftelindex << ", Right: " << rightelindex << ", this: " << this->Index() << endl;
		DebugStop(); 
	}
	
	int leftside = fLeftElSide.Side();
	int rightside = fRightElSide.Side();
	
	buf.Write(&leftelindex,1);
	buf.Write(&leftside,1);
	buf.Write(&rightelindex,1);
	buf.Write(&rightside,1);
	buf.Write(fCenterNormal);
}

/**
 Read the element data from a stream
 */
void TPZInterfaceElement::Read(TPZStream &buf, void *context)
{
	TPZCompEl::Read(buf,context);
	if (this->Reference()){
		this->Reference()->IncrementNumInterfaces();
	}
	else{
		PZError << "ERROR at " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " - this->Reference() is NULL\n";
		DebugStop();
	}
	int64_t leftelindex;
	int64_t rightelindex;
	int leftside, rightside;
	//  int matid;
	buf.Read(&leftelindex,1);
	buf.Read(&leftside,1);
	buf.Read(&rightelindex,1);
	buf.Read(&rightside,1);
	this->fLeftElSide.SetElement ( Mesh()->ElementVec()[leftelindex]  );
	this->fRightElSide.SetElement( Mesh()->ElementVec()[rightelindex] );
	this->fLeftElSide.SetSide( leftside );
	this->fRightElSide.SetSide( rightside );
	
	buf.Read(fCenterNormal);
}

void TPZInterfaceElement::InitializeElementMatrix(TPZElementMatrix &ef){
	
	TPZInterpolationSpace * left = dynamic_cast<TPZInterpolationSpace*>(this->LeftElement());
	TPZInterpolationSpace * right = dynamic_cast<TPZInterpolationSpace*>(this->RightElement());
	
#ifdef PZDEBUG
	if (!left || !right){
		PZError << "\nError at TPZInterfaceElement::InitializeElementMatrix null neighbour\n";
		ef.Reset();
		DebugStop();
		return;
	}
	if(!left->Material() || !right->Material()){
		PZError << "\n Error at TPZInterfaceElement::InitializeElementMatrix null material\n";
		ef.Reset();
		DebugStop();
		return;
	}
#endif
	
    TPZMaterial *mat = Material();
	const int numdof = mat->NStateVariables();
	ef.fNumStateVars = numdof;
	
	int nshapel = left ->NShapeF();
	int nshaper = right->NShapeF();
	const unsigned int nstatel = left->Material()->NStateVariables();
	const unsigned int nstater = right->Material()->NStateVariables();
	
	TPZManVector<TPZConnect*> ConnectL, ConnectR;
	TPZManVector<int64_t> ConnectIndexL, ConnectIndexR;
	
	this->GetConnects( this->LeftElementSide(),  ConnectL, ConnectIndexL );
	this->GetConnects( this->RightElementSide(), ConnectR, ConnectIndexR );
	const int64_t ncon = ConnectL.NElements() + ConnectR.NElements();
	const int neql = nshapel * nstatel;
	const int neqr = nshaper * nstater;
	const int neq = neql + neqr;
    const int numloadcases = mat->NumLoadCases();
	ef.fMat.Redim(neq,numloadcases);
	ef.fBlock.SetNBlocks(ncon);
	ef.fConnect.Resize(ncon);
	
	int64_t ic = 0;
	int64_t n = ConnectL.NElements();
	for(int64_t i = 0; i < n; i++) {
        TPZConnect &c = left->Connect(i);
		const unsigned int nshape = left->NConnectShapeF(i,c.Order());
		const int con_neq = nstatel * nshape;
#ifdef PZDEBUG
        if(c.NShape() != nshape || c.NState() != nstatel)
        {
            DebugStop();
        }
#endif
		ef.fBlock.Set(ic,con_neq);
		(ef.fConnect)[ic] = ConnectIndexL[i];
		ic++;
	}
	n = ConnectR.NElements();
	for(int64_t i = 0; i < n; i++) {
        TPZConnect &c = right->Connect(i);
		const unsigned int nshape = right->NConnectShapeF(i,c.Order());
		const int con_neq = nstater * nshape;
#ifdef PZDEBUG
        if (c.NShape() != nshape || c.NState() != nstater) {
            DebugStop();
        }
#endif
        
		ef.fBlock.Set(ic,con_neq);
		(ef.fConnect)[ic] = ConnectIndexR[i];
		ic++;
	}
	ef.fBlock.Resequence();
	
}

void TPZInterfaceElement::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef){
	
	TPZInterpolationSpace * left = dynamic_cast<TPZInterpolationSpace*>(this->LeftElement());
	TPZInterpolationSpace * right = dynamic_cast<TPZInterpolationSpace*>(this->RightElement());
	
#ifdef PZDEBUG
	if (!left || !right){
		PZError << "\nError at TPZInterfaceElement::InitializeElementMatrix null neighbour\n";
		ek.Reset();
		ef.Reset();
		DebugStop();
		return;
	}
	if(!left->Material() || !right->Material()){
		PZError << "\n Error at TPZInterfaceElement::InitializeElementMatrix null material\n";
		ek.Reset();
		ef.Reset();
		DebugStop();
		return;
	}
#endif
	
    ek.fMesh = Mesh();
    ek.fType = TPZElementMatrix::EK;
    ef.fType = TPZElementMatrix::EF;
    ef.fMesh = ek.fMesh;
    TPZMaterial *mat = Material();
	const int numdof = mat->NStateVariables();
	ek.fNumStateVars = numdof;
	ef.fNumStateVars = numdof;
	
	int nshapel = left ->NShapeF();
	int nshaper = right->NShapeF();
	const unsigned int nstatel = left->Material()->NStateVariables();
	const unsigned int nstater = right->Material()->NStateVariables();
	
	TPZManVector<TPZConnect*> ConnectL, ConnectR;
	TPZManVector<int64_t> ConnectIndexL, ConnectIndexR;
	
	this->GetConnects( this->LeftElementSide(),  ConnectL, ConnectIndexL );
	this->GetConnects( this->RightElementSide(), ConnectR, ConnectIndexR );
	const int64_t ncon = ConnectL.NElements() + ConnectR.NElements();
	const int neql = nshapel * nstatel;
	const int neqr = nshaper * nstater;
	const int neq = neql + neqr;
    const int numloadcases = mat->NumLoadCases();
	ek.fMat.Redim(neq,neq);
	ef.fMat.Redim(neq,numloadcases);
	ek.fBlock.SetNBlocks(ncon);
	ef.fBlock.SetNBlocks(ncon);
	ek.fConnect.Resize(ncon);
	ef.fConnect.Resize(ncon);
	
#ifdef PZDEBUG
	TPZStack<STATE> solutionvec;
#endif
	
	int64_t ic = 0;
	int64_t n = ConnectL.NElements();
	for(int i = 0; i < n; i++) {
        TPZConnect &c = left->Connect(i);
		const unsigned int nshape = left->NConnectShapeF(i,c.Order());
		const int con_neq = nstatel * nshape;
#ifdef PZDEBUG
        if(c.NShape() != nshape || c.NState() != nstatel)
        {
            DebugStop();
        }
#endif
		ek.fBlock.Set(ic,con_neq );
		ef.fBlock.Set(ic,con_neq);
		(ef.fConnect)[ic] = ConnectIndexL[i];
		(ek.fConnect)[ic] = ConnectIndexL[i];
#ifdef PZDEBUG
		TPZConnect &con = Mesh()->ConnectVec()[ConnectIndexL[i]];
		int64_t seqnum = con.SequenceNumber();
		int blsize = Mesh()->Block().Size(seqnum);
		int64_t pos = Mesh()->Block().Position(seqnum);
		for (int ip=0; ip<blsize; ip++) 
		{
			solutionvec.Push(Mesh()->Solution()(pos+ip)*(STATE)1.e15);
		}
#endif
		ic++;
	}
	n = ConnectR.NElements();
	for(int64_t i = 0; i < n; i++) {
        TPZConnect &c = right->Connect(i);
		const unsigned int nshape = right->NConnectShapeF(i,c.Order());
		const int con_neq = nstater * nshape;
#ifdef PZDEBUG
        if (c.NShape() != nshape || c.NState() != nstater) {
            DebugStop();
        }
#endif
		ek.fBlock.Set(ic,con_neq );
		ef.fBlock.Set(ic,con_neq);
		(ef.fConnect)[ic] = ConnectIndexR[i];
		(ek.fConnect)[ic] = ConnectIndexR[i];
#ifdef PZDEBUG
		TPZConnect &con = Mesh()->ConnectVec()[ConnectIndexR[i]];
		int64_t seqnum = con.SequenceNumber();
		int blsize = Mesh()->Block().Size(seqnum);
		int64_t pos = Mesh()->Block().Position(seqnum);
		for (int ip=0; ip<blsize; ip++) 
		{
			solutionvec.Push(Mesh()->Solution()(pos+ip)*(STATE)1.e15);
		}
#endif
		ic++;
	}
#ifdef PZDEBUG
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
		std::stringstream sout;
		sout.precision(15);
		sout << "ElementSolution = { " << solutionvec << " }/1000000000000000.;" << std::endl;
		LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
#endif
	ek.fBlock.Resequence();
	ef.fBlock.Resequence();
}

void TPZInterfaceElement::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "elemento de interface " << Index() << " Indice deste Material--> " <<this->Material()->Id()<< std::endl;
		
		LOGPZ_DEBUG(logger, sout.str().c_str());
	}
#endif
	
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(Material());
	if(!mat || mat->Name() == "no_name"){
		PZError << "TPZInterfaceElement::CalcStiff interface material null, do nothing\n";
		ek.Reset();
		ef.Reset();
		return;
	}
	
	TPZInterpolationSpace * left = dynamic_cast<TPZInterpolationSpace*>(this->LeftElement());
	TPZInterpolationSpace * right = dynamic_cast<TPZInterpolationSpace*>(this->RightElement());
	
	if (!left || !right){
		PZError << "\nError at TPZInterfaceElement::CalcStiff null neighbour\n";
		ek.Reset();
		ef.Reset();
		return;
	}
	if(!left->Material() || !right->Material()){
		PZError << "\n Error at TPZInterfaceElement::CalcStiff null material\n";
		ek.Reset();
		ef.Reset();
		return;
	}
	
	//TPZMaterialData data;
	const int dim = this->Dimension();
	const int diml = left->Dimension();
	const int dimr = right->Dimension();
//	int nshapel = left ->NShapeF();
//	int nshaper = right->NShapeF();
//	const int nstatel = left->Material()->NStateVariables();
//	const int nstater = right->Material()->NStateVariables();
	
	TPZMaterialData dataright;
	TPZMaterialData dataleft;
    TPZMaterialData data;
    
    mat->FillDataRequirementsInterface(data);
	
	InitMaterialData(dataleft,left);
	InitMaterialData(dataright,right);
    InitMaterialData(data);
	
	
#ifdef PZDEBUG
	if( !dataleft.x.size() ||!dataright.x.size() ){
		PZError << "\n Error at TPZInterfaceElement::CalcStiff null interface\n";
		ek.Reset();
		ef.Reset();
		DebugStop();
	}
#endif
	
    InitializeElementMatrix(ek, ef);
    /*
	TPZManVector<TPZConnect*> ConnectL, ConnectR;
	TPZManVector<int64_t> ConnectIndexL, ConnectIndexR;
	
	this->GetConnects( this->LeftElementSide(),  ConnectL, ConnectIndexL );
	this->GetConnects( this->RightElementSide(), ConnectR, ConnectIndexR );
	const int64_t ncon = ConnectL.NElements() + ConnectR.NElements();
	const int neql = nshapel * nstatel;
	const int neqr = nshaper * nstater;
	const int neq = neql + neqr;
    const int numloadcases = mat->NumLoadCases();
	ek.fMat.Redim(neq,neq);
	ef.fMat.Redim(neq,numloadcases);
	ek.fBlock.SetNBlocks(ncon);
	ef.fBlock.SetNBlocks(ncon);
	ek.fConnect.Resize(ncon);
	ef.fConnect.Resize(ncon);
	
	int64_t ic = 0;
	int64_t n = ConnectL.NElements();
	for(int64_t i = 0; i < n; i++) {
		const int nshape = left->NConnectShapeF(i);
		const int con_neq = nstatel * nshape;
		ek.fBlock.Set(ic,con_neq );
		ef.fBlock.Set(ic,con_neq);
		(ef.fConnect)[ic] = ConnectIndexL[i];
		(ek.fConnect)[ic] = ConnectIndexL[i];
		ic++;
	}
	n = ConnectR.NElements();
	for(int64_t i = 0; i < n; i++) {
		const int nshape = right->NConnectShapeF(i);
		const int con_neq = nstater * nshape;
		ek.fBlock.Set(ic,con_neq );
		ef.fBlock.Set(ic,con_neq);
		(ef.fConnect)[ic] = ConnectIndexR[i];
		(ek.fConnect)[ic] = ConnectIndexR[i];
		ic++;
	}
	ek.fBlock.Resequence();
	ef.fBlock.Resequence();
	
    */
    
    //LOOKING FOR MAX INTERPOLATION ORDER
	int leftmaxp = left->MaxOrder();
	int rightmaxp = right->MaxOrder();
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) 
	{
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "ordem maxima na esquerda " << leftmaxp<<std::endl;
            sout << "ordem maxima na direita " << rightmaxp<<std::endl;
            LOGPZ_DEBUG(logger, sout.str().c_str());
        }
	}
#endif
	
	
	//Max interpolation order
	const int p = (leftmaxp > rightmaxp) ? leftmaxp : rightmaxp;
	
	TPZGeoEl *ref = Reference();
	TPZIntPoints *intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, 2*(p) );

    
	const int npoints = intrule->NPoints();
	
	//integration points in left and right elements: making transformations to neighbour elements
	TPZTransform<> TransfLeft, TransfRight;
	this->ComputeSideTransform(this->LeftElementSide(), TransfLeft);
	this->ComputeSideTransform(this->RightElementSide(), TransfRight);
	
	TPZManVector<REAL,3> intpoint(dim), LeftIntPoint(diml), RightIntPoint(dimr);
	REAL weight;
	//LOOP OVER INTEGRATION POINTS
	for(int ip = 0; ip < npoints; ip++){
		
		intrule->Point(ip,intpoint,weight);
		ref->Jacobian( intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
		weight *= fabs(data.detjac);
		
		TransfLeft.Apply( intpoint, LeftIntPoint );
		TransfRight.Apply( intpoint, RightIntPoint );
		
#ifdef PZDEBUG
		this->CheckConsistencyOfMappedQsi(this->LeftElementSide(), intpoint, LeftIntPoint);
		this->CheckConsistencyOfMappedQsi(this->RightElementSide(), intpoint, RightIntPoint);
#endif


		this->ComputeRequiredData(dataleft, left, LeftIntPoint);
		this->ComputeRequiredData(dataright, right, RightIntPoint);
        data.intLocPtIndex = ip;
		this->ComputeRequiredData(data,intpoint);
        
        // dataleft.x nao esta preenchido!!
        data.x = dataleft.x;
        data.p = dataright.p;
        if (dataleft.p > data.p) {
            data.p = dataleft.p;
        }
        
		mat->ContributeInterface(data,dataleft,dataright, weight, ek.fMat, ef.fMat);
		
	}//loop over integration points
	
	delete intrule;
}


void TPZInterfaceElement::InitializeIntegrationRule(){
    
    TPZInterpolationSpace * left = dynamic_cast<TPZInterpolationSpace*>(this->LeftElement());
    TPZInterpolationSpace * right = dynamic_cast<TPZInterpolationSpace*>(this->RightElement());
    
    //LOOKING FOR MAX INTERPOLATION ORDER
    int leftmaxp = left->MaxOrder();
    int rightmaxp = right->MaxOrder();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "ordem maxima na esquerda " << leftmaxp<<std::endl;
            sout << "ordem maxima na direita " << rightmaxp<<std::endl;
            LOGPZ_DEBUG(logger, sout.str().c_str());
            }
            }
#endif
            
            
            //Max interpolation order
            const int p = (leftmaxp > rightmaxp) ? leftmaxp : rightmaxp;
            
            TPZGeoEl *ref = Reference();
            TPZIntPoints *intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, 2*(p) );
    
    this->SetIntegrationRule(intrule);
}

void TPZInterfaceElement::GetConnects(TPZCompElSide &elside, TPZVec<TPZConnect*> &connects, TPZVec<int64_t> &connectindex){
	
	TPZCompEl * el = elside.Element();
    TPZCompMesh *comm = Mesh()->CommonMesh(el->Mesh());
    TPZCompMesh *elmesh = el->Mesh();
    TPZCompMesh *thismesh = Mesh();
	
	if(el){
		int ncon = el->NConnects();
		connects.Resize(ncon);
		connects.Fill(NULL);
		connectindex.Resize(ncon);
		connectindex.Fill(-1);
		int64_t index;
		for(int i = 0; i < ncon; i++){
			index = el->ConnectIndex(i);
            if (elmesh != thismesh) {
                int64_t superind = elmesh->PutinSuperMesh(index, comm);
                index = thismesh->GetFromSuperMesh(superind, comm);
            }
			connectindex[i] = index;
			connects[i] = &(fMesh->ConnectVec()[ index ]);
		}//for
		
	}
	else{   //if (!el)
		connects.Resize(0);
		connectindex.Resize(0);
	}
	
}//end of method

void TPZInterfaceElement::EvaluateInterfaceJump(TPZSolVec &jump, int opt){
	
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(Material());
	if(!mat || mat->Name() == "no_name"){
		PZError << "TPZInterfaceElement::EvaluateInterfaceJump interface material null, do nothing\n";
		DebugStop();
		return;
	}
	
	TPZMaterialData datal, datar, data;
	const int nstatel = this->LeftElement()->Material()->NStateVariables();
	const int nstater = this->RightElement()->Material()->NStateVariables();
	const int njump = (nstatel > nstater) ? nstatel : nstater;
	
	//LOOKING FOR MAX INTERPOLATION ORDER
	TPZGeoEl *ref = Reference();
	TPZAutoPointer<TPZIntPoints> intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, 2 );
	TPZManVector<int> order(3);
	intrule->GetOrder(order);
	int maxorder = intrule->GetMaxOrder();
	order.Fill(maxorder);
	intrule->SetOrder(order);
	const int npoints = intrule->NPoints();
	TPZManVector<REAL> intpoint(3);
	data.x.Resize(3);
	REAL weight;
	
    int64_t numbersol = jump.size();
    for (int64_t is=0; is<numbersol; is++) {
        jump[is].Resize(njump);
        jump[is].Fill(0.);
    }
    //LOOP OVER INTEGRATION POINTS
	for(int ip = 0; ip < npoints; ip++){
		intrule->Point(ip,intpoint,weight);
		ref->Jacobian( intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
		ref->X(intpoint,data.x);
		weight *= fabs(data.detjac);
		
		//method NeighbourSolution will compute the transformation in this->MapQsi every time it is called
		//(which means for all integration point). Instead of calling NeighbourSolution the whole method
		//may be written here but keeping the transformation computed at first integration point (as done in CalcStiff).
		this->NeighbourSolution(this->LeftElementSide(), intpoint, datal.sol, datal.dsol, datal.axes);
		this->NeighbourSolution(this->RightElementSide(), intpoint, datar.sol, datar.dsol, datar.axes);
		
		if (LeftElement()->NConnects() == 0){
			datal.sol.Resize(0);
			datal.dsol.Resize(0);
		}
		
		if (RightElement()->NConnects() == 0){
			datar.sol.Resize(0);
			datar.dsol.Resize(0);
		}
		TPZSolVec localjump(numbersol);
        for (int64_t is=0; is<numbersol; is++) {
            localjump[is].Resize(njump,0.);
        }
		mat->InterfaceJump(data.x, datal.sol, datar.sol, localjump);
		
        for (int64_t is=0; is<numbersol; is++) {
            if(opt == 0){
                for(int64_t ier = 0; ier < njump; ier++){
                    jump[is][ier] += localjump[is][ier]*localjump[is][ier]*(STATE)weight;
                }
            }
            if(opt == 1){
                for(int64_t ier = 0; ier < njump; ier++){
                    if( fabs(localjump[is][ier]) > fabs(jump[is][ier]) ){
                        jump[is][ier] = fabs( localjump[is][ier] );
                    }//if
                }//for
            }//if
        }		
	}//loop over integration points
	
	//Norma sobre o elemento
    for (int64_t is=0; is<numbersol; is++) {
        if(opt == 0){
            for(int64_t ier = 0; ier < njump; ier++){
                jump[is][ier] = sqrt(jump[is][ier]);
            }//for
        }//if
    }
	
}//method

void TPZInterfaceElement::ComputeErrorFace(int errorid,
										   TPZVec<STATE> &errorL,
										   TPZVec<STATE> &errorR){
	
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(Material());
	if(!mat){
		PZError << "TPZInterfaceElement::ComputeErrorFace interface material null, do nothing\n";
		DebugStop();
		return;
	}
	
	
	TPZInterpolationSpace * left = dynamic_cast<TPZInterpolationSpace*>(this->LeftElement());
	TPZInterpolationSpace * right = dynamic_cast<TPZInterpolationSpace*>(this->RightElement());
	
	if (!left || !right){
		PZError << "\nError at TPZInterfaceElement::ComputeErrorFace null neighbour\n";
		DebugStop();
		return;
	}
	if(!left->Material() || !right->Material()){
		PZError << "\n Error at TPZInterfaceElement::ComputeErrorFace null material\n";
		DebugStop();
		return;
	}
	
	TPZMaterialData data, datal, datar;
	data.SetAllRequirements(true);
	datal.SetAllRequirements(true);
	datar.SetAllRequirements(true);
	this->InitMaterialData(datal,left);
	this->InitMaterialData(datar,right);
	this->InitMaterialData(data);
   
	const int dim = this->Dimension();
	const int diml = left->Dimension();
	const int dimr = right->Dimension();
	
	//LOOKING FOR MAX INTERPOLATION ORDER
	datal.p = left->MaxOrder();
	datar.p = right->MaxOrder();
	//Max interpolation order
	const int p = (datal.p > datar.p) ? datal.p : datar.p;
	
	TPZGeoEl *ref = Reference();
    int intorder = mat->IntegrationRuleOrder(p);
	TPZAutoPointer<TPZIntPoints> intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, intorder );
    if(mat->HasForcingFunction())
    {
        intorder = intrule->GetMaxOrder();
        TPZManVector<int,3> order(ref->Dimension(),intorder);
        intrule->SetOrder(order);
    }
	//   mat->SetIntegrationRule(intrule, p, dim);
	const int npoints = intrule->NPoints();
	
	//integration points in left and right elements: making transformations to neighbour elements
	TPZTransform<> TransfLeft, TransfRight;
	this->ComputeSideTransform(this->LeftElementSide(), TransfLeft);
	this->ComputeSideTransform(this->RightElementSide(), TransfRight);
	
	TPZManVector<REAL,3> intpoint(dim), LeftIntPoint(diml), RightIntPoint(dimr);
	REAL weight;
	//LOOP OVER INTEGRATION POINTS
	for(int ip = 0; ip < npoints; ip++){
		
		intrule->Point(ip,intpoint,weight);
		ref->Jacobian( intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
		weight *= fabs(data.detjac);
		
		this->Normal(data.axes,data.normal);
		
		TransfLeft.Apply( intpoint, LeftIntPoint );
		TransfRight.Apply( intpoint, RightIntPoint );
		
#ifdef PZDEBUG
		this->CheckConsistencyOfMappedQsi(this->LeftElementSide(), intpoint, LeftIntPoint);
		this->CheckConsistencyOfMappedQsi(this->RightElementSide(), intpoint, RightIntPoint);
#endif
		
		//left->ComputeShape(LeftIntPoint, datal.x, datal.jacobian, datal.axes, datal.detjac, datal.jacinv, datal.phi, datal.dphix);
		//right->ComputeShape(RightIntPoint, datar.x, datar.jacobian, datar.axes, datar.detjac, datar.jacinv, datar.phi, datar.dphix);
		
		this->ComputeRequiredData(datal, left, LeftIntPoint);
		this->ComputeRequiredData(datar, right, RightIntPoint);
        //data.SetAllRequirements(true);
        data.fNeedsHSize=true;
        data.intLocPtIndex=ip;
		this->ComputeRequiredData(data,intpoint);
		
		mat->ContributeInterfaceErrors(data, datal, datar, weight,errorL,errorR,errorid);
		
	}//loop over integration points
	
}

void TPZInterfaceElement::Integrate(int variable, TPZVec<STATE> & value){
	const int varsize = this->Material()->NSolutionVariables(variable);
	value.Resize(varsize);
	value.Fill(0.);
}

void TPZInterfaceElement::IntegrateInterface(int variable, TPZVec<REAL> & value){
	TPZMaterial * material = Material();
	if(!material){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " : no material for this element\n";
		DebugStop();
		return;
	}
	if (!this->Reference()){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " : no reference element\n";
		DebugStop();
		return;
	}
	const int dim = this->Dimension();
	TPZInterpolationSpace *left = dynamic_cast<TPZInterpolationSpace*>(this->LeftElement());
	TPZInterpolationSpace *right = dynamic_cast<TPZInterpolationSpace*>(this->RightElement());
	if (!left || !right){
		PZError << "\nError at TPZInterfaceElement::IntegrateInterface null neighbour\n";
		DebugStop();
		return;
	}
	if(!left->Material() || !right->Material()){
		PZError << "\n Error at TPZInterfaceElement::IntegrateInterface null material\n";
		DebugStop();
		return;
	}
    TPZDiscontinuousGalerkin *discgal = dynamic_cast<TPZDiscontinuousGalerkin *>(material);
	
	//local variables
	REAL weight;
	TPZMaterialData data, datal, datar;
	this->InitMaterialData(datal,left);
	this->InitMaterialData(datar,right);
	this->InitMaterialData(data);
	TPZGeoEl *ref = Reference();
	datal.p = left->MaxOrder();
	datar.p = right->MaxOrder();
	TPZManVector<REAL, 3> intpoint(dim,0.);
	const int varsize = material->NSolutionVariables(variable);
	//Max interpolation order
	const int p = (datal.p > datar.p) ? datal.p : datar.p;
	
	//Integration rule
    int intorder = material->IntegrationRuleOrder(p);
	TPZAutoPointer<TPZIntPoints> intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, intorder);
    if(material->HasForcingFunction())
    {
        intorder = intrule->GetMaxOrder();
        TPZManVector<int,3> order(ref->Dimension(),intorder);
        intrule->SetOrder(order);
    }
	//  material->SetIntegrationRule(intrule, p, ref->Dimension());
	
	//loop over integration points
	const int npoints = intrule->NPoints();
	int ip, iv;
	value.Resize(varsize);
	value.Fill(0.);
	TPZManVector<STATE> locval(varsize);
	for(ip=0;ip<npoints;ip++){
		intrule->Point(ip,intpoint,weight);
		ref->Jacobian(intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
		weight *= fabs(data.detjac);
		this->Normal(data.axes,data.normal);
//		this->ComputeSolution(intpoint, data.phi, data.dphix, data.axes, data.sol, data.dsol);
		this->NeighbourSolution(this->LeftElementSide(), intpoint, datal.sol, datal.dsol, datal.axes);
		this->NeighbourSolution(this->RightElementSide(), intpoint, datar.sol, datar.dsol, datar.axes);
		discgal->SolutionDisc(data, datal, datar, variable, locval);
		for(iv = 0; iv < varsize; iv++) {
#ifdef STATE_COMPLEX
			value[iv] += locval[iv].real()*weight;
#else
            value[iv] += locval[iv]*weight;
#endif
		}//for iv
	}//for ip
	
}//method

void TPZInterfaceElement::ComputeSideTransform(TPZCompElSide &Neighbor, TPZTransform<> &transf){
	TPZGeoEl * neighel = Neighbor.Element()->Reference();
	const int dim = this->Dimension();
	TPZTransform<> LocalTransf(dim);
	TPZGeoElSide thisgeoside(this->Reference(), this->Reference()->NSides()-1);
	TPZGeoElSide neighgeoside(neighel, Neighbor.Side());
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "thisgeoside = \n";
        thisgeoside.Print(sout);
        sout << "neighgeoside = \n";
        neighgeoside.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	thisgeoside.SideTransform3(neighgeoside, LocalTransf);
	
	TPZGeoElSide highdim(neighel, neighel->NSides()-1);
	transf = neighgeoside.SideToSideTransform(highdim).Multiply(LocalTransf);
}//ComputeSideTransform

void TPZInterfaceElement::MapQsi(TPZCompElSide &Neighbor, TPZVec<REAL> &qsi, TPZVec<REAL> &NeighIntPoint){
	TPZTransform<> Transf;
	this->ComputeSideTransform(Neighbor, Transf);
	Transf.Apply( qsi, NeighIntPoint );
#ifdef PZDEBUG
	this->CheckConsistencyOfMappedQsi(Neighbor, qsi, NeighIntPoint);
#endif
}//MapQsi

bool TPZInterfaceElement::CheckConsistencyOfMappedQsi(TPZCompElSide &Neighbor, TPZVec<REAL> &qsi, TPZVec<REAL>&NeighIntPoint){
    REAL tol = 0.;
    ZeroTolerance(tol);
    tol *= 100.;
	TPZManVector<REAL,3> FaceXPoint(3), XPoint(3);
	this->Reference()->X( qsi, FaceXPoint);
	Neighbor.Element()->Reference()->X( NeighIntPoint, XPoint);
	int i, n = FaceXPoint.NElements();
	if (n != XPoint.NElements() ){
		PZError << __PRETTY_FUNCTION__ << std::endl
		<< "Face X point and LeftElement X point have not same dimension." << std::endl;
		return false;
	}
	REAL erro = 0.;
	for(i = 0; i < n; i++){
		erro += (XPoint[i] - FaceXPoint[i])*(XPoint[i] - FaceXPoint[i]);
	}
	erro = sqrt(erro);
	if (erro > tol){
		PZError << __PRETTY_FUNCTION__ << std::endl
		<< "Face X point and LeftElement X point are not same." << std::endl;
		this->Print(cout);
		return false;
	}
	return true;
}//void

void TPZInterfaceElement::NeighbourSolution(TPZCompElSide & Neighbor, TPZVec<REAL> & qsi,
                                            TPZSolVec &sol, TPZGradSolVec &dsol,
                                            TPZFMatrix<REAL> &NeighborAxes){
	TPZGeoEl * neighel = Neighbor.Element()->Reference();
	const int neighdim = neighel->Dimension();
	TPZManVector<REAL,3> NeighIntPoint(neighdim);
	this->MapQsi(Neighbor, qsi, NeighIntPoint);
	Neighbor.Element()->ComputeSolution(NeighIntPoint, sol, dsol, NeighborAxes);
}

void TPZInterfaceElement::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
                                          const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol){
	sol.Resize(0);
	dsol.Resize(0);
}

void TPZInterfaceElement::ComputeSolution(TPZVec<REAL> &qsi,
                                          TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes){
	sol.Resize(0);
	dsol.Resize(0);
	axes.Zero();
}

/**
 * @brief Computes solution and its derivatives in the local coordinate qsi.
 * @param qsi master element coordinate
 * @param data contains all elements to compute the solution
 */
void TPZInterfaceElement::ComputeSolution(TPZVec<REAL> &qsi,
                             TPZMaterialData &data)
{
    DebugStop();
}


void TPZInterfaceElement::ComputeSolution(TPZVec<REAL> &qsi,
										  TPZVec<REAL> &normal,
										  TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix<REAL> &leftaxes,
										  TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes){
	this->Normal(qsi, normal);
	this->NeighbourSolution(this->fLeftElSide, qsi, leftsol, dleftsol, leftaxes);
	this->NeighbourSolution(this->fRightElSide, qsi, rightsol, drightsol, rightaxes);
}//method

void TPZInterfaceElement::InitMaterialData(TPZMaterialData &data, TPZInterpolationSpace *elem){
	TPZMaterial *matorig = Material();
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(matorig);
	if (!mat){
		PZError << "FATAL ERROR AT "  << __PRETTY_FUNCTION__ << "\n";
		DebugStop();
	}
    elem->InitMaterialData(data);
	mat->FillDataRequirementsInterface(data);
    if (data.fNeedsNeighborSol) {
        data.fNeedsSol = true;
    }
//
//	//this values needs to be computed only once
	if(data.fNeedsNeighborCenter){
		TPZManVector<REAL,3> qsi(3);
		data.XCenter.Resize(3);
		TPZGeoEl * gel = elem->Reference();
		gel->CenterPoint(gel->NSides()-1,qsi);
		gel->X(qsi,data.XCenter);
	}
	
	data.gelElId = elem->Reference()->Id();
	
}//void


void TPZInterfaceElement::InitMaterialData(TPZMaterialData &data){
//	TPZMaterial *matorig = Material().operator->();
    const int dim = this->Dimension();
    data.axes.Redim(dim,3);
    data.jacobian.Redim(dim,dim);
	data.jacinv.Redim(dim,dim);
	data.x.Resize(3);
	data.normal = this->fCenterNormal;
}

void TPZInterfaceElement::ComputeRequiredData(TPZMaterialData &data,
                                              TPZInterpolationSpace *elem,
                                              TPZVec<REAL> &IntPoint){
    
    data.intGlobPtIndex = -1;
    //elem->ComputeShape(IntPoint, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphix);
    elem->ComputeShape(IntPoint,data);
    
	if (data.fNeedsNeighborSol){
		elem->ComputeSolution(IntPoint, data.phi, data.dphix, data.axes, data.sol, data.dsol );
	}
	
	data.p = elem->MaxOrder();

}
void TPZInterfaceElement::ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &xi)
{

	if (data.fNeedsSol){
        // the interface elements have no approximation space!!
        DebugStop();
	}
	
	if (data.fNeedsHSize){
		const int dim = this->Dimension();
		REAL faceSize;
		if (dim == 0){//it means I am a point
            DebugStop();
            faceSize = 1.;
		}
		else{
			faceSize = 2.*this->Reference()->ElementRadius();//Igor Mozolevski's suggestion. It works well for elements with small aspect ratio
		}
		data.HSize = faceSize;
	}
	
	if(data.fNeedsNormal )//&& !this->Reference()->IsLinearMapping())
    {
		this->ComputeNormal(data.axes,data.normal);
	}
	
}//void

/** @brief adds the connect indexes associated with base shape functions to the set */
void TPZInterfaceElement::BuildCornerConnectList(std::set<int64_t> &connectindexes) const
{
    TPZCompEl *left = LeftElement();
    TPZCompEl *right = RightElement();
#ifdef PZDEBUG
    if (!left || !right ) {
        DebugStop();
    }
#endif
    left->BuildCornerConnectList(connectindexes);
    right->BuildCornerConnectList(connectindexes);
}

