/**
 * @file
 * @brief Contains the implementation of the TPZCompElDisc methods.
 */

#include "pztransfer.h"
#include "pzelmat.h"
#include "pzmatrix.h"
#include "pzelmat.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzerror.h"
#include "pzconnect.h"
#include "TPZMaterial.h"
#include "pzmanvector.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"
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
#include "pzgraphmesh.h"

#include <sstream>
#include <cmath>
#include <algorithm>

#include "time.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include <math.h>
#include <stdio.h>
#include "TPZMaterialDataT.h"
#include "TPZMatSingleSpace.h"
#include "tpzautopointer.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzcompeldisc");
#endif


using namespace pzshape;
using namespace std;

void TPZCompElDisc::SetTensorialShape(){
	this->fShapefunctionType = pzshape::TPZShapeDisc::ETensorial;
	this->SetDegree( this->Degree() );
}

void TPZCompElDisc::SetTotalOrderShape(){
	this->fShapefunctionType = pzshape::TPZShapeDisc::EOrdemTotal;
	this->SetDegree( this->Degree() );
}

void TPZCompElDisc::SetTensorialShapeFull(){
	this->fShapefunctionType = pzshape::TPZShapeDisc::ETensorialFull;
	this->SetDegree( this->Degree() );
}

void TPZCompElDisc::SetTotalOrderShapeFull(){
	this->fShapefunctionType = pzshape::TPZShapeDisc::EOrdemTotalFull;
	this->SetDegree( this->Degree() );
}

void TPZCompElDisc::SetTensorialShape(TPZCompMesh * cmesh){
	if(!cmesh) return;
	int64_t nel = cmesh->NElements();
	for(int64_t iel = 0; iel < nel; iel++){
		TPZCompEl * cel = cmesh->ElementVec()[iel];
		if(!cel) continue;
		TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cel);
		if(!disc) continue;
		disc->SetTensorialShape();
	}
}

void TPZCompElDisc::SetTotalOrderShape(TPZCompMesh * cmesh){
	if(!cmesh) return;
	int64_t nel = cmesh->NElements();
	for(int64_t iel = 0; iel < nel; iel++){
		TPZCompEl * cel = cmesh->ElementVec()[iel];
		if(!cel) continue;
		TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cel);
		if(!disc) continue;
		disc->SetTotalOrderShape();
	}
}

TPZCompElDisc::~TPZCompElDisc() {
	TPZGeoEl * ref = this->Reference();
	if (ref){
		if(ref->Reference() == this) ref->ResetReference();
	}//if (ref)
}

TPZCompElDisc::TPZCompElDisc() : TPZRegisterClassId(&TPZCompElDisc::ClassId),
TPZInterpolationSpace(), fConnectIndex(-1), fExternalShape(), fCenterPoint(3,0.)
{
	this->fShapefunctionType = pzshape::TPZShapeDisc::ETensorial;
	this->fIntRule = NULL;
    fUseQsiEta = true;
    fConstC = 1.;

}

TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh) :
TPZRegisterClassId(&TPZCompElDisc::ClassId),TPZInterpolationSpace(mesh,0), fConnectIndex(-1), fExternalShape(), fCenterPoint(3,0.)
{
	this->fShapefunctionType = pzshape::TPZShapeDisc::ETensorial;  
	this->fIntRule = this->CreateIntegrationRule();
    this->SetTrueUseQsiEta();
}

TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh, const TPZCompElDisc &copy) :
TPZRegisterClassId(&TPZCompElDisc::ClassId),
TPZInterpolationSpace(mesh,copy), fConnectIndex(copy.fIndex), fCenterPoint(copy.fCenterPoint) {
    fShapefunctionType = copy.fShapefunctionType;
    //criando nova malha computacional
    Reference()->SetReference(this);
    //TPZMaterial * mat = copy.Material();
    fConstC = copy.fConstC;
    CreateMidSideConnect();
    this->SetDegree( copy.Degree() );
    //as interfaces foram clonadas
    if (copy.fIntRule){
        this->fIntRule = copy.GetIntegrationRule().Clone();
    }
    else{
        this->fIntRule = NULL;
    }
    this->SetExternalShapeFunction(copy.fExternalShape);
    fUseQsiEta = copy.fUseQsiEta;
}

TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,
                             const TPZCompElDisc &copy,
                             std::map<int64_t,int64_t> &gl2lcConMap,
                             std::map<int64_t,int64_t> &gl2lcElMap) : 
TPZRegisterClassId(&TPZCompElDisc::ClassId),TPZInterpolationSpace(mesh,copy), fConnectIndex(-1), fCenterPoint(copy.fCenterPoint)
{
	fShapefunctionType = copy.fShapefunctionType;
	//TPZMaterial * mat = copy.Material();
	gl2lcElMap[this->fIndex] = copy.fIndex;
	
	if (copy.fIntRule){
		this->fIntRule = copy.GetIntegrationRule().Clone();
	}
	else{
		this->fIntRule = NULL;
	}
	
	this->SetExternalShapeFunction(copy.fExternalShape);
    fUseQsiEta = copy.fUseQsiEta;
}

TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,TPZGeoEl *ref) :
TPZRegisterClassId(&TPZCompElDisc::ClassId),TPZInterpolationSpace(mesh,ref), fConnectIndex(-1), fExternalShape(), fCenterPoint(this->Dimension(),0.)
{
	this->fShapefunctionType = pzshape::TPZShapeDisc::ETensorial;  
	ref->SetReference(this);
	CreateMidSideConnect();
	this->SetDegree( fMesh->GetDefaultOrder() );
	
    //criando os elementos interface
    //Now, the interface is created by calling the CreateIntefaces method in class TPZApproxSpace.
	//CreateInterfaces();
	
	this->fIntRule = this->CreateIntegrationRule();
    SetTrueUseQsiEta();
	
}

REAL TPZCompElDisc::NormalizeConst()
{
	TPZGeoEl *ref = Reference();
	//greater distance between the point and the vertices of the inner element
	int nnodes = ref->NNodes(),i;
	if(nnodes == 1) return 1.0;//point element
	REAL maxdist,dist;
	int64_t inode = ref->NodeIndex(0);//first node of the element
	TPZGeoNode node = ref->Mesh()->NodeVec()[inode];
	maxdist = pow(node.Coord(0)-fCenterPoint[0],(REAL)2.)+pow(node.Coord(1)-fCenterPoint[1],(REAL)2.);   // I do casting because the developer use double exponent (2.) - I think better using int (2) - Jorge
	maxdist += pow(node.Coord(2)-fCenterPoint[2],(REAL)2.);
	maxdist = sqrt(maxdist);
	for(i=1;i<nnodes;i++){
		inode = ref->NodeIndex(i);//subsequent node
		node = ref->Mesh()->NodeVec()[inode];
		dist = pow(node.Coord(0)-fCenterPoint[0],(REAL)2.)+pow(node.Coord(1)-fCenterPoint[1],(REAL)2.);
		dist += pow(node.Coord(2)-fCenterPoint[2],(REAL)2.);
		dist = sqrt(dist);
		if(maxdist < dist) maxdist = dist;
	}
	return maxdist;
}

/**
 * @brief Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */
void TPZCompElDisc::InitMaterialData(TPZMaterialData &data)
{
    TPZInterpolationSpace::InitMaterialData(data);
    data.XCenter = fCenterPoint;
}

void TPZCompElDisc::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                                 TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                                 REAL &detjac, TPZFMatrix<REAL> &jacinv,
                                 TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphix){
	TPZGeoEl * ref = this->Reference();
	if (!ref){
		PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
		return;
	}//if
    if(fUseQsiEta==true){
        ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
        this->Shape(intpoint,phi,dphi);
        this->Convert2Axes(dphi,jacinv,dphix);
        ref->X(intpoint, X);
    }else{
        dphi.Zero();
        if(Degree() == 0 && !fExternalShape)
        {
            phi(0,0) = 1.;
            return;
        }
        ref->X(intpoint, X);
        this->ShapeX(X,phi,dphix);
        //axes is identity in discontinuous elements
        axes.Resize(dphix.Rows(), 3);
        axes.Zero();
        for(int i = 0; i < axes.Rows(); i++) axes(i,i) = 1.;
    }
	
}

void TPZCompElDisc::ComputeShape(TPZVec<REAL> &intpoint,TPZMaterialData &data){
    this->ComputeShape(intpoint, data.x, data.jacobian, data.axes,data.detjac, data.jacinv, data.phi, data.fDPhi, data.dphix);
}


void TPZCompElDisc::Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
    if(fUseQsiEta==true){
        this->ShapeX(qsi,phi,dphi);
    }else{

        TPZManVector<REAL,4> x(3);
        this->Reference()->X(qsi,x);
        this->ShapeX(x,phi,dphi);
    }
}


void TPZCompElDisc::ShapeX(TPZVec<REAL> &X, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi){
	const int Degree = this->Degree();
	if(Degree < 0) return;
	const int dim = this->Dimension();
	TPZShapeDisc::Shape(dim, fConstC,fCenterPoint,X,Degree,phi,dphi,fShapefunctionType);
	
	//now appending external shape functions
	this->AppendExternalShapeFunctions(X,phi,dphi);
}//method

void TPZCompElDisc::AppendExternalShapeFunctions(TPZVec<REAL> &X, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi){
	
	//adding external shape functions whether they exist
	if(!this->fExternalShape.operator ->()) return;
	
	TPZManVector<STATE> extPhi;
	TPZFNMatrix<100,STATE> extDPhi;
    TPZFNMatrix<100> ThisPhi, ThisDPhi;
	
	//computing external shape functions
	this->fExternalShape->Execute(X, extPhi, extDPhi);
	
	//now appending all shape functions
	{
		
		const int ndiscphi = TPZShapeDisc::NShapeF(this->Degree(),this->Dimension(),fShapefunctionType);
		const int nextphi = this->fExternalShape->NFunctions();
		
#ifdef PZDEBUG
		if(phi.Cols() != 1){
			PZError << "\nError at " << __PRETTY_FUNCTION__ << "\n";
			DebugStop();
		}
#endif  
		
		ThisPhi = phi;
		phi.Resize(ndiscphi+nextphi,1);
		phi.Zero();
		for(int i = 0; i < ndiscphi; i++) {
			phi(i,0) = ThisPhi(i,0);
		}
		for(int i = 0; i < nextphi; i++) {
#ifdef STATE_COMPLEX
			phi(i+ndiscphi,0) = extPhi[i].real();
#else
			phi(i+ndiscphi,0) = (REAL)extPhi[i];
#endif
		}
		
	} 
	
	{
#ifdef PZDEBUG
		if(dphi.Rows() > extDPhi.Rows()){
			PZError << "\nError at " << __PRETTY_FUNCTION__ << "\n";
			DebugStop();
		}
#endif
		
		ThisDPhi = dphi;
		const int ndiscdphi = TPZShapeDisc::NShapeF(this->Degree(),this->Dimension(),fShapefunctionType);
		const int nextdphi = this->fExternalShape->NFunctions();
		const int64_t nderiv = ThisDPhi.Rows();
		dphi.Resize(nderiv, ndiscdphi+nextdphi);
		dphi.Zero();
		for(int64_t i = 0; i < nderiv; i++){
			for(int j = 0; j < ndiscdphi; j++){
				dphi(i,j) = ThisDPhi(i,j);
			}
		}
		for(int64_t i = 0; i < nderiv; i++){
			for(int j = 0; j < nextdphi; j++){
#ifdef STATE_COMPLEX
				dphi(i,j+ndiscdphi) = extDPhi(i,j).real();
#else
                dphi(i,j+ndiscdphi) = (REAL)extDPhi(i,j);
#endif
			}
		}
	}
}

void TPZCompElDisc::Print(std::ostream &out) const{
	
	out << "\nDiscontinous element : \n";
	TPZCompEl::Print(out);
	if(Reference()) out << "\tGeometric reference index : " << Reference()->Index() << endl;
	out << "\tMaterial id : " << Reference()->MaterialId() << endl
	<< "\tDegree of interpolation : " <<  this->Degree() << endl
	<< "\tConnect index : " << fConnectIndex << endl
    << "\tUsing qsieta : " << (int)fUseQsiEta << endl
	<< "\tNormalizing constant : " << fConstC << endl
	<< "\tCenter point of the element : ";
    if(fUseQsiEta == false)
    {
        int size = fCenterPoint.NElements(),i;
        for(i=0;i<size-1;i++) out << fCenterPoint[i] << " , ";
        out << "stored xy : " << fCenterPoint[i] << endl;
    }
    else
    {
        TPZGeoEl *Ref = Reference();
        if (Ref) {
            TPZManVector<REAL,3> xcenter(3),loccenter(fCenterPoint);
            Ref->X(loccenter, xcenter);
            out << "stored qsieta " << fCenterPoint << " computed xy " << xcenter << std::endl;
        }
    }
	out << "\tDimension : " << this->Dimension() << endl;
}

int64_t TPZCompElDisc::ConnectIndex(int side) const{
	return fConnectIndex;
}

int TPZCompElDisc::NConnects() const{
	return (fConnectIndex !=-1);
}

int64_t TPZCompElDisc::CreateMidSideConnect(){
	// primeiro sao criados os elementos de volume depois os elementos BC associados aos seus lados
	// num estagio inicial o elemento BC eh acoplado ao elemento ELV de volume de tal forma
	// que ambos sao vizinhos
	// o elemento BC nao pode ser dividido se o elemento ELV associado nao for dividido primeiro
	// caso o elemento ELV for dividido, entao o elemento BC associado deveria ser dividido
	// tambem para manter a CC consistente com a malha
	// caso ELV for dividido e BC nao eh entao ELV sera LowerLevelElement do elemento BC
	TPZMaterial * material = Material();
	if(!material){
		PZError << "\nTPZCompElDisc::CreateMidSideConnect Material nulo\n";
		return -1;
	}
	
	TPZGeoEl *ref = Reference();
	TPZStack<TPZCompElSide> list;
	int nsides = ref->NSides();
	int dimgrid = Mesh()->Dimension();
	int dim = Dimension();
	int existsconnect = 0;
	
	if(dimgrid == dim){
		//este eh um elemento de volume
		//procura-se elemento superposto
		TPZCompElSide(this,nsides-1).EqualLevelElementList(list,0,0);
		int64_t size = list.NElements(), i;
		for(i=0;i<size;i++){
			int dimel = list[i].Element()->Reference()->Dimension();
			if(dimel == dimgrid){
				int64_t connectindex = list[i].Element()->ConnectIndex(0);
				list[i].Element()->SetConnectIndex(0,connectindex);
				existsconnect = 1;
				break;
			}
		}
	}
	
	if(dim != dimgrid/* - 1*/){ //dimgrid - 1 = interface dimension
		// o atual eh um elemento BC
		fConnectIndex = -1;//=> return NshapeF() = 0
		return fConnectIndex;
	}
	
	if(!existsconnect){
		//o atual eh um elemento de volume e
		//nao achou-se um elemento superposto
		int nvar = Material()->NStateVariables();
		int nshape = 1; /// if the CompelDisc has no connect NShapeF() return -1
        /// therefore nshape is initialize to 1
		int64_t newnodeindex = Mesh()->AllocateNewConnect(nshape,nvar,0);
		SetConnectIndex(0,newnodeindex);
		TPZConnect &newnod = Mesh()->ConnectVec()[newnodeindex];
        newnod.IncrementElConnected();
        nshape = this->NShapeF();
        newnod.SetNShape(nshape);
		int64_t seqnum = newnod.SequenceNumber();
		Mesh()->Block().Set(seqnum,nvar*nshape);
		Mesh()->ConnectVec()[fConnectIndex].IncrementElConnected();
	}
	
	return fConnectIndex;
}

int TPZCompElDisc::NShapeF() const {
	if(fConnectIndex == -1) return 0;
	//deve ter pelo menos um connect
	
	int nExtShape = 0;
	if(fExternalShape) nExtShape = fExternalShape->NFunctions();
	
	int dim = Dimension();
	const int discShape = TPZShapeDisc::NShapeF(this->Degree(),dim,fShapefunctionType);
	return (discShape + nExtShape);  
	
}

int TPZCompElDisc::NConnectShapeF(int inod, int order) const {
#ifdef PZDEBUG2
	if (inod != 0){
		PZError << "\nFATAL ERROR AT " << __PRETTY_FUNCTION__
		<< " - TPZCompElDisc has only one connect and inod = " << inod << "\n";
	}
#endif
    if (order != Degree()) {
        DebugStop();
    }
	return this->NShapeF();
}

void TPZCompElDisc::InternalPoint(TPZVec<REAL> &point){
	//ponto deformado
	point.Resize(3,0.);
	point[0] = fCenterPoint[0];
	point[1] = fCenterPoint[1];
	point[2] = fCenterPoint[2];
}

REAL TPZCompElDisc::SizeOfElement()
{
	TPZGeoEl *ref = Reference();
	
	int dim = ref->Dimension();
	int side = ref->NSides()-1;
	if(dim == 2) ref->SideArea(side);
	if(!dim || dim > 2){
		PZError << "TPZCompElDisc::SizeOfElement case not permited\n";
		return 0.;
	}
	if(dim == 1){
		TPZGeoNode node0 = Mesh()->Reference()->NodeVec()[ref->NodeIndex(0)];
		TPZGeoNode node1 = Mesh()->Reference()->NodeVec()[ref->NodeIndex(1)];
		TPZVec<REAL> no0(3),no1(3);
		for(int i=0;i<3;i++){
			no0[i] = node0.Coord(i);
			no1[i] = node1.Coord(i);
		}
		return ref->Distance(no0,no1);
	}
	PZError << "TPZCompElDisc::SizeOfElement this in case that it is not contemplated\n";
	return 0.;
}

void TPZCompElDisc::Divide(int64_t index,TPZVec<int64_t> &subindex,int interpolatesolution){
	
	Mesh()->SetAllCreateFunctions(*this);
	
	if (Mesh()->ElementVec()[index] != this) {
		PZError << "TPZInterpolatedElement::Divide index error";
		subindex.Resize(0);
		return;
	}
	TPZMaterial * material = Material();
	if(!material)
	{
		PZError << __PRETTY_FUNCTION__ << " no material\n";
		return;
	}
	
	TPZGeoEl *ref = Reference();
	RemoveInterfaces();
	
#ifdef PZDEBUG2
	if(0){//TESTE
		ofstream mesh("MALHADIV.out");//TESTE
		Mesh()->Reference()->Print(mesh);//TESTE
		Mesh()->Print(mesh);//TESTE
		mesh.flush();  //TESTE
		mesh.close();//TESTE
	}//TESTE
#endif
	
	//divide o elemento geom�rico
	int nsubs = ref->NSubElements();
	subindex.Resize(nsubs);
	TPZManVector<TPZGeoEl *> geosubs(nsubs);
	ref->Divide(geosubs);
	if(!geosubs.NElements()) {
		subindex.Resize(0);
		return;
	}
    
    // The refinement pattern may be adjusted during the division process. 02/25/2013
    nsubs = ref->NSubElements();
	subindex.Resize(nsubs);
	
	this->Mesh()->ElementVec()[index] = NULL;
	ref->ResetReference();
	TPZCompElDisc *discel;
	int i,deg;
	deg = this->Degree();
	
	for (i=0;i<nsubs;i++){
		TPZCompEl* cel = Mesh()->CreateCompEl(geosubs[i]);//aqui
        subindex[i] = cel->Index();
		//new TPZCompElDisc(*Mesh(),geosubs[i],subindex[i]);
		discel = dynamic_cast<TPZCompElDisc *> (Mesh()->ElementVec()[subindex[i]]);
		if (!discel){
			std::stringstream mess;
			mess << __PRETTY_FUNCTION__ << " - discel is NULL ";
			LOGPZ_ERROR(logger, mess.str() );
			continue;
		}
		discel->fShapefunctionType = this->fShapefunctionType;
		discel->SetDegree(deg);
	}
	
	if (interpolatesolution){
		Mesh()->ExpandSolution();
		for(i=0; i<nsubs; i++) {
			discel = dynamic_cast<TPZCompElDisc *> ( Mesh()->ElementVec()[subindex[i]] );
			if (!discel){
				std::stringstream mess;
				mess << __PRETTY_FUNCTION__ << " - discel is NULL ";
				LOGPZ_ERROR(logger, mess.str() );
				continue;
			}
			if(discel->Dimension() < material->Dimension()) continue;//elemento BC
			discel->InterpolateSolution(*this);
		}
	}//if interpolate
	
	delete this;
}

template<class TVar>
void TPZCompElDisc::SolutionXInternal(TPZVec<REAL> &x, TPZVec<TVar> &uh){
	TPZCompMesh *finemesh = Mesh();
	TPZBlock &fineblock = finemesh->Block();
	int nstate = Material()->NStateVariables();
	TPZFMatrix<TVar> &FineMeshSol = finemesh->Solution();
	int matsize = NShapeF(),dim = Dimension();
	TPZFMatrix<REAL> phix(matsize,1,0.);
	TPZFMatrix<REAL> dphix(dim,matsize,0.);
	ShapeX(x,phix,dphix);
	TPZConnect *df = &Connect(0);
	int64_t dfseq = df->SequenceNumber();
	int dfvar = fineblock.Size(dfseq);
	int64_t pos   = fineblock.Position(dfseq);
	int iv = 0,d;
	uh.Fill(0.);
	for(d=0; d<dfvar; d++) {
		uh[iv%nstate] += (TVar)phix(iv/nstate,0)*FineMeshSol(pos+d,0);
		iv++;
	}
}

void TPZCompElDisc::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension)
{
	TPZGeoEl *ref = Reference();
	int matid = Material()->Id();
	int nsides = ref->NSides();
    bool to_postpro = grmesh.Material_Is_PostProcessed(matid);
	if(dimension == 2 && to_postpro){
		if(nsides == 9){
			new TPZGraphElQ2dd(this,&grmesh);
			return;
		}
		if(nsides == 7){
			new TPZGraphElT2dMapped(this,&grmesh);
			return;
		}
	}//2d
	
	if(dimension == 3 && to_postpro){
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
	
	if(dimension == 1 && this->NConnects() /* && mat > 0 */){
		if(nsides == 3){
		  new TPZGraphEl1dd(this,&grmesh);
			return;
		}
	}//1d
}

int TPZCompElDisc::NSides(){
	
	return Reference()->NSides();
}

int TPZCompElDisc::NInterfaces(){
	
	int nsides = this->NSides();
	
	switch( nsides )
    {
		case 3: //line
			return 2;
			//      break;
			
		case 7: //triangle
			return 3;
			//      break;
			
		case 9: //square
			return 4;
			
		case 15: // Tetrahedra.
			return 4;
			//      break;
			
		case 19: // Prism.
			return 5;
			//     break;
			
		case 21: // Pyramid.
			return 6;
			//     break;
			
		case 27: // Hexaedra.
			return 8;
			//     break;
			
		default:
			PZError << "TPZCompElDisc::NFaces() - Unknown element shape!" << endl;
			DebugStop();
	}
	return 0;
}

//#include "TPZAgglomerateEl.h"
void TPZCompElDisc::AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight){
	
	int i,npoints;
	TPZVec<REAL> pt(3),x(3,0.0);
	TPZFMatrix<REAL> jacobian(3,3),jacinv(3,3),axes(3,3);
	REAL detjac,wt;
	
	TPZGeoEl *subgel = Reference();
	if(!subgel) PZError << "TPZCompElDisc::AccumulateIntegrationRule data error, null geometric reference\n";
	TPZIntPoints *rule = subgel->CreateSideIntegrationRule(subgel->NSides()-1,degree);
	npoints = rule->NPoints();
	
	for(i=0;i<npoints;i++){
		
		rule->Point(i,pt,wt);
		subgel->Jacobian(pt,jacobian,axes,detjac,jacinv);
		subgel->X(pt, x);
		
		point.Push(x[0]);
		point.Push(x[1]);
		point.Push(x[2]);
		
		weight.Push(wt * fabs(detjac));
	}
	delete rule;
}


void TPZCompElDisc::CenterPoint(TPZVec<REAL> &center){
	
	TPZGeoEl *ref = Reference();
	if(ref || Type() == EDiscontinuous){
		ref->CenterPoint(ref->NSides()-1,center);
		return;
	} else {//aglomerado
		//     TPZStack<TPZCompEl *> elvec;
		//     dynamic_cast<TPZAgglomerateElement *>(this)->ListOfDiscEl(elvec);
		//     TPZGeoEl *ref = elvec[0]->Reference();
		//     ref->CenterPoint(ref->NSides()-1,center);
		PZError << "TPZCompElDisc::CenterPoint center points not exists!\n";
	}
}

void TPZCompElDisc::BuildTransferMatrix(TPZCompElDisc &coarsel, TPZTransfer<STATE> &transfer){
	// accumulates the transfer coefficients between the current element and the
	// coarse element into the transfer matrix, using the transformation t
	
	int locnshape = NShapeF();
	int cornshape = coarsel.NShapeF();
	
	// compare interpolation orders
	// the interpolation order of this >= that interpolation order of coarse
	int locdeg = Degree(), coarsedeg = coarsel.Degree();
	if(coarsedeg > locdeg) {
		SetDegree(coarsedeg);
	}
	
	TPZFNMatrix<500,STATE> loclocmat(locnshape,locnshape);
	TPZFNMatrix<200,STATE> loccormat(locnshape,cornshape);
	loclocmat.Zero();
	loccormat.Zero();
	
	TPZGeoEl *ref = Reference();
	int integdeg = locdeg >= coarsedeg ? locdeg : coarsedeg;
	TPZIntPoints *intrule = ref->CreateSideIntegrationRule(ref->NSides()-1,2*integdeg);
	int dimension = Dimension();
	
	TPZFNMatrix<50> locphi(locnshape,1);
	TPZFNMatrix<150> locdphi(dimension,locnshape);
	locphi.Zero();
	locdphi.Zero();
	// derivative of the shape function
	// in the master domain
	
	TPZFMatrix<REAL> corphi(cornshape,1);
	TPZFMatrix<REAL> cordphi(dimension,cornshape);
	// derivative of the shape function
	// in the master domain
	
	TPZManVector<REAL> int_point(dimension);
	TPZFNMatrix<9> jacobian(dimension,dimension);
	TPZFMatrix<REAL> jacinv(dimension,dimension);
	TPZFNMatrix<9> axes(3,3);
	TPZManVector<REAL> x(3);
	
	int_point.Fill(0.,0);
	REAL jac_det = 1.;
	ref->Jacobian( int_point, jacobian , axes, jac_det, jacinv);
	REAL multiplier = 1./jac_det;
	
	int numintpoints = intrule->NPoints();
	REAL weight;
	int lin,ljn,cjn;
	
	for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {
		
		intrule->Point(int_ind,int_point,weight);
		ref->Jacobian( int_point, jacobian , axes, jac_det, jacinv);
		ref->X(int_point, x);
		Shape(int_point,locphi,locdphi);
		weight *= jac_det;
		corphi.Zero();
		cordphi.Zero();
		coarsel.Shape(int_point,corphi,cordphi);
		
		for(lin=0; lin<locnshape; lin++) {
			for(ljn=0; ljn<locnshape; ljn++) {
				loclocmat(lin,ljn) += weight*locphi(lin,0)*locphi(ljn,0)*multiplier;
			}
			for(cjn=0; cjn<cornshape; cjn++) {
				loccormat(lin,cjn) += weight*locphi(lin,0)*corphi(cjn,0)*multiplier;
			}
		}
		jacobian.Zero();
	}
	loclocmat.SolveDirect(loccormat,ELDLt);
	
	
	int64_t locblockseq = Connect(0).SequenceNumber();
	TPZStack<int> globblockvec;
	int numnonzero = 0;
	int cind = coarsel.ConnectIndex(0);
	TPZConnect &con = coarsel.Mesh()->ConnectVec()[cind];
	int corblockseq = con.SequenceNumber();
	if(locnshape == 0 || cornshape == 0)
		PZError << "TPZCompElDisc::BuilTransferMatrix error I\n";
	TPZFMatrix<STATE> smalll(locnshape,cornshape,0.);
	loccormat.GetSub(0,0,locnshape,cornshape,smalll);
	REAL tol = Norm(smalll);
	if(tol >= 1.e-10) {
		globblockvec.Push(corblockseq);
		numnonzero++;
	}
	if(!numnonzero)
		PZError << "TPZCompElDisc::BuilTransferMatrix error II\n";
	if(transfer.HasRowDefinition(locblockseq))
		PZError << "TPZCompElDisc::BuilTransferMatrix error III\n";
	transfer.AddBlockNumbers(locblockseq,globblockvec);
	if(cornshape == 0 || locnshape == 0)
		PZError << "TPZCompElDisc::BuilTransferMatrix error IV\n";
	loccormat.GetSub(0,0,locnshape,cornshape,smalll);
	transfer.SetBlockMatrix(locblockseq,globblockvec[0],smalll);
	
	SetDegree(locdeg);
	delete intrule;
}

void TPZCompElDisc::AccumulateVertices(TPZStack<TPZGeoNode *> &nodes) {
	TPZGeoEl *geo = Reference();
	
	//Code isnt place to chat
	//#warning "Este metodo nao funciona para aglomerados contendo aglomerados"
	if(!geo) {
		PZError <<  "TPZCompElDisc::AccumulateVertices null reference\n";
		return;
	}
	int nvertices = geo->NNodes();
	int l;
	for(l=0;l<nvertices;l++) nodes.Push( geo->NodePtr(l) );
}

void TPZCompElDisc::SetDegree(int degree) {
	this->fPreferredOrder = degree;
	if (fConnectIndex < 0) return;
	TPZConnect &c = Mesh()->ConnectVec()[fConnectIndex];
	c.SetOrder(degree,fConnectIndex);
	int64_t seqnum = c.SequenceNumber();
	int nvar = 1;
	TPZMaterial * mat = Material();
	if(mat) nvar = mat->NStateVariables();
	int nshapef = this->NShapeF();
    c.SetNShape(nshapef);
    c.SetNState(nvar);
	Mesh()->Block().Set(seqnum,nshapef*nvar);
}

void TPZCompElDisc::SetExternalShapeFunction(TPZAutoPointer<TPZFunction<STATE> > externalShapes){
	this->fExternalShape = externalShapes;
	//in order of ajust block size because NShapeF may have changed
	if (fConnectIndex < 0) return;
	TPZConnect &c = Mesh()->ConnectVec()[fConnectIndex];
	int64_t seqnum = c.SequenceNumber();
	int nvar = 1;
	TPZMaterial * mat = Material();
	if(mat) nvar = mat->NStateVariables();
	int nshapef = this->NShapeF();
    c.SetNShape(nshapef);
    c.SetNState(nvar);
	Mesh()->Block().Set(seqnum,nshapef*nvar);
}

bool TPZCompElDisc::HasExternalShapeFunction(){
	return (this->fExternalShape);
}

/**
 * returns the unique identifier for reading/writing objects to streams
 */
int TPZCompElDisc::ClassId() const{
    return Hash("TPZCompElDisc") ^ TPZInterpolationSpace::ClassId() << 1;
}

#ifndef BORLAND
template class
TPZRestoreClass< TPZCompElDisc>;
#endif

/**
 Save the element data to a stream
 */
void TPZCompElDisc::Write(TPZStream &buf, int withclassid) const
{
	TPZInterpolationSpace::Write(buf,withclassid);
	buf.Write(fCenterPoint);
	buf.Write(&fConnectIndex,1);
	buf.Write(&fConstC,1);
	int matid = Material()->Id();
	buf.Write(&matid,1);
	int shapetype = fShapefunctionType;
	buf.Write(&shapetype,1);
	if(this->fExternalShape.operator ->()){
		int um = 1;
		buf.Write(&um,1);
		this->fExternalShape->Write(buf,withclassid);
	}
	else{
		int zero = 0;
		buf.Write(&zero,1);
	}
	int hasIntRule = this->fIntRule ? 1 : 0;
    buf.Write(&hasIntRule,1);
	if( this->fIntRule){
		TPZManVector<int> pOrder(3);
		this->fIntRule->GetOrder(pOrder);
		buf.Write(pOrder);
		TPZGeoEl *gel = this->Reference();
		int hasAssociatedGel = gel ? 1 : 0;
		buf.Write(&hasAssociatedGel);
		if(gel){
            TPZPersistenceManager::WritePointer(gel, &buf);
		}
	}
}

/**
 Read the element data from a stream
 */
void TPZCompElDisc::Read(TPZStream &buf, void *context)
{
	TPZInterpolationSpace::Read(buf,context);
	buf.Read<3>(fCenterPoint);
	buf.Read(&fConnectIndex,1);
	buf.Read(&fConstC,1);
	int matid;
	buf.Read(&matid,1);
	//  fMaterial = Mesh()->FindMaterial(matid);
	int shapetype;
	buf.Read(&shapetype,1);
	fShapefunctionType = (TPZShapeDisc::MShapeType) shapetype;
	int hasExternalShape;
	buf.Read(&hasExternalShape,1);
	if(hasExternalShape == 1){
// #warning Como faz?
// #warning    this->fExternalShape->
	}	
	int hasIntRule;
	buf.Read(&hasIntRule,1);
	if( hasIntRule ){
		TPZManVector<int> pOrder(3);
		buf.Read(pOrder);
		int hasAssociatedGel;
		buf.Read(&hasAssociatedGel,1);
		if(hasAssociatedGel){
            TPZGeoEl *gel = dynamic_cast<TPZGeoEl *>(TPZPersistenceManager::GetInstance(&buf));
			TPZAutoPointer<TPZIntPoints> result = gel->CreateSideIntegrationRule(gel->NSides()-1, 0);
			result->SetOrder(pOrder);
			this->fIntRule = result;
		}
		else{
			this->fIntRule = nullptr;
		}
	}
	else{
		this->fIntRule = nullptr;
	}
	
}

TPZAutoPointer<TPZIntPoints> TPZCompElDisc::CreateIntegrationRule() const{
	TPZGeoEl * gel = this->Reference();
	if(gel){
		const int integ = Max( 2 * this->Degree()+1, 0);
		TPZAutoPointer<TPZIntPoints> result = gel->CreateSideIntegrationRule(gel->NSides()-1,integ);
		return result;
	}
	else{
		return NULL;
	}
}

const TPZIntPoints &TPZCompElDisc::GetIntegrationRule() const {
	if(this->fIntRule == 0){
		DebugStop();
	}
  if (fIntegrationRule) {
    return *fIntegrationRule;
  }
  else
  {
    return *(fIntRule.operator->());
  }
}

TPZIntPoints &TPZCompElDisc::GetIntegrationRule() {
	if(this->fIntRule == 0){
		DebugStop();
	}
  if (fIntegrationRule)
  {
    return *fIntegrationRule;
  }
  else
  {
    return *(fIntRule.operator->());
  }
}

int TPZCompElDisc::MaxOrderExceptExternalShapes(){
	return TPZInterpolationSpace::MaxOrder();
}

int TPZCompElDisc::MaxOrder(){
	int result = TPZInterpolationSpace::MaxOrder();
	if(this->fExternalShape.operator ->()){
		int extOrder = this->fExternalShape->PolynomialOrder();
		if(extOrder > result) result = extOrder;
	}
	return result;
}

/** @brief adds the connect indexes associated with base shape functions to the set */
void TPZCompElDisc::BuildCornerConnectList(std::set<int64_t> &connectindexes) const
{
    connectindexes.insert(ConnectIndex(0));
}

void TPZCompElDisc::SetCreateFunctions(TPZCompMesh* mesh) {
    mesh->SetAllCreateFunctionsDiscontinuous();
}

int TPZCompElDisc::Degree() const {
    if (fConnectIndex < 0) {
        return -1;
    }
    return this->Connect(0).Order();
}

#define INSTANTIATE_METHODS(TVar) \
template \
 void TPZCompElDisc::SolutionXInternal<TVar>(TPZVec<REAL> &x, TPZVec<TVar> &uh);

INSTANTIATE_METHODS(STATE)
INSTANTIATE_METHODS(CSTATE)
#undef INSTANTIATE_METHODS
