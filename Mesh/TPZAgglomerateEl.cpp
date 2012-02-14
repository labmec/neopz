/**
 * @file
 * @brief Contains the implementation of the TPZAgglomerateElement methods.
 */
//$Id: TPZAgglomerateEl.cpp,v 1.53 2011-05-13 20:46:50 phil Exp $

#include "TPZAgglomerateEl.h"
#include "TPZInterfaceEl.h"
//#include "TPZEulerConsLaw.h"
//#include "TPZConservationLaw.h"
#include "pzdiscgal.h"
#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzquad.h"
#include "pzelmat.h"
#include "pzgraphel.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "pzgraphel1d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pztrigraph.h"
#include "pzgraphel.h"
#include "pzerror.h"
#include "pzmeshid.h"
#include "tpzagglomeratemesh.h"

using namespace std;

TPZAgglomerateElement::TPZAgglomerateElement(int nummat,int &index,TPZCompMesh &aggcmesh,TPZCompMesh *finemesh) :
TPZCompElDisc(aggcmesh,index) {
	
	/**
	 * o algomerado aponta para nulo mas o elemento computacional
	 * que ele agrupa aponta para o geom�rico original
	 * a copia do material na nova malha (malha aglomerada) neste
	 * ponto j�devia existir
	 */
	fMotherMesh = finemesh;
	
	//<!> Initialized as bignumber because I expect the radius will be lesser than that.
	fInnerRadius = TPZMaterial::gBigNumber;
	
	//<!>It is set up as zero. It must be initialized after.
	fNFaces = 0;
	
	fMaterialId = nummat;
	CreateMidSideConnect();
}

TPZAgglomerateElement::TPZAgglomerateElement() : TPZCompElDisc(), fIndexes(),
fMotherMesh(0),fInnerRadius(-1.),fNFaces(-1),fMaterialId(999)
{
}

void TPZAgglomerateElement::AddSubElementIndex(TPZCompMesh *aggcmesh,int subel,int father){
	
	TPZAgglomerateElement *agg = dynamic_cast<TPZAgglomerateElement *>(aggcmesh->ElementVec()[father]);
	if(!agg){
		for(int i=0;i<20;i++)
			PZError << "TPZAgglomerateElement::AddSubElementIndex null element\n";
		return;
	}
	agg->fIndexes.Push(subel);
	//a partir daqui os sub-elementos conhecem o aglomerado
	//depois a malha fina recupera as refer�cias originais
	//  agg->SetReference2(agg->NIndexes()-1);
}

void TPZAgglomerateElement::InitializeElement() {
	
	int indsize = NIndexes();
	//verificar se os materiais dos sub-elementos s� iguais
	int i,maxdeg = -1,mat,mat2;//bc degree
	mat = SubElement(0)->Material()->Id();
	for(i=1;i<indsize;i++){
		mat2 = SubElement(i)->Material()->Id();
		if(mat2 != mat){
			for(int k=0;k<10;k++)
				PZError << "TPZAgglomerateElement::TPZAgglomerateElement data error, distinct material\n";
			DebugStop();
		}
	}
	if(indsize < 1){
		PZError << "TPZAgglomerateElement::TPZAgglomerateElement empty list\n";
		return;
	}
	//tomando o grau como o m�imo grau dos sub-elementos
	for(i=0;i<indsize;i++){
		int deg = dynamic_cast<TPZCompElDisc *>(SubElement(i))->Degree();
		if(deg > maxdeg) maxdeg = deg;
	}
	SetDegree(maxdeg);
	CenterPoint();
	REAL cons = NormalizeConst();
	SetConstC(cons);
	//caso o conjunto de elementos sendo aglomerado preenchen totalmente
	//um nico elemento geom�rico esse ser�a referencia
	// ISTO E PURA PERFUMARIA, O elemento deve funcionar sem isto
	//  TPZGeoEl *ref = CalculateReference();
	//  TPZCompElDisc::SetReference(ref);
	//  if(ref) ref->SetReference(this);
}

void TPZAgglomerateElement::AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight){
	//acumula as regras de integra�o dos elementos aglomerados (descont�uos)
	int nsubs = NIndexes(),i;
	for(i=0; i<nsubs; i++){
		TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(SubElement(i));
		if(!disc){
			PZError << "TPZAgglomerateElement::AccumulateIntegrationRule, null index element\n";
			DebugStop();
		}
		//acumule integration rule
		disc->AccumulateIntegrationRule(degree,point,weight);
	}
}

void TPZAgglomerateElement::CenterPoint(){
	
	int indsize = NIndexes(),i;
	TPZVec<REAL> volumes(indsize);
	for(i=0;i<indsize;i++){
		//    TPZCompEl *cel = FineElement(i);
		TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(SubElement(i));
		if(!disc) {
			volumes[i] = 0.;
		} else {
			volumes[i] = disc->VolumeOfEl();
		}
	}
	REAL voltot = 0.0;
	for(i=0;i<indsize;i++) voltot += volumes[i];
	REAL centx=0.0,centy=0.0,centz=0.0;
	for(i=0;i<indsize;i++){
		//o decont�uo tem fCenterPoint
		//    TPZCompElDisc *cel = dynamic_cast<TPZCompElDisc *>(FineElement(i));
		TPZCompElDisc *cel = dynamic_cast<TPZCompElDisc *>(SubElement(i));
		if(!cel)  continue;
		centx += cel->CenterPoint(0) * volumes[i];
		centy += cel->CenterPoint(1) * volumes[i];
		centz += cel->CenterPoint(2) * volumes[i];
	}
	SetCenterPoint(0,centx/voltot);
	SetCenterPoint(1,centy/voltot);
	SetCenterPoint(2,centz/voltot);
}

void TPZAgglomerateElement::CenterPoint(TPZVec<REAL> &center){
	for(int i = 0; i < 3; i++)
		center[i] = fCenterPoint[i];
}

REAL TPZAgglomerateElement::VolumeOfEl(){
	
	/**
	 * um elemento aglomerado �formado de grupos de sub-elementos aglomerados
	 * ou de sub-elementos descont�uos n� aglomerados
	 */
	REAL volume = 0.0;
	int nindex = NIndexes(),i;
	for(i=0;i<nindex;i++){
		TPZCompEl *cel = dynamic_cast<TPZCompElDisc *>(SubElement(i));
		if(!cel) continue;
		volume += cel->VolumeOfEl();
	}
	return volume;
}


void TPZAgglomerateElement::CalcResidual(TPZFMatrix &Rhs,TPZCompElDisc *el){
	
	PZError << "TPZAgglomerateElement::CalcResidual DEVE SER IMPLEMENTADO";
}

void TPZAgglomerateElement::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){
	
	if(Reference()) return TPZCompElDisc::CalcStiff(ek,ef);
	
	if(!Material()){
		cout << "TPZCompElDisc::CalcStiff : no material for this element\n";
		ek.Reset();
		ef.Reset();
		return;
	}
	int ncon = NConnects();
	int dim = Dimension();
	int nstate = Material()->NStateVariables();
	int nshape = NShapeF();
	TPZBlock &block = Mesh()->Block();
	TPZFMatrix &MeshSol = Mesh()->Solution();
    int numbersol = MeshSol.Cols();
	int numeq = nshape * nstate;
	
	
	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,1);
	if(ncon){//pode ser no m�imo ncon = 1
		ek.fBlock.SetNBlocks(ncon);
		ef.fBlock.SetNBlocks(ncon);
		ek.fBlock.Set(0,NShapeF()*nstate);
		ef.fBlock.Set(0,NShapeF()*nstate);
	}
	ek.fConnect.Resize(ncon);
	ef.fConnect.Resize(ncon);
	for(int i=0;i<ncon;i++){
		(ef.fConnect)[i] = ConnectIndex(i);
		(ek.fConnect)[i] = ConnectIndex(i);
	}
	if(ncon==0) return;//elemento CC no passa
	REAL weight;
	int integ = 2*Degree();
	TPZStack<REAL> points,weights;
	//acumula as regras dos obtidos por aglomera�o
	AccumulateIntegrationRule(integ,points,weights);//integra fi*fj
	int ip,npoints = weights.NElements();
	
	TPZMaterialData data;
	this->InitMaterialData(data);
	for(ip=0;ip<npoints;ip++){
		data.x[0] = points[3*ip];
		data.x[1] = points[3*ip+1];
		data.x[2] = points[3*ip+2];
		weight = weights[ip];
		ShapeX(data.x,data.phi,data.dphix);
		//solucao da iteracao anterior
        for (int is=0; is<numbersol; is++) {
            data.sol[is].Fill(0.);
            data.dsol[is].Zero();
        }
		for(int in=0; in<ncon; in++) {
			TPZConnect *df = &Connect(in);
			int dfseq = df->SequenceNumber();
			int dfvar = block.Size(dfseq);
			int pos = block.Position(dfseq);
			int iv = 0,d;
			for(int jn=0; jn<dfvar; jn++) {
                for (int is=0; is<numbersol; is++) {
                    data.sol[is][iv%nstate] += data.phi(iv/nstate,0)*MeshSol(pos+jn,is);
                    for(d=0; d<dim; d++) data.dsol[is](d,iv%nstate) += data.dphix(d,iv/nstate)*MeshSol(pos+jn,is);
                }
				iv++;
			}
		}
		Material()->Contribute(data,weight,ek.fMat,ef.fMat);
	}
}


TPZCompEl *TPZAgglomerateElement::SubElement(int sub) const{
	
	int nsubs = NIndexes();
	
#ifdef DEBUG
	if(sub < 0  || sub > nsubs){
		PZError << "TPZAgglomerateElement::SubElement sub-element out of range\n";
		return NULL;
	}
#endif
	
	return fMotherMesh->ElementVec()[fIndexes[sub]];
}


REAL TPZAgglomerateElement::NormalizeConst() {
	//maior distancia entre o ponto interior e os v�tices do elemento
	REAL maxsub,maxall = -1.0;
	int nindex = NIndexes(),i;
	for(i=0;i<nindex;i++){
		TPZCompElDisc *cel = dynamic_cast<TPZCompElDisc *>(SubElement(i));
		if(!cel) continue;
		maxsub = cel->NormalizeConst();
		if(maxall < maxsub) maxall = maxsub;
	}
	return maxall;
}


int TPZAgglomerateElement::CreateMidSideConnect() {	
	if(!Material())
	{
		PZError << "\nTPZCompElDisc::CreateMidSideConnect Material nulo\n";
		DebugStop();
	}
	
	TPZStack<TPZCompElSide> list;
	int dim = Dimension();
	int InterfaceDimension = this->Material()->Dimension() - 1;
	
	if(dim == InterfaceDimension){
		// o atual �um elemento BC
		SetConnectIndex(0, -1);
		SetDegree(-1);//=> nshape = 0
	} else {
		int nvar = Material()->NStateVariables();
        int nshape = NShapeF();
        int order = 0;
		int newnodeindex = Mesh()->AllocateNewConnect(nshape, nvar, order);
		SetConnectIndex(0,newnodeindex);
		TPZConnect &newnod = Mesh()->ConnectVec()[newnodeindex];
		int seqnum = newnod.SequenceNumber();
		Mesh()->Block().Set(seqnum,nvar*nshape);
		Mesh()->ConnectVec()[ConnectIndex()].IncrementElConnected();
	}

	return ConnectIndex();
}

int TPZAgglomerateElement::Dimension() const {
	return this->Material()->Dimension();
}

void TPZAgglomerateElement::Print(std::ostream &out) const {
	
	out << "\nTPZAgglomerateElement element : \n";
	out << "\tComputational mesh : " << fMotherMesh << endl;
	out << "\tAgglomerate elements indexes : ";
	int naggel = NIndexes(),i;
	for(i=0;i<naggel;i++) out << fIndexes[i] << " ";
	out << "\n\tMaterial id : " << fMotherMesh->ElementVec()[fIndexes[0]]->Material()->Id() << endl
	<< "\tDegrau of interpolation : " <<  Degree() << endl
	<< "\tConnect index : " << ConnectIndex() << endl
	<< "\tNormalizing constant : " << ConstC() << endl
	<< "\tCenter point of the element : ";
	for(i=0;i<2;i++) out << TPZCompElDisc::CenterPoint(i) << " , ";
	out << TPZCompElDisc::CenterPoint(i) << endl;
}

void TPZAgglomerateElement::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
	if(!Reference()) return;
	int mat = Material()->Id();
	int nsides = NSides();
	
	if(dimension == 2 && mat > 0){
		if(nsides == 9){
			new TPZGraphElQ2dd(this,&grmesh);
			return;
		}
		if(nsides == 7){
			new TPZGraphElTd(this,&grmesh);
			return;
		}
	}
	if(dimension == 3 && mat > 0){
		new TPZGraphElQ3dd(this,&grmesh);
	}
	if(dimension == 1 && mat > 0){
		new TPZGraphEl1dd(this,&grmesh);
	}
}

int TPZAgglomerateElement::NSides(){
	
	if(Reference()) return Reference()->NSides();
	return -1;
}

void TPZAgglomerateElement::ListOfDiscEl(TPZStack<TPZCompEl *> &elvec){
	//agrupa na lista todos os elementos discontinuos
	//agrupados no elemento aglomerado atual
	int indsize = NIndexes(),i;
	for(i=0;i<indsize;i++){
		//    TPZCompEl *cel = FineElement(i);
		TPZCompElDisc *cel = dynamic_cast<TPZCompElDisc *>(SubElement(i));
		if(!cel) continue;
		if(cel->Type() == EAgglomerate){
			dynamic_cast<TPZAgglomerateElement *>(cel)->ListOfDiscEl(elvec);
		} else if(cel->Type() == EDiscontinuous){
			elvec.Push(cel);//guarda todos os descont�uos
		} else {
			PZError << "TPZAgglomerateElement::NSides unknow type element\n";
		}
	}
}

//#warning "TPZAgglomerateElement::IndexesDiscSubEls is inconsistent"
void TPZAgglomerateElement::IndexesDiscSubEls(TPZStack<int> &elvec){
	//agrupa na lista todos os indexes dos elementos discontinuos
	//agrupados no elemento aglomerado atual
	int indsize = NIndexes(),i;
	for(i=0;i<indsize;i++){
		//    TPZCompEl *cel = FineElement(i);
		TPZCompEl *cel = SubElement(i);
		if(cel->Type() == EAgglomerate){
			dynamic_cast<TPZAgglomerateElement *>(cel)->IndexesDiscSubEls(elvec);
		} else if(cel->Type() == EDiscontinuous){
			elvec.Push(cel->Index());//guarda todos os descont�uos
		} else {
			PZError << "TPZAgglomerateElement::NSides unknow type element\n";
		}
	}
}

/*
 void TPZAgglomerateElement::CreateMaterialCopy(TPZCompMesh &aggcmesh){
 
 //criando copias dos materiais
 int nmats = fMotherMesh->MaterialVec().NElements(),i;
 for(i=0;i<nmats;i++){//achando material de volume
 TPZMaterial *mat = fMotherMesh->MaterialVec()[i];
 if(!mat) continue;
 if(mat->Id() > 0){
 //       if( !strcmp(mat->Name(),"TPZEulerConsLaw") ){
 // 	TPZEulerConsLaw *euler = dynamic_cast<TPZEulerConsLaw *>(mat);
 // 	mat = euler->NewMaterial();//mat = new TPZEulerConsLaw(*euler);//copia
 // 	aggcmesh.InsertMaterialObject(mat);
 //       }
 TPZDiscontinuousGalerkin * consmat = dynamic_cast<TPZDiscontinuousGalerkin *>(mat);
 if ( consmat ){
 mat = consmat->NewMaterial();
 aggcmesh.InsertMaterialObject(mat);
 } else {
 PZError << "TPZAgglomerateElement::CreateMaterialCopy material not defined, implement now (Tiago)!\n";
 }
 }
 }
 for(i=0;i<nmats;i++){//achando material de CC
 TPZMaterial *mat = fMotherMesh->MaterialVec()[i];
 if(!mat) continue;
 if(mat->Id() < 0){//CC: id < 0
 TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
 if(!bc) PZError << "TPZCompElDisc::CreateAgglomerateMesh null bc material\n";
 int nummat = bc->Material()->Id();//mat. vol.
 TPZMaterial *material = aggcmesh.FindMaterial(nummat);
 if(!material) PZError << "TPZCompElDisc::CreateAgglomerateMesh volume material not exists"
 << " (implement Tiago)\n";
 TPZBndCond *copy = new TPZBndCond(*bc,material);
 aggcmesh.InsertMaterialObject(copy);
 }
 }
 }
 */

/*
 TPZGeoEl *TPZAgglomerateElement::CalculateReference(){
 
 //return TPZCompElDisc::Reference();
 //porenquanto interfaces n� s� aglomeradas
 //if(Dimension() == gInterfaceDimension) return NULL;
 
 TPZStack<TPZCompEl *> elvec;
 ListOfDiscEl(elvec);
 int size = elvec.NElements(),i,nlevels = 1;
 TPZGeoEl *ref0 = elvec[0]->Reference();
 if(size == 1) return ref0;
 
 TPZGeoEl *fat = FatherN(ref0,nlevels);
 while(fat){
 for(i=1;i<size;i++){
 TPZGeoEl *ref = elvec[i]->Reference();
 if( FatherN(ref,nlevels) != fat ) break;
 }
 if( i < size ){
 nlevels++;
 fat = FatherN(ref0,nlevels);
 } else {
 break;//o pai maior foi achado
 }
 }
 if( fat && size == NSubCompEl(fat) ) return fat;
 return NULL;
 }
 
 int TPZAgglomerateElement::NSubsOfLevels(TPZGeoEl *father,int nlevels){
 
 //quantos sub-elementos father cont� at�nlevels n�eis abaixo dele
 if(!father) return 0;
 if(nlevels == 1) return father->NSubElements();
 int numberels = 0,i,lev = 1;
 while(lev < nlevels){
 int nsubs = father->NSubElements();
 numberels = nsubs;
 for(i=0;i<nsubs;i++){
 TPZGeoEl *sub = father->SubElement(i);
 if(!sub->Reference()) NSubsOfLevels(sub,lev);//s�computacionais
 numberels += NSubsOfLevels(sub,nlevels-1);
 }
 lev++;
 }
 return numberels;
 }
 
 int TPZAgglomerateElement::NSubCompEl(TPZGeoEl *father){
 
 //quantos sub-elementos computacionais father aglomera
 if(!father) return 0;
 int numberels = 0,i,nsubs = father->NSubElements();
 for(i=0;i<nsubs;i++){
 TPZGeoEl *sub = father->SubElement(i);
 if(!sub) return 0;
 if(!sub->Reference()){
 numberels += NSubCompEl(sub);//s�computacionais
 } else {
 numberels++;
 }
 }
 return numberels;
 }
 
 TPZGeoEl *TPZAgglomerateElement::FatherN(TPZGeoEl *sub,int n){
 
 //procura-se o ancestral n n�eis acima de sub
 if(!sub) return NULL;
 if(n == 0) return sub;
 int niv = 1;
 TPZGeoEl *fat;
 while( niv < (n+1) ){
 fat = sub->Father();
 if(!fat) return NULL;
 sub = fat;
 niv++;
 }
 return fat;
 }
 */

//int Level(TPZGeoEl *gel);
void TPZAgglomerateElement::ListOfGroupings(TPZCompMesh *finemesh,TPZVec<int> &accumlist,
											int nivel,int &numaggl,int dim){
	//todo elemento de volume deve ser agrupado nem que for para um s�elemento
	finemesh->SetDimModel(dim);
	cout << "\nTPZAgglomerateElement::ListOfGroupings para malha 2D\n";
	cout << "Este metodo somente funciona para agrupar elementos descontinuos \n";
	int nel = finemesh->NElements(),i;
	//n� todo index �sub-elemento
	accumlist.Resize(nel,-1);
	int mdim = finemesh->Dimension();
	int sdim = mdim - 1;//dimens� superficial
	for(i=0;i<nel;i++){
		TPZCompEl *cel = finemesh->ElementVec()[i];
		if(!cel) continue;
		TPZGeoEl *gel = cel->Reference();
		if(!gel) {
			cout << "TPZAgglomerateElement::ListOfGroupings nao funciona para esta malha\n";
			return;
		}
		int type = cel->Type(),eldim = cel->Dimension();
		//agrupando elementos computacionais de volume
		if(type == EInterface) continue;//interface ser�clonada
		if(eldim == sdim) continue;//discontinuo BC ser�clonado
		TPZGeoEl *father = gel->Father();
		if(!father) continue;
		while(father && father->Level() != nivel) father = father->Father();//nova
		if (!father) continue;//nova
		//if(Level(father) != nivel) continue;//antiga
		int fatid = father->Id();
		accumlist[i] = fatid;
	}
	//reordena a lista por ordem crescente do pai
	TPZVec<int> list(accumlist);
	int j;
	for(i=0;i<nel;i++){
		for(j=i+1;j<nel;j++){
			if(list[i] > list[j]){
				int aux = list[i];
				list[i] = list[j];
				list[j] = aux;
			}
		}
	}
	//conta o nmero de elementos obtidos por aglomera�
	numaggl = 0;
	int act2 = -1;
	for(i=0;i<nel;i++){
		int act = list[i];
		if(act == act2) continue;
		for(j=i+1;j<nel;j++){
			int next = list[j];
			if(next != act2){
				numaggl++;
				j = nel;
				act2 = act;
			}
		}
	}
	//reformula o index do pai de 0 a nmax
	TPZVec<int> newlist(accumlist);
	int newfat = 0;
	for(i=0;i<nel;i++){
		int fatid1 = newlist[i];
		if(fatid1 < 0) continue;
		accumlist[i] = newfat;
		newlist[i] = -1;
		for(j=i+1;j<nel;j++){
			int fatid2 = newlist[j];
			if(fatid2 == fatid1){
				accumlist[j] = newfat;
				newlist[j] = -1;
			}
		}
		newfat++;
	}
	if(newfat != numaggl) cout << "TPZAgglomerateElement::ListOfGroupings nmero de pais n� confere\n";
	if(!newfat && !numaggl) cout << "TPZAgglomerateElement::ListOfGroupings lista de elementos aglomerados vacia\n";
}

void TPZAgglomerateElement::Print(TPZStack<int> &listindex){
	
	cout << "TPZAgglomerateElement::Print agrupamento de indexes: saida AGRUPAMENTO.out\n";
	ofstream out("AGRUPAMENTO.out");
	int size = listindex.NElements(),i,father = 0;
	out << "\n\t\t\t* * * AGRUPAMENTO DE INDEXES * * *\n\n";
	int exists = 0;
	for(i=0;i<size;i++){
		if(listindex[i] == father){
			out << i << ' ';
			exists = 1;
		}
		if( (i+1) == size && exists ){
			out << "-> father = " << father << endl;
			i = -1;
			father++;
			exists = 0;
		}
	}
}

/*
 int Level(TPZGeoEl *gel){
 //retorna o n�el do elemento gel
 if(!gel) return -1;
 TPZGeoEl *fat = gel->Father();
 if(!fat) return 0;
 int niv = 0;
 while(fat){
 fat = fat->Father();
 niv++;
 }
 return niv;
 }
 */

void TPZAgglomerateElement::ProjectSolution(TPZFMatrix &projectsol){
	
	//projeta a solu�o dos elementos contidos na aglomera�o
	//no elemento por eles aglomerado
	
	int dimension = Dimension();
	int aggmatsize = NShapeF();
	int nvar = Material()->NStateVariables();
	TPZFMatrix aggmat(aggmatsize,aggmatsize,0.);
	TPZFMatrix loadvec(aggmatsize,nvar,0.);
	//verificar que o grau do grande �<= que o grau de cada um dos pequenos ???
	TPZStack<int> elvec;
	//  IndexesDiscSubEls(elvec);
	//int size = elvec.NElements(),i;
	int size = NIndexes(),i;
	int mindegree = 1000,coarsedeg = Degree(),maxdegree = 0;
	for(i=0;i<size;i++){
		TPZCompElDisc *disc =
		//      dynamic_cast<TPZCompElDisc *>(FineElement(i));
		dynamic_cast<TPZCompElDisc *>(SubElement(i));
		int degree = disc->Degree();
		if(mindegree > degree) mindegree = degree;
		if(maxdegree < degree) maxdegree = degree;
	}
	if(coarsedeg > mindegree){
		//PZError << "TPZAgglomerateElement::RestrictionOperator incompatible degrees\n";
		//return;
		cout << "TPZAgglomerateElement::RestrictionOperator MUDANDO A ORDEM DO ELEMENTO\n";
		SetDegree(mindegree);
		coarsedeg = mindegree;
	}
	//para integrar uh * fi sobre cada o sub-elemento:
	//fi base do grande, uh solu�o sobre os pequenos
	//   TPZIntPoints *intrule =
	//     ref->CreateSideIntegrationRule(ref->NSides()-1,maxdegree + coarsedeg);
	//eventualmente pode criar uma regra para cada sub-elemento dentro do la�
	//de integra�o para reduzir o nmero de pontos caso os grus s� distintos
	TPZFMatrix aggphix(aggmatsize,1);
	TPZFMatrix aggdphix(dimension,aggmatsize);
	TPZStack<REAL> accintpoint, accweight;
	TPZFMatrix jacobian(dimension,dimension),jacinv(dimension,dimension);
	TPZFMatrix axes(3,3,0.);
	TPZManVector<REAL,3> x(3,0.);
	TPZVec<REAL> uh(nvar);
	REAL weight;
	int in,jn,kn,ip,ind;
	TPZCompElDisc *finedisc;
	
	for(ind=0;ind<size;ind++){
		//    finedisc = dynamic_cast<TPZCompElDisc *>(FineElement(ind));
		finedisc = dynamic_cast<TPZCompElDisc *>(SubElement(ind));
		finedisc->AccumulateIntegrationRule(maxdegree+coarsedeg,accintpoint,accweight);
		int npoints = accweight.NElements();
		
		for(ip=0;ip<npoints;ip++){
			for(in=0; in<3; in++) x[in] = accintpoint[3*ip+in];
			weight = accweight[ip];
			ShapeX(x,aggphix,aggdphix);
			finedisc->SolutionX(x,uh);
			//projetando a solu�o fina no elemento aglomerado
			//a soma dos detjac dos sub-elementos d�o detjac do grande
			//a geometria do grande �a soma das geometrias dos pequenos
			for(in=0; in<aggmatsize; in++) {
				for(jn=0; jn<aggmatsize; jn++) {
					//ordem do grande <= a menor ordem dos pequenos
					// => a regra  dos pequenos integra ok
					aggmat(in,jn) += weight*aggphix(in,0)*aggphix(jn,0);
				}
				//a soma das regras dos pequenos cobre a geometria do grande
				for(kn=0; kn<nvar; kn++) {
					loadvec(in,kn) += weight*aggphix(in,0)*uh[kn];
				}
			}
		}//fim for ip
	}
	//achando a solu�o restringida
	aggmat.SolveDirect(loadvec,ELDLt);//ELU
	//transferindo a solu�o obtida por restri�o
	int iv=0,d;
	TPZBlock &block = Mesh()->Block();
	TPZConnect *df = &Connect(0);
	int dfseq = df->SequenceNumber();
	int dfvar = block.Size(dfseq);
	int pos   = block.Position(dfseq);
	for(d=0; d<dfvar; d++) {
		//block(dfseq,0,d,0) = loadvec(iv/nvar,iv%nvar);
		projectsol(pos+d,0) = loadvec(iv/nvar,iv%nvar);
		iv++;
	}
}

REAL TPZAgglomerateElement::LesserEdgeOfEl(){
	
	
	int j,l,nvertices;
	TPZVec<REAL> point0(3),point1(3);
	//TPZGeoNode *node0,*node1;
	//procura-se a maior distancia entre dois nodos do conjunto de elementos
	//em cada dire�o X, Y ou Z, retorna-se a menor dessas 3 distancias
	TPZStack<TPZGeoNode *> nodes;
	AccumulateVertices(nodes);
	nvertices = nodes.NElements();
	REAL maxX=-1.,maxY=-1.,maxZ=-1.;
	REAL distX,distY,distZ;
	for(l=0;l<nvertices;l++){
		for(j=l+1;j<nvertices;j++){
			//node0 = nodes[l];
			//node1 = nodes[j];
			distX = fabs(nodes[l]->Coord(0) - nodes[j]->Coord(0));
			distY = fabs(nodes[l]->Coord(1) - nodes[j]->Coord(1));
			distZ = fabs(nodes[l]->Coord(2) - nodes[j]->Coord(2));
			if(distX > maxX) maxX = distX;
			if(distY > maxY) maxY = distY;
			if(distZ > maxZ) maxZ = distZ;
		}
	}
	if(maxX < 1.e-8 && maxY < 1.e-8  && maxZ < 1.e-8) {
		PZError << "TPZCompElDisc::LesserEdgeOfEl degenerate element\n";
	}
	if(maxX < 1.e-8) maxX = 1000000.0;
	if(maxY < 1.e-8) maxY = 1000000.0;
	if(maxZ < 1.e-8) maxZ = 1000000.0;
	maxX = maxX < maxY ? maxX : maxY;
	maxX = maxX < maxZ ? maxX : maxZ;
	return maxX;
}


/*
 if(0){
 //procurando a menor distancia entre dois nodos de qualquer um dos sub-elementos???!!!
 for(i=0;i<size;i++){
 TPZCompEl *comp = This->MotherMesh()->ElementVec()[elvec[i]];
 TPZGeoEl *geo = comp->Reference();
 if(!geo) PZError <<  "TPZAgglomerateElement::Solution null reference\n";
 int nvertices = geo->NNodes();
 for(l=0;l<nvertices;l++){
 for(j=l+1;j<nvertices;j++){
 node0 = geo->NodePtr(l);
 node1 = geo->NodePtr(j);
 for(k=0;k<3;k++){
 point0[k] = node0->Coord(k);
 point1[k] = node1->Coord(k);
 }
 dist = TPZGeoEl::Distance(point0,point1);
 if(dist < mindist) mindist = dist;
 }
 }
 }
 }
 return mindist;
 */


void TPZAgglomerateElement::AccumulateVertices(TPZStack<TPZGeoNode *> &nodes) {
	int nsubs = fIndexes.NElements();
	int in;
	for(in=0; in<nsubs; in++) {
		TPZCompEl *cel = fMotherMesh->ElementVec()[fIndexes[in]];
		if(!cel) continue;
		TPZCompElDisc *cdisc = dynamic_cast<TPZCompElDisc *>(cel);
		if(!cdisc) continue;
		cdisc->AccumulateVertices(nodes);
	}
}



TPZAgglomerateMesh *TPZAgglomerateElement::CreateAgglomerateMesh(TPZCompMesh *finemesh, TPZVec<int> &accumlist,int numaggl){
	
	/**
	 * somente s� aglomerados elementos de volume
	 * elementos de interface n� s� aglomerados, cada interface deve conhecer
	 * o elemento aglomerado esquerdo e direito
	 * elementos fantasmas ou BC n� s� aglomerados, a cada elemento BC
	 * corresponde um elemento interface de igual tamanho ou n�el
	 * todo elemento interface e todo elemento BC deve ser clonado para a malha aglomerada
	 * todo elemento de volume deve ter associado um agrupamento podendo ser um nico elemento
	 * a posi� K de accumlist indica o index K do elemento computacional que ser�acumulado,
	 * o inteiro guardado nessa posi� indica o elemento ao qual sera
	 * aglomerado, assim si accumlist[8] = 4 ent� o elemento computacional
	 * de index 8 ser�agrupado para formar o elemento aglomerado de index 4
	 */
	int nlist = accumlist.NElements();
	if(numaggl < 1 || nlist < 2){
		PZError << "TPZCompElDisc::CreateAgglomerateMesh number agglomerate elements"
		<< " out of range\nMALHA AGGLOMERADA N� CRIADA\n";
		//    aggmesh = new TPZCompMesh(*finemesh);
		return 0;
	}
	//TPZFlowCompMesh aggmesh(finemesh->Reference());
	int index,nel = finemesh->NElements(),nummat,size = accumlist.NElements(),i;
	TPZAgglomerateMesh *aggmesh = new TPZAgglomerateMesh(finemesh);
	//copiando materiais para nova malha
	//finemesh->fMaterials eh copiado para aggmesh
	finemesh->CopyMaterials(*aggmesh);
	
	TPZVec<int> IdElNewMesh(accumlist.NElements());
	IdElNewMesh.Fill(-1);
	
	//criando agrupamentos iniciais
	for(i=0;i<numaggl;i++){
		int k = 0;
		//o elemento de index k tal que accumlist[k] == i vai aglomerar no elemento de index/id i
		while( accumlist[k] != i && k < size) k++;
		
#ifdef DEBUG
		if(k == size){
			PZError << "TPZCompElDisc::CreateAgglomerateMesh material not found\n";
			DebugStop();
		}
#endif
		
		nummat = finemesh->ElementVec()[k]->Material()->Id();
		//criando elemento de index/id i e inserindo na malha aggmesh
		new TPZAgglomerateElement(nummat,index,*aggmesh,finemesh);
		IdElNewMesh[k] = index;
	}
	
	//inicializando os aglomerados e clonando elementos BC
	int meshdim = finemesh->Dimension(),father,type,eldim;
	int surfacedim = meshdim-1;
	for(i=0;i<nel;i++){
		TPZCompEl *cel = finemesh->ElementVec()[i];
		if(!cel) continue;
		father = accumlist[i];
		type = cel->Type();
		eldim = cel->Dimension();
		if(type == EDiscontinuous || type == EAgglomerate){//elemento descont�uo de volume ou aglomerado
			if(eldim == meshdim){
				if(father == -1){
					PZError << "TPZCompElDisc::CreateAgglomerateMesh element null father\n";
					if(aggmesh) delete aggmesh;
					return 0;
				}
				//incorporando o index do sub-elemento
				//o elemento de index i vai para o pai father = accumlist[i]
				TPZAgglomerateElement::AddSubElementIndex(aggmesh,i,father);
				IdElNewMesh[i] = father;
			} else if(eldim == surfacedim){
				//clonando o elemento descont�uo BC, o clone existira na malha aglomerada
				TPZCompEl * AggCell = dynamic_cast<TPZCompElDisc *>(cel)->Clone(*aggmesh,index);
				IdElNewMesh[i] = AggCell->Index();
			}
		}
	}
	
	//agora os geom�ricos apontam para os respectivos aglomerados
	//computacionais recuperaram as refer�cias com TPZCompMesh::LoadReferences()
	for(i=0;i<nel;i++){
		TPZCompEl *cel = finemesh->ElementVec()[i];
		if(!cel) continue;
		if(cel->Type() == EInterface){//elemento interface
			TPZInterfaceElement *interf = dynamic_cast<TPZInterfaceElement *>(cel);
			TPZCompEl *leftel = interf->LeftElement(),*rightel = interf->RightElement();
			int indleft = leftel->Index();
			int indright = rightel->Index();
			if(accumlist[indleft] == -1 && accumlist[indright] == -1){
				//no m�imo um elemento pode ser BC
				PZError << "TPZCompElDisc::CreateAgglomerateMesh data error\n";
				if(aggmesh) delete aggmesh;
				return 0;
			}
			int index;
			//interfaces com esquerdo e direito iguais n� s� clonadas
			if(accumlist[indleft] == accumlist[indright]) continue;
			
			//clonando o elemento interface: a interface existir�na malha aglomerada
			//leftel and rightel
			int leftid = interf->LeftElement()->Index();
			int leftside = interf->LeftElementSide().Side();
			leftid = IdElNewMesh[leftid];
#ifdef DEBUG
			if (leftid == -1){
				PZError <<  "\nLeftid cannot be -1." << endl;
			}
#endif
			TPZCompElDisc * leftagg = dynamic_cast<TPZCompElDisc*> (aggmesh->ElementVec()[leftid]) ;
			
			int rightid = interf->RightElement()->Index();
			int rightside = interf->RightElementSide().Side();
			rightid = IdElNewMesh[rightid];
#ifdef DEBUG
			if (rightid == -1){
				PZError <<  "\nRightid cannot be -1." << endl;
			}
#endif
			TPZCompElDisc * rightagg = dynamic_cast<TPZCompElDisc*> (aggmesh->ElementVec()[rightid]) ;
			
			TPZCompElSide leftcompelside(leftagg, leftside), rightcompelside(rightagg, rightside);
			interf->CloneInterface(*aggmesh,index, /*leftagg*/leftcompelside, /*rightagg*/rightcompelside);
			
		}
	}
	
	nel = aggmesh->ElementVec().NElements();
	
	//inizializando elementos aglomerados preenchendo os demais dados
	for(i=0;i<nel;i++){
		TPZCompEl *cel = aggmesh->ElementVec()[i];
		if(cel->Type() == EAgglomerate){//agglomerate element
			TPZAgglomerateElement *agg = dynamic_cast<TPZAgglomerateElement *>(cel);
			agg->InitializeElement();
			//it is set up as zero. It will be computed below.
			agg->SetNInterfaces(0);
		}
	}
	
	int dim = aggmesh->Dimension();
	//<!>Loop over all elements to compute for each interface the distance between interface's center point to volume's center point.
	//The lessest distance is adopted as the element's inner radius (fInnerRadius).
	//It is also computed, the number of interfaces of an agglomerate element (fNFaces).
	for(i = 0; i < nel; i++){
		
		TPZCompEl *cel = aggmesh->ElementVec()[i];
		
		TPZInterfaceElement * interface = dynamic_cast<TPZInterfaceElement * > (cel);
		
		if (interface) {
			
			TPZCompElDisc * LeftEl   = dynamic_cast<TPZCompElDisc*>( interface -> LeftElement()  );
			TPZCompElDisc * RightEl  = dynamic_cast<TPZCompElDisc*>( interface -> RightElement() );
			
			TPZManVector<REAL, 3> InterfaceCenter(3), RefInterfaceCenter(3), LeftCenter(3), RightCenter(3);
			
			LeftEl->CenterPoint(LeftCenter);
			RightEl->CenterPoint(RightCenter);
			
			
			REAL LeftDistance = -1., RightDistance = -1.;
			int nsides = interface->Reference()->NSides();
			
			
			for (int irib = 0; irib < nsides; irib++){
				REAL left = 0., right = 0.;
				//Getting the node coordinate
				interface->Reference()->CenterPoint(irib, RefInterfaceCenter);
				interface->Reference()->X(RefInterfaceCenter, InterfaceCenter);
				
				for (int isum = 0; isum < dim; isum++)
				{
					left  += (LeftCenter[isum] - InterfaceCenter[isum]) * (LeftCenter[isum] - InterfaceCenter[isum]);
					right += (RightCenter[isum] - InterfaceCenter[isum]) * (RightCenter[isum] - InterfaceCenter[isum]);
				}
				left  = sqrt(left);
				right = sqrt(right);
				if(left < LeftDistance || LeftDistance < 0.) LeftDistance = left;
				if(right < RightDistance || RightDistance < 0.) RightDistance = right;
			}
			
			
			//if the stored inner radius is bigger than the computed one, the computed one takes its place, because of we are computing the INNER radius.
			//only Agglomerate must be used here. CompElDisc have Reference()->ElementRadius()
			TPZAgglomerateElement *LeftAgg = dynamic_cast <TPZAgglomerateElement*> (LeftEl);
			
			if (LeftAgg){
				if ( LeftDistance < LeftAgg->InnerRadius2() )
					LeftAgg->SetInnerRadius(LeftDistance);
				LeftAgg ->SetNInterfaces(LeftAgg->NInterfaces() + 1);
			}
			
			TPZAgglomerateElement *RightAgg = dynamic_cast <TPZAgglomerateElement*> (RightEl);
			
			if (RightAgg){
				if ( RightDistance < RightAgg->InnerRadius2() )
					RightAgg->SetInnerRadius(RightDistance);
				RightAgg->SetNInterfaces(RightAgg->NInterfaces() + 1);
			}
			
			
		}//end of if
		
	}//end of loop over elements
	return aggmesh;
	
}//end of method

void TPZAgglomerateElement::ComputeNeighbours(TPZCompMesh *mesh, map<TPZCompElDisc *,set<TPZCompElDisc *> > &neighbour) {
	
	int nelem = mesh->NElements();
	int iel;
	//  int meshdim = mesh->Dimension();
	for(iel=0; iel<nelem; iel++) {
		TPZCompEl *cel = mesh->ElementVec()[iel];
		TPZInterfaceElement *inter = dynamic_cast<TPZInterfaceElement *>(cel);
		if(!inter) continue;
		TPZCompElDisc * LeftEl  =  dynamic_cast<TPZCompElDisc*>( inter->LeftElement() );
		TPZCompElDisc * RightEl =  dynamic_cast<TPZCompElDisc*>( inter->RightElement());
		if(LeftEl->Dimension() == RightEl->Dimension()) {
			neighbour[LeftEl].insert(RightEl);
			neighbour[RightEl].insert(LeftEl);
		}
	}
}


/**
 * returns the unique identifier for reading/writing objects to streams
 */
int TPZAgglomerateElement::ClassId() const
{
	return TPZAGGLOMERATEELID;
}

template class
TPZRestoreClass< TPZAgglomerateElement, TPZAGGLOMERATEELID>;

/**
 Save the element data to a stream
 */
void TPZAgglomerateElement::Write(TPZStream &buf, int withclassid)
{
	TPZCompElDisc::Write(buf,withclassid);
	TPZAgglomerateMesh *mesh = dynamic_cast<TPZAgglomerateMesh *> (Mesh());
	if(!mesh)
	{
		cout << "TPZAgglomerateEl::Write built from non agglomeratemesh at " << __FILE__ << ":" << __LINE__ << endl;
	}
	WriteObjects(buf,fIndexes);
	buf.Write(&fInnerRadius,1);
	buf.Write(&fNFaces,1);
}

/**
 Read the element data from a stream
 */
void TPZAgglomerateElement::Read(TPZStream &buf, void *context)
{
	TPZCompElDisc::Read(buf,context);
	TPZAgglomerateMesh *mesh = dynamic_cast<TPZAgglomerateMesh *> (Mesh());
	if(!mesh)
	{
		cout << "TPZAgglomerateEl::Read built from non agglomeratemesh at " << __FILE__ << ":" << __LINE__ << endl;
		fMotherMesh = 0;
	} else {
		fMotherMesh = mesh->FineMesh();
	}
	ReadObjects(buf,fIndexes);
	buf.Read(&fInnerRadius,1);
	buf.Read(&fNFaces,1);
}

