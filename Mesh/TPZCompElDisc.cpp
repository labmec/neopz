//$Id: TPZCompElDisc.cpp,v 1.45 2004-02-04 20:30:24 tiago Exp $

// -*- c++ -*- 

#include "pztransfer.h"
#include "pzelmat.h"
#include "pzelgc3d.h"
#include "pzelgt3d.h"
#include "pzelgpi3d.h"
#include "pzelgpr3d.h"
#include "pzelgt2d.h"
#include "pzelgq2d.h"
#include "pzmatrix.h"
#include "pzelmat.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzerror.h"
#include "pzconnect.h"
#include "pzshapequad.h"
#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzmat1dlin.h"
#include "pztempmat.h"
#include "pzmanvector.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"
#include "TPZConservationLaw.h"
#include "TPZEulerConsLaw.h"
#include "pzgraphel.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "pzgraphel1d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pztrigraph.h"
#include "pzgraphel.h"
//#include "TPZFlowCMesh.h"

#include "time.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include <math.h>
#include <stdio.h>

int TPZCompElDisc::gDegree = 0;

//construtor do elemento aglomerado
TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,int &index) :
		TPZCompEl(mesh,index), fCenterPoint(3) {
  fDegree = gDegree;
  fReference = NULL;
  fMaterial = NULL;
  fShapefunctionType = TPZShapeDisc::EOrdemTotal;
}

TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh, const TPZCompElDisc &copy) :
		TPZCompEl(mesh,copy), fCenterPoint(copy.fCenterPoint) {
  fDegree = copy.fDegree;
  fShapefunctionType = copy.fShapefunctionType;
  fReference = copy.fReference;
  TPZMaterial *mat = copy.Material();
  if(mat) {
    int materialid = mat->Id();
    fMaterial = mesh.FindMaterial(materialid);
  } else {
    fMaterial = 0;
  }
}

TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh, const TPZCompElDisc &copy,int &index) :
		TPZCompEl(mesh,copy,index), fCenterPoint(copy.fCenterPoint) {
  fDegree = copy.fDegree;
  fShapefunctionType = copy.fShapefunctionType;
  fReference = copy.fReference;
  //criando nova malha computacional
  fReference->SetReference(this);
  TPZMaterial *mat = copy.Material();
  if(mat) {
    int materialid = mat->Id();
    fMaterial = mesh.FindMaterial(materialid);
  } else {
    fMaterial = 0;//não deveria acontecer
  }
  fConstC = copy.fConstC;
  CreateMidSideConnect();
  //as interfaces foram clonadas
}

//construtor do elemento descontínuo
TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,TPZGeoEl *ref,int &index) :
		TPZCompEl(mesh,index), fCenterPoint(3) {
  fDegree = gDegree;
  fShapefunctionType = TPZShapeDisc::EOrdemTotal;
  switch(ref->Type()) {
  case EQuadrilateral:
  case ECube:
  case EPrisma:
    fShapefunctionType = TPZShapeDisc::ETensorial;
  }
  fReference = ref;
  ref->SetReference(this);
  //fMesh = &mesh;
  int materialid = ref->MaterialId();
  fMaterial = mesh.FindMaterial(materialid);
  CreateMidSideConnect();
  ref->CenterPoint(ref->NSides()-1,fCenterPoint);
  TPZVec<REAL> csi(fCenterPoint);
  ref->X(csi,fCenterPoint);
  fConstC = NormalizeConst();
  //criando os elementos interface
  CreateInterfaces();
}

void TPZCompElDisc::CreateInterfaces(){
  //não verifica-se caso o elemento de contorno
  //é maior em tamanho que o interface associado
  //caso AdjustBoundaryElement não for aplicado
  //a malha é criada consistentemente
  int nsides = fReference->NSides();
  int InterfaceDimension = fMaterial->Dimension() - 1;
  int side;
  nsides--;//last face
  for(side=nsides;side>=0;side--){
    if(fReference->SideDimension(side) != InterfaceDimension) continue;
    TPZCompElSide thisside(this,side);
    if(ExistsInterface(thisside.Reference())) {
      cout << "TPZCompElDisc::CreateInterface inconsistent: interface already exists\n";
      continue;
    }
    TPZStack<TPZCompElSide> highlist;
    thisside.HigherLevelElementList(highlist,0,1);
    //a interface se cria uma vez só: quando existem ambos 
    //elementos esquerdo e direito (computacionais)
    if(!highlist.NElements()) {
      CreateInterface(side);//só tem iguais ou grande => pode criar a interface
    } else {
      int ns = highlist.NElements();
      int is;
      for(is=0; is<ns; is++) {//existem pequenos ligados ao lado atual 
	if(highlist[is].Reference().Dimension() != InterfaceDimension) continue;
	TPZCompElDisc *del = dynamic_cast<TPZCompElDisc *> (highlist[is].Element());
	if(!del) continue;
	del->CreateInterface(highlist[is].Side());
      }
    }
  }
}

void TPZCompElDisc::CreateInterface(int side){

  TPZCompElSide thisside(this,side);
  TPZStack<TPZCompElSide> list;
  list.Resize(0);
  thisside.EqualLevelElementList(list,0,1);//retorna distinto ao atual ou nulo
  int size = list.NElements();
  //espera-se ter os elementos computacionais esquerdo e direito 
  //já criados antes de criar o elemento interface
  if(size){
    //neste caso existem ambos: esquerdo e direito, this e list[0], e são vizinhos
/*    int matid1 = fMaterial->Id();
    int matid2 = list[0].Element()->Material()->Id(); 
    int matid = (matid1 < matid2) ? matid1 : matid2; */
    
//Antes (comentado acima) a interface tinha o material de menor id entre esquerdo e direito
//Agora, a interface tem o material do elemento de menor dimensão (i.e. condicao de contorno)
//Seria interessante testar se o material do contorno eh filho do TPZBndCond, mas traria problemas ao cfdk.
    int matid;
    int thisdim = this->Dimension();
    int neighbourdim = list[0].Element()->Dimension();
    if (thisdim == neighbourdim) 
      matid = this->Material()->Id();
    else {
      if (thisdim < neighbourdim)
	matid = this->Material()->Id();
      else
	matid = list[0].Element()->Material()->Id();
    }


    TPZGeoEl *gel = fReference->CreateBCGeoEl(side,matid);
    //isto acertou as vizinhan¢as da interface geométrica com o atual
    int index;
    TPZCompElDisc *list0 = dynamic_cast<TPZCompElDisc *>(list[0].Element());
    if(Dimension() > list0->Dimension()){
      //o de volume é o direito caso um deles seja BC
      //a normal aponta para fora do contorno
      new TPZInterfaceElement(*fMesh,gel,index,this,list0/*,side*/);
    } else {
      //caso contrário ou caso ambos sejam de volume 
      new TPZInterfaceElement(*fMesh,gel,index,list0,this/*,list[0].Side()*/);
    }
    return;
  }
  //aqui não existe igual: só pode existir lower 
  //(pois isso foi verificado no CreateInterfaces())
  TPZCompElSide lower = thisside.LowerLevelElementList(0);
  if(lower.Exists()){
    int matid1 = fMaterial->Id();
    int matid2 = lower.Element()->Material()->Id();
    int matid = (matid1 < matid2) ? matid1 : matid2;
    //existem esquerdo e direito: this e lower
    TPZGeoEl *gel = fReference->CreateBCGeoEl(side,matid);
    int index;
    TPZCompElDisc *lowcel = dynamic_cast<TPZCompElDisc *>(lower.Element());
    if(Dimension() > lowcel->Dimension()){
      //para que o elemento esquerdo seja de volume
      new TPZInterfaceElement(*fMesh,gel,index,this,lowcel/*,side*/);
    } else {
      new TPZInterfaceElement(*fMesh,gel,index,lowcel,this/*,lower.Side()*/);
    }
    return;
  }
}

REAL TPZCompElDisc::NormalizeConst(){
  //maior distancia entre o ponto interior e os vértices do elemento
  int nnodes = fReference->NNodes(),i;
  if(nnodes == 1) return 1.0;//elemento ponto
  REAL maxdist,dist;
  int inode = fReference->NodeIndex(0);//primeiro nó do elemento
  TPZGeoNode node = fReference->Mesh()->NodeVec()[inode];
  maxdist = pow(node.Coord(0)-fCenterPoint[0],2.)+pow(node.Coord(1)-fCenterPoint[1],2.);
  maxdist += pow(node.Coord(2)-fCenterPoint[2],2.);
  maxdist = sqrt(maxdist);
  for(i=1;i<nnodes;i++){
    inode = fReference->NodeIndex(i);//nós sub-seguintes
    node = fReference->Mesh()->NodeVec()[inode];
    dist = pow(node.Coord(0)-fCenterPoint[0],2.)+pow(node.Coord(1)-fCenterPoint[1],2.);
    dist += pow(node.Coord(2)-fCenterPoint[2],2.);
    dist = sqrt(dist);
    if(maxdist < dist) maxdist = dist;
  }
  return maxdist;  
}

void TPZCompElDisc::Shape(TPZVec<REAL> X, TPZFMatrix &phi, TPZFMatrix &dphi) {

  if(fDegree < 0) return;
  if(Dimension()==0)
    TPZShapeDisc::Shape0D(fConstC,fCenterPoint,X,fDegree,phi,dphi);
  else
    if(Dimension()==1)
      TPZShapeDisc::Shape1D(fConstC,fCenterPoint,X,fDegree,phi,dphi);
    else
      if(Dimension()==2)
	TPZShapeDisc::Shape2D(fConstC,fCenterPoint,X,fDegree,phi,dphi,fShapefunctionType);
      else
	if(Dimension()==3)
	  TPZShapeDisc::Shape3D(fConstC,fCenterPoint,X,fDegree,phi,dphi,fShapefunctionType);  
}

void TPZCompElDisc::Print(ostream &out) {

  out << "\nDiscontinous element : \n";
  //out << "\tGeometric reference id : " << fReference->Id() << endl
  out << "\tMaterial id : " << fReference->MaterialId() << endl
      << "\tDegrau of interpolation : " <<  fDegree << endl
      << "\tConnect index : " << fConnectIndex << endl
      << "\tNormalizing constant : " << fConstC << endl
      << "\tCenter point of the element : ";
  int size = fCenterPoint.NElements(),i;
  for(i=0;i<size-1;i++) out << fCenterPoint[i] << " , ";
  out << fCenterPoint[i] << endl;
}

int TPZCompElDisc::ConnectIndex(int side) {

  return fConnectIndex;
}

int TPZCompElDisc::NConnects(){

  
  return (fConnectIndex !=-1);

}

int TPZCompElDisc::CreateMidSideConnect(){
  // primeiro são criados os elementos de volume depois os elementos BC associados aos seus lados
  // num estágio inicial o elemento BC é acoplado ao elemento ELV de volume de tal forma 
  // que ambos são vizinhos
  // o elemento BC não pode ser dividido se o elemento ELV associado não for dividido primeiro
  // caso o elemento ELV é dividido, então o elemento BC associado deveria ser dividido
  // tambem para manter a CC consistente com a malha
  // caso ELV é dividido e BC não é então ELV é LowerLevelElement do elemento BC
  if(!Material())
    PZError << "\nTPZCompElDisc::CreateMidSideConnect Material nulo\n";

  TPZStack<TPZCompElSide> list;
  int nsides = fReference->NSides();
  int dimgrid = fMaterial->Dimension();
  int dim = Dimension();
  int existsconnect = 0;

  if(dimgrid == dim){
    //este é um elemento de volume
    //procura-se elemento superposto
    TPZCompElSide(this,nsides-1).EqualLevelElementList(list,0,0);
    int size = list.NElements(),i;
    for(i=0;i<size;i++){
      int dimel = list[i].Element()->Reference()->Dimension();
      if(dimel == dimgrid){
	int connectindex = list[i].Element()->ConnectIndex(0);
	list[i].Element()->SetConnectIndex(0,connectindex);
	existsconnect = 1;
	break;
      }
    }   
  }
  
  if(dim == dimgrid - 1){ //dimgrid - 1 = interface dimension
    // o atual é um elemento BC
    fConnectIndex = -1;//=> return NshapeF() = 0
    fDegree = -1;
    return fConnectIndex;
  }

  if(!existsconnect){
    //o atual é um elemento de volume e
    //não achou-se um elemento superposto
    int nvar = Material()->NStateVariables();
    int newnodeindex = Mesh()->AllocateNewConnect();
    SetConnectIndex(0,newnodeindex);
    TPZConnect &newnod = Mesh()->ConnectVec()[newnodeindex];
    int seqnum = newnod.SequenceNumber();
    Mesh()->Block().Set(seqnum,nvar*NShapeF());
    Mesh()->ConnectVec()[fConnectIndex].IncrementElConnected();
  }

  return fConnectIndex;
}

int TPZCompElDisc::NShapeF(){


  if(fConnectIndex == -1) return 0;
  //deve ter pelo menos um connect

  int dim = Dimension();
  return TPZShapeDisc::NShapeF(fDegree,dim,fShapefunctionType);
}


void TPZCompElDisc::InternalPoint(TPZVec<REAL> &point){
  //ponto deformado
  point.Resize(3,0.);
  point[0] = fCenterPoint[0];
  point[1] = fCenterPoint[1];
  point[2] = fCenterPoint[2];
}


void TPZCompElDisc::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){

  if(fMaterial == NULL){
    cout << "TPZCompElDisc::CalcStiff : no material for this element\n";
    return;
  }
  int ncon = NConnects();
  int dim = Dimension();
  int nstate = fMaterial->NStateVariables();
  int nshape = NShapeF();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix &MeshSol = Mesh()->Solution();
  int numeq = nshape * nstate;

  // clean ek and ef
  if(!ek.fMat) ek.fMat = new TPZFMatrix();
  if(!ef.fMat) ef.fMat = new TPZFMatrix();
  if(!ek.fBlock) ek.fBlock = new TPZBlock(ek.fMat);
  if(!ef.fBlock) ef.fBlock = new TPZBlock(ef.fMat);

  ek.fMat->Redim(numeq,numeq);
  ef.fMat->Redim(numeq,1);
  if(ncon){//pode serr no máximo ncon = 1
    ek.fBlock->SetNBlocks(ncon);
    ef.fBlock->SetNBlocks(ncon); 
    ek.fBlock->Set(0,NShapeF()*nstate);
    ef.fBlock->Set(0,NShapeF()*nstate);
  }
  if( !ek.fMat || !ef.fMat || !ek.fBlock || !ef.fBlock){
    cout << "TPZInterpolatedElement.calc_stiff : not enough storage for local stifness"
      " matrix \n";
    Print(cout);
    if(ek.fMat)   delete ek.fMat;
    if(ek.fBlock) delete ek.fBlock;
    if(ef.fMat)   delete ef.fMat;
    if(ef.fBlock) delete ef.fBlock;
    ek.fMat=  NULL;
    ek.fBlock = NULL;
    ef.fMat = NULL;
    ef.fBlock = NULL;
    return;
  }
  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);
  for(int i=0;i<ncon;i++){
    (ef.fConnect)[i] = ConnectIndex(i);
    (ek.fConnect)[i] = ConnectIndex(i);
  }
  if(ncon==0) return;//elemento CC no passa
  TPZFMatrix phix(nshape,1),dphix(dim,nshape);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL detjac,weight;
  int integ = 2*Degree();
  TPZIntPoints *intrule = Reference()->CreateSideIntegrationRule(Reference()->NSides()-1,integ);
  int npoints = intrule->NPoints(),ip;                                              //integra fi*fj
  TPZVec<REAL> sol(nstate,0.);
  TPZFMatrix dsol(dim,nstate,0.);

  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    fReference->X(intpoint, x);
    weight *= fabs(detjac);
    Shape(x,phix,dphix);
    axes.Identity();
    //solu¢ão da itera¢ão anterior
    sol.Fill(0.);
    dsol.Zero();
    for(int in=0; in<ncon; in++) {
      TPZConnect *df = &Connect(in);
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      int pos = block.Position(dfseq);
      int iv = 0,d;
      for(int jn=0; jn<dfvar; jn++) {
	sol[iv%nstate] += phix(iv/nstate,0)*MeshSol(pos+jn,0);
	for(d=0; d<dim; d++)
	  dsol(d,iv%nstate) += dphix(d,iv/nstate)*MeshSol(pos+jn,0);
	iv++;
      }
    }
    fMaterial->Contribute(x,jacinv,sol,dsol,weight,axes,phix,dphix,*ek.fMat,*ef.fMat);
  }
}

void TPZCompElDisc::CalcResidual(TPZElementMatrix &ef){

  if(fMaterial == NULL){
    cout << "TPZCompElDisc::CalcStiff : no material for this element\n";
    return;
  }
  int ncon = NConnects();
  int dim = Dimension();
  int nstate = fMaterial->NStateVariables();
  int nshape = NShapeF();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix &MeshSol = Mesh()->Solution();
  int numeq = nshape * nstate;

  // clean ef
  if(!ef.fMat) ef.fMat = new TPZFMatrix();
  if(!ef.fBlock) ef.fBlock = new TPZBlock(ef.fMat);

  ef.fMat->Redim(numeq,1);
  if(ncon){//pode serr no máximo ncon = 1
    ef.fBlock->SetNBlocks(ncon); 
    ef.fBlock->Set(0,NShapeF()*nstate);
  }
  if( !ef.fMat || !ef.fBlock){
    cout << "TPZInterpolatedElement.calc_stiff : not enough storage for local stifness"
      " matrix \n";
    Print(cout);
    if(ef.fMat)   delete ef.fMat;
    if(ef.fBlock) delete ef.fBlock;
    ef.fMat = NULL;
    ef.fBlock = NULL;
    return;
  }
  ef.fConnect.Resize(ncon);
  for(int i=0;i<ncon;i++){
    (ef.fConnect)[i] = ConnectIndex(i);
  }
  if(ncon==0) return;//elemento CC no passa
  TPZFMatrix phix(nshape,1),dphix(dim,nshape);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL detjac,weight;
  int integ = 2*Degree();
  TPZIntPoints *intrule = Reference()->CreateSideIntegrationRule(Reference()->NSides()-1,integ);
  int npoints = intrule->NPoints(),ip;                                              //integra fi*fj
  TPZVec<REAL> sol(nstate,0.);
  TPZFMatrix dsol(dim,nstate,0.);

  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    fReference->X(intpoint, x);
    weight *= fabs(detjac);
    Shape(x,phix,dphix);
    axes.Identity();
    //solu¢ão da itera¢ão anterior
    sol.Fill(0.);
    dsol.Zero();
    for(int in=0; in<ncon; in++) {
      TPZConnect *df = &Connect(in);
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      int pos = block.Position(dfseq);
      int iv = 0,d;
      for(int jn=0; jn<dfvar; jn++) {
	sol[iv%nstate] += phix(iv/nstate,0)*MeshSol(pos+jn,0);
	for(d=0; d<dim; d++)
	  dsol(d,iv%nstate) += dphix(d,iv/nstate)*MeshSol(pos+jn,0);
	iv++;
      }
    }
    fMaterial->Contribute(x,jacinv,sol,dsol,weight,axes,phix,dphix,*ef.fMat);
  }
}

REAL TPZCompElDisc::SizeOfElement(){

  int dim = fReference->Dimension();
  int side = fReference->NSides()-1;
  if(dim == 2) fReference->SideArea(side);
  if(!dim || dim > 2){
    PZError << "TPZCompElDisc::SizeOfElement case not permited\n";
    return 0.;
  }
  if(dim == 1){
    TPZGeoNode node0 = Mesh()->Reference()->NodeVec()[fReference->NodeIndex(0)];
    TPZGeoNode node1 = Mesh()->Reference()->NodeVec()[fReference->NodeIndex(1)];
    TPZVec<REAL> no0(3),no1(3);
    for(int i=0;i<3;i++){
      no0[i] = node0.Coord(i);
      no1[i] = node1.Coord(i);
    }
    return fReference->Distance(no0,no1);
  }
  PZError << "TPZCompElDisc::SizeOfElement this in case that it is not contemplated\n";
  return 0.;
}

void TPZCompElDisc::Divide(int index,TPZVec<int> &subindex,int degree){

  if (fMesh->ElementVec()[index] != this) {
    PZError << "TPZInterpolatedElement::Divide index error";
    subindex.Resize(0);
    return;
  }
  if(Type() == EInterface){
    PZError << "TPZInterpolatedElement::Divide element interface cannot be divided!\n";
    exit(-1);
  }

  RemoveInterfaces();

  if(0){//TESTE
    ofstream mesh("MALHADIV.out");//TESTE
    Mesh()->Reference()->Print(mesh);//TESTE
    Mesh()->Print(mesh);//TESTE
    mesh.flush();  //TESTE
    mesh.close();//TESTE
    return;//TESTE
  }//TESTE

  //divide o elemento geométrico
  int nsubs = fReference->NSubElements();
  subindex.Resize(nsubs);
  TPZManVector<TPZGeoEl *> geosubs(nsubs);
  fReference->Divide(geosubs);
  if(!geosubs.NElements()) {
    subindex.Resize(0);
    return;
  }

  fReference->ResetReference();
  TPZCompElDisc *discel;
  int i,deg;
  if(degree) deg = degree;
  else deg = Degree();
  for (i=0;i<nsubs;i++){
    new TPZCompElDisc(*Mesh(),geosubs[i],subindex[i]);
    discel = (TPZCompElDisc *) fMesh->ElementVec()[subindex[i]];
    discel->SetDegree(deg);
  }

  Mesh()->ExpandSolution();
  for(i=0; i<nsubs; i++) {
    discel = (TPZCompElDisc *) fMesh->ElementVec()[subindex[i]];
    if(discel->Dimension() < fMaterial->Dimension()) continue;//elemento BC
    discel->InterpolateSolution(*this);
  }
  delete this;
}

void TPZCompElDisc::InterpolateSolution(TPZCompElDisc &coarsel){
  // accumulates the transfer coefficients between the current element and the
  // coarse element into the transfer matrix, using the transformation t
  TPZTransform t(Dimension());
  Reference()->BuildTransform2(fReference->NSides()-1,coarsel.Reference(),t);

  int locmatsize = NShapeF();
  int cormatsize = coarsel.NShapeF();
  int nvar = fMaterial->NStateVariables();
  int dimension = Dimension();

  TPZFMatrix loclocmat(locmatsize,locmatsize,0.);
  TPZFMatrix projectmat(locmatsize,nvar,0.);

  TPZVec<int> prevorder(dimension),order(dimension);
  TPZIntPoints *intrule = Reference()->CreateSideIntegrationRule(Reference()->NSides()-1,Degree());
  intrule->GetOrder(prevorder);
  int i;
  for(i=0;i<dimension;i++) order[i] = Degree()*2;
  intrule->SetOrder(order);

  TPZFMatrix locphi(locmatsize,1);
  TPZFMatrix locdphi(dimension,locmatsize);	// derivative of the shape function
  // in the master domain

  TPZFMatrix corphi(cormatsize,dimension);
  TPZFMatrix cordphi(dimension,cormatsize);	// derivative of the shape function
  // in the master domain

  TPZVec<REAL> int_point(dimension),coarse_int_point(dimension);
  TPZFMatrix jacobian(dimension,dimension),jacinv(dimension,dimension);
  TPZFMatrix axes(3,3,0.);
  REAL zero = 0.;
  TPZVec<REAL> x(3,zero);
  TPZVec<REAL> u(nvar);

  int numintpoints = intrule->NPoints();
  REAL weight;
  int lin,ljn,cjn;
  TPZConnect *df;
  TPZBlock &block = Mesh()->Block();

  for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {

    intrule->Point(int_ind,int_point,weight);
    REAL jac_det = 1.;
    Reference()->Jacobian( int_point, jacobian , axes,jac_det,jacinv);
    Reference()->X(int_point, x);
    Shape(x,locphi,locdphi);
    axes.Identity();
#warning "Este codigo esta errado!!"

    weight *= jac_det;
    t.Apply(int_point,coarse_int_point);
    coarsel.Shape(coarse_int_point,corphi,cordphi);
    u.Fill(0.);
    int iv = 0;
    
    df = &coarsel.Connect(0);
    int dfseq = df->SequenceNumber();
    
    int dfvar = block.Size(dfseq);
    for(ljn=0; ljn<dfvar; ljn++) {
      u[iv%nvar] += corphi(iv/nvar,0)*block(dfseq,0,ljn,0);
      iv++;
    }

    for(lin=0; lin<locmatsize; lin++) {
      for(ljn=0; ljn<locmatsize; ljn++) {
	loclocmat(lin,ljn) += weight*locphi(lin,0)*locphi(ljn,0);
      }
      for(cjn=0; cjn<nvar; cjn++) {
	projectmat(lin,cjn) += weight*locphi(lin,0)*u[cjn];
      }
    }
    jacobian.Zero();
  }
  loclocmat.SolveDirect(projectmat,ELU);
  // identify the non-zero blocks for each row
  int iv=0;

  df = &Connect(0);
  int dfseq = df->SequenceNumber();
  int dfvar = block.Size(dfseq);
  for(ljn=0; ljn<dfvar; ljn++) {
    block(dfseq,0,ljn,0) = projectmat(iv/nvar,iv%nvar);
    iv++;
  }

  intrule->SetOrder(prevorder);
  
}

int TPZCompElDisc::ExistsInterface(TPZGeoElSide geosd){

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

void TPZCompElDisc::RemoveInterfaces(){

  int nsides = Reference()->NSides();
  int InterfaceDimension = fMaterial->Dimension();
  int is;
  TPZStack<TPZCompElSide> list,equal;
  for(is=0;is<nsides;is++){
    TPZCompElSide thisside(this,is);
    if(thisside.Reference().Dimension() != InterfaceDimension) continue;
    // procurar na lista de elementos iguais
    list.Resize(0);// o lado atual é uma face
    //thisside.EqualLevelElementList(list,0,0);// monta a lista de elementos iguais
    RemoveInterface(is);// chame remove interface do elemento atual (para o side atual)
    thisside.HigherLevelElementList(list,0,0);// procurar na lista de elementos menores (todos)
    int size = list.NElements(),i;            // 'isto pode incluir elementos interfaces'
    //tirando os elementos de interface da lista
    for(i=0;i<size;i++){
      if(list[i].Element()->Type() == EInterface) list[i] = TPZCompElSide();//tirando interface
    }
    for(i=0;i<size;i++){// percorre os elementos menores
      if(!list[i].Element()) continue;
      TPZGeoElSide geolist = list[i].Reference();//TESTE
      if(geolist.Dimension() != InterfaceDimension) continue;
      equal.Resize(0);// para cada elemento menor e' preciso verificar a dimensao,
      list[i].EqualLevelElementList(equal,0,0);//montar a lista de elementos iguais (todos)
      equal.Push(list[i]);//não é incorporado no método anterior
      int neq = equal.NElements(),k=-1;
      while(++k < neq) if(equal[k].Element()->Type() != EInterface) break;//procurando elemento descontínuo cujo
      if(!neq || k == neq){                               //lado faz parte da parti¢ão do lado side do this
	PZError << "TPZCompElDisc::RemoveInterfaces inconsistency of data";
	exit(-1);//elemento descontínuo não achado: ERRO
      }// chame removeinterface do elemento descontinuo menor
      ((TPZCompElDisc *)equal[k].Element())->RemoveInterface(equal[k].Side());
    }
  }

}

void TPZCompElDisc::RemoveInterface(int side) {
  
  TPZStack<TPZCompElSide> list;
  list.Resize(0);
  TPZCompElSide thisside(this,side);
  thisside.EqualLevelElementList(list,0,0);// monta a lista de elementos iguais
  int size = list.NElements(),i=-1;
  while(++i < size) if(list[i].Element()->Type() == EInterface) break;// procura aquele que e derivado de TPZInterfaceEl
  if(!size || i == size){
    //PZError << "\nTPZCompElDisc::RemoveInterface interface element not found (no problems)\n";
    return;// nada a ser feito
  }
  // aqui existe a interface
  TPZCompEl *cel = list[i].Element();
  TPZGeoEl *gel = cel->Reference();
  gel->RemoveConnectivities();// deleta o elemento das vizinhancas
  TPZGeoMesh *gmesh = Mesh()->Reference();
  int index = gmesh->ElementIndex(gel);// identifica o index do elemento
  gmesh->ElementVec()[index] = NULL;
  delete cel;
  delete gel;// deleta o elemento
  gmesh->ElementVec().SetFree(index);// Chame SetFree do vetor de elementos da malha geometrica para o index
}


void TPZCompElDisc::Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol) {
  //#ifdef _AUTODIFF
  //  TPZConservationLaw2 *mat = dynamic_cast<TPZConservationLaw2 *>(fMaterial);
  //#else
  //  TPZConservationLaw *mat = dynamic_cast<TPZConservationLaw *>(fMaterial);
  //#endif

  if(var >= 100) {
    TPZCompEl::Solution(qsi,var,sol);
    return;
  }
  int nshape = NShapeF();
  int dim = Dimension();
  int ncon = NConnects();
  if(var == 99) {
    sol[0] = Degree();
    return;
  }
  TPZBlock &block = fMesh->Block();
  TPZFMatrix &Sol = fMesh->Solution();

  if(fMaterial == NULL){
    PZError << "TPZIntEl::Solution : no Material for this element\n";
    Print(PZError);
    return;
  }

  int numdof = fMaterial->NStateVariables();
  TPZFNMatrix<220> phi(nshape,1);
  TPZFNMatrix<660> dphi(dim,nshape);
  TPZManVector<REAL> u(numdof);
  TPZFMatrix du(dim,numdof,0.);
  TPZFMatrix axes(3,3,0.);
  REAL jacstore[10],jacinvstore[10];
  TPZFMatrix jacobian(dim,dim,jacstore,10);
  TPZFMatrix jacinv(dim,dim,jacinvstore,10);
  TPZManVector<REAL> x(3);
  //REAL detjac;
  //fReference->Jacobian(qsi,jacobian,axes,detjac,jacinv);//(calcula axes)
  if(var >= 0){
    fReference->X(qsi,x);
  } else if(var < 0){
    //neste caso 0 ponto qsi está no elemento deformado
    var *= -1;//recuperando var
    for(int i=0;i<3;i++) x[i] = qsi[i];
  }
  Shape(x,phi,dphi);
  int iv=0,in,jn,d;
  TPZConnect *df;
  u.Fill(0.);
  du.Zero();
  for(in=0; in<ncon; in++) {
    df = &Connect(in);
    int dfseq = df->SequenceNumber();
    int dfvar = block.Size(dfseq);
    int pos = block.Position(dfseq);
    for(jn=0; jn<dfvar; jn++) {
      u[iv%numdof] += phi(iv/numdof,0)*Sol(pos+jn,0);
      for(d=0; d<dim; d++){
	du(d,iv%numdof) += dphi(d,iv/numdof)*Sol(pos+jn,0);
      }
      iv++;
    }
  }
  fMaterial->Solution(u,du,axes,var,sol);
}

void TPZCompElDisc::Solution(TPZVec<REAL> &x, TPZVec<REAL> &uh){

  TPZCompMesh *finemesh = Mesh();
  TPZBlock &fineblock = finemesh->Block();
  int nstate = Material()->NStateVariables();
  TPZFMatrix &FineMeshSol = finemesh->Solution();
  int matsize = NShapeF(),dim = Dimension();
  TPZFMatrix phix(matsize,1,0.);
  TPZFMatrix dphix(dim,matsize,0.);
  Shape(x,phix,dphix);
  TPZConnect *df = &Connect(0);
  int dfseq = df->SequenceNumber();
  int dfvar = fineblock.Size(dfseq);
  int pos   = fineblock.Position(dfseq);
  int iv = 0,d;
  uh.Fill(0.);
  for(d=0; d<dfvar; d++) {
    uh[iv%nstate] += phix(iv/nstate,0)*FineMeshSol(pos+d,0);
    iv++;
  }
}



void TPZCompElDisc::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
  int mat = Material()->Id();
  int nsides = fReference->NSides();

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

int TPZCompElDisc::NSides(){

  return fReference->NSides();
}

int TPZCompElDisc::NInterfaces(){

  int nsides = this->NSides();

  switch( nsides )
    {
    case 3: //line
      return 2;
      break;

    case 7: //triangle
      return 3;
      break;

    case 9: //square
      return 4;      
      
    case 15: // Tetrahedra.
      return 4;
      break;
	
    case 19: // Prism.
      return 5;
      break;	

    case 21: // Pyramid.
      return 6;
      break;

    case 27: // Hexaedra.
      return 8;
      break;

    default:
      PZError << "TPZCompElDisc::NFaces() - Unknown element shape!" << endl;
      exit (-1);
    }
}

//#include "TPZAgglomerateEl.h"
void TPZCompElDisc::AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight){

  int i,npoints;
  TPZVec<REAL> pt(3),x(3,0.0);
  TPZFMatrix jacobian(3,3),jacinv(3,3),axes(3,3);
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
}


void TPZCompElDisc::CenterPoint(TPZVec<REAL> &center){

  if(Reference() || Type() == EDiscontinuous){
    fReference->CenterPoint(fReference->NSides()-1,center);
    return;
  } else {//aglomerado
//     TPZStack<TPZCompEl *> elvec;
//     dynamic_cast<TPZAgglomerateElement *>(this)->ListOfDiscEl(elvec);
//     TPZGeoEl *ref = elvec[0]->Reference();
//     ref->CenterPoint(ref->NSides()-1,center);
    PZError << "TPZCompElDisc::CenterPoint center points not exists!\n";
  }
}



void TPZCompElDisc::EvaluateError(  void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
				    TPZVec<REAL> &errors,TPZBlock * /*flux */) {

  int NErrors = this->Material()->NEvalErrors();
  errors.Resize(NErrors);
  errors.Fill(0.);
  if(fMaterial == NULL){
    PZError << "TPZInterpolatedElement::EvaluateError : no material for this element\n";
    Print(PZError);
    return;
  }
  int problemdimension = Mesh()->Dimension();
  if(Reference()->Dimension() < problemdimension) return;
  // Adjust the order of the integration rule
  int nsides = Reference()->NSides();
  TPZIntPoints *intrule = Reference()->CreateSideIntegrationRule(nsides-1,20);
  int dim = Dimension();

  int ndof = fMaterial->NStateVariables();
  int nflux = fMaterial->NFluxes();
  int nshape = NShapeF();
  //suficiente para ordem 5 do cubo
  TPZFNMatrix<220> phi(nshape,1);
  TPZFNMatrix<660> dphi(dim,nshape),dphix(dim,nshape);
  TPZFNMatrix<9> jacobian(dim,dim),jacinv(dim,dim);
  TPZFNMatrix<9> axes(3,3);
  TPZManVector<REAL,3> x(3);//TPZVec<REAL> x(3,0.);
  TPZManVector<REAL,6> u_exact(ndof);
  TPZFNMatrix<90> du_exact(dim,ndof);
  TPZManVector<REAL,3> intpoint(3), values(NErrors);
  values.Fill(0.0);
  REAL detjac,weight;
  TPZManVector<REAL,6> u(ndof);
  TPZFNMatrix<90> dudx(dim,ndof);
  TPZManVector<REAL,9> flux_el(nflux,0.);
  TPZMaterial *matp = fMaterial;
  int ncon = NConnects();
  TPZBlock &block = Mesh()->Block();

  for(int nint=0; nint<intrule->NPoints(); nint++) {

    intrule->Point(nint,intpoint,weight);
    fReference->Jacobian( intpoint , jacobian, axes, detjac , jacinv);
    fReference->X( intpoint , x);
    Shape(x,phi,dphix);
    axes.Identity();
    weight *= fabs(detjac);
    int iv=0,in,jn,d;
    TPZConnect *df;
    u.Fill(0.);
    dudx.Zero();
    for(in=0; in<ncon; in++) {
      df = &Connect(in);
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      for(jn=0; jn<dfvar; jn++) {
	u[iv%ndof] += phi(iv/ndof,0)*block(dfseq,0,jn,0);
	for(d=0; d<dim; d++)
	  dudx(d,iv%ndof) += dphix(d,iv/ndof)*block(dfseq,0,jn,0);
	iv++;
      }
    }//solucao calculculada no sistema local : elementos 2d
    //contribuções dos erros
    if(fp) {
      fp(x,u_exact,du_exact);
      matp->Errors(x,u,dudx,axes,flux_el,u_exact,du_exact,values);
      for(int ier = 0; ier < NErrors; ier++)
	errors[ier] += values[ier]*weight;
    }
  }//fim for : integration rule
   //Norma sobre o elemento
  for(int ier = 0; ier < NErrors; ier++)
    errors[ier] = sqrt(errors[ier]);
  delete intrule;
}


void TPZCompElDisc::BuildTransferMatrix(TPZCompElDisc &coarsel, TPZTransfer &transfer){
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

  TPZFNMatrix<500> loclocmat(locnshape,locnshape);
  TPZFMatrix loccormat(locnshape,cornshape);
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

  TPZFMatrix corphi(cornshape,1);
  TPZFMatrix cordphi(dimension,cornshape);
  // derivative of the shape function
  // in the master domain

  TPZManVector<REAL> int_point(dimension);
  TPZFNMatrix<9> jacobian(dimension,dimension);
  TPZFMatrix jacinv(dimension,dimension);
  TPZFNMatrix<9> axes(3,3);
  TPZManVector<REAL> x(3);

  int_point.Fill(0.,0);
  REAL jac_det = 1.;
  fReference->Jacobian( int_point, jacobian , axes, jac_det, jacinv);
  REAL multiplier = 1./jac_det;

  int numintpoints = intrule->NPoints();
  REAL weight;
  int lin,ljn,cjn;

  for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {

    intrule->Point(int_ind,int_point,weight);
    fReference->Jacobian( int_point, jacobian , axes, jac_det, jacinv);
    fReference->X(int_point, x);
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


  int locblockseq = Connect(0).SequenceNumber();
  TPZStack<int> globblockvec;
  int numnonzero = 0;
  int cind = coarsel.ConnectIndex(0);
  TPZConnect &con = coarsel.Mesh()->ConnectVec()[cind];
  int corblockseq = con.SequenceNumber();
  if(locnshape == 0 || cornshape == 0)
    PZError << "TPZCompElDisc::BuilTransferMatrix error I\n";
  TPZFMatrix small(locnshape,cornshape,0.);
  loccormat.GetSub(0,0,locnshape,cornshape,small);
  REAL tol = Norm(small);
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
  loccormat.GetSub(0,0,locnshape,cornshape,small);
  transfer.SetBlockMatrix(locblockseq,globblockvec[0],small);

  SetDegree(locdeg);
}

void TPZCompElDisc::AccumulateVertices(TPZStack<TPZGeoNode *> &nodes) {
  TPZGeoEl *geo = Reference();
#warning "Este metodo nao funciona para aglomerados contendo aglomerados"
  if(!geo) {
    PZError <<  "TPZCompElDisc::AccumulateVertices null reference\n";
    return;
  }
  int nvertices = geo->NNodes();
  int l;
  for(l=0;l<nvertices;l++) nodes.Push( geo->NodePtr(l) );
}
