// _*_ c++ _*_ 

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
#include "pzgraphel.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "pzgraphel1d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pztrigraph.h"
#include "pzgraphel.h"
#include "time.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include <math.h>
#include <stdio.h>

int TPZCompElDisc::gDegree = 0;
int TPZCompElDisc::gInterfaceDimension = 2;//default

//construtor do elemento aglomerado
TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,int &index) :
		TPZCompEl(mesh,index), fCenterPoint(3) {
  fDegree = gDegree;
  fReference = NULL;
  fMaterial = NULL;
}

//construtor do elemento descontínuo
TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,TPZGeoEl *ref,int &index) :
		TPZCompEl(mesh,index), fCenterPoint(3) {
  fDegree = gDegree;
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
  
  int nsides = fReference->NSides();
  int side;
  nsides--;//last face
  for(side=nsides;side>=0;side--){
    if(fReference->SideDimension(side) != gInterfaceDimension) continue;
    TPZCompElSide thisside(this,side);
    if(ExistsInterface(thisside.Reference())) {
      int stop;
      cout << "TPZCompElDisc::CreateInterface inconsistent: interface already exists\n";
      cin >> stop;
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
	if(highlist[is].Reference().Dimension() != gInterfaceDimension) continue;
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
    int matid = fMaterial->Id();
    if(matid < 0){
      matid = list[0].Element()->Material()->Id();
    }
    TPZGeoEl *gel = fReference->CreateBCGeoEl(side,matid);
    //isto acertou as vizinhanas da interface geométrica com o atual
    int index;
    TPZCompElDisc *list0 = dynamic_cast<TPZCompElDisc *>(list[0].Element());
    if(Dimension() > list0->Dimension()){
      //o de volume é o direito caso um deles seja BC
      //a normal aponta para fora do contorno
      new TPZInterfaceElement(*fMesh,gel,index,this,list0,side);
    } else {
      //caso contrário ou caso ambos sejam de volume 
      new TPZInterfaceElement(*fMesh,gel,index,list0,this,list[0].Side());
    }
    return;
  }
  //aqui não existe igual: só pode existir lower 
  //(pois isso foi verificado no CreateInterfaces())
  TPZCompElSide lower = thisside.LowerLevelElementList(0);
  if(lower.Exists()){
    //existem esquerdo e direito: this e lower
    TPZGeoEl *gel = fReference->CreateBCGeoEl(side,fMaterial->Id());
    int index;
    TPZCompElDisc *lowcel = dynamic_cast<TPZCompElDisc *>(lower.Element());
    if(Dimension() > lowcel->Dimension()){
      //para que o elemento esquerdo seja de volume
      new TPZInterfaceElement(*fMesh,gel,index,this,lowcel,side);
    } else {
      new TPZInterfaceElement(*fMesh,gel,index,lowcel,this,lower.Side());
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
	TPZShapeDisc::Shape2D(fConstC,fCenterPoint,X,fDegree,phi,dphi);
      else
	if(Dimension()==3)
	  TPZShapeDisc::Shape3D(fConstC,fCenterPoint,X,fDegree,phi,dphi);  
}

void TPZCompElDisc::Print(ostream &out) {

  if(Type() == 15){//descontínuo
    out << "\nDiscontinous element : \n";
    out << "\tGeometric reference id : " << fReference->Id() << endl;
  }
  out << "\tMaterial id : " << Material()->Id() << endl
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
  int dimgrid = 1+gInterfaceDimension;
  int dim = fReference->Dimension();
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
  
  if(dim == gInterfaceDimension){
    // o atual é um elemento BC
    fConnectIndex = -1;
    fDegree = -1;//=> nshape = 0
    return fConnectIndex;
  }

  if(!existsconnect){
    //o atual é um elemento de volume e
    //não achou-se um elemento superposto
    int nvar = Material()->NStateVariables();
    int newnodeindex = Mesh()->AllocateNewConnect();
    TPZConnect &newnod = Mesh()->ConnectVec()[newnodeindex];
    int seqnum = newnod.SequenceNumber();
    Mesh()->Block().Set(seqnum,nvar*NShapeF());
    SetConnectIndex(0,newnodeindex);
    Mesh()->ConnectVec()[fConnectIndex].IncrementElConnected();
  }

  return fConnectIndex;
}

int TPZCompElDisc::NShapeF(){

  if(fConnectIndex == -1) return 0;
  //deve ter pelo menos um connect
  int dim = Dimension();
  int i,sum=0;
  switch(dim){  
    case 0:
      //return 1;
    case 1:
      return (fDegree+1);
    case 2:
      return (fDegree+1)*(fDegree+2)/2;
    case 3:
      for(i=0;i<(fDegree+1);i++) sum += (i+1)*(i+2)/2;;
      return sum;
    default:
      PZError << "TPZCompElDisc::NShapeF case not exists\n";
      return -1;
  }
}


void TPZCompElDisc::InternalPoint(TPZVec<REAL> &point){

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
  if(Type() == 16){
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
    Shape(int_point,locphi,locdphi);
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
    if(neighcompside.Element()->Type() == 16) 
      return 1;
  }
  return 0;
}

void TPZCompElDisc::RemoveInterfaces(){

  int nsides = Reference()->NSides();
  int is;
  TPZStack<TPZCompElSide> list,equal;
  for(is=0;is<nsides;is++){
    TPZCompElSide thisside(this,is);
    if(thisside.Reference().Dimension() != gInterfaceDimension) continue;
    // procurar na lista de elementos iguais
    list.Resize(0);// o lado atual é uma face
    //thisside.EqualLevelElementList(list,0,0);// monta a lista de elementos iguais
    RemoveInterface(is);// chame remove interface do elemento atual (para o side atual)
    thisside.HigherLevelElementList(list,0,0);// procurar na lista de elementos menores (todos)
    int size = list.NElements(),i;            // 'isto pode incluir elementos interfaces'
    //tirando os elementos de interface da lista
    for(i=0;i<size;i++){
      if(list[i].Element()->Type() == 16) list[i] = TPZCompElSide();//tirando interface
    }
    for(i=0;i<size;i++){// percorre os elementos menores
      if(!list[i].Element()) continue;
      TPZGeoElSide geolist = list[i].Reference();//TESTE
      if(geolist.Dimension() != gInterfaceDimension) continue;
      equal.Resize(0);// para cada elemento menor e' preciso verificar a dimensao,
      list[i].EqualLevelElementList(equal,0,0);//montar a lista de elementos iguais (todos)
      equal.Push(list[i]);//não é incorporado no método anterior
      int neq = equal.NElements(),k=-1;
      while(++k < neq) if(equal[k].Element()->Type() != 16) break;//procurando elemento descontínuo cujo
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
  while(++i < size) if(list[i].Element()->Type() == 16) break;// procura aquele que e derivado de TPZInterfaceEl
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

  TPZConservationLaw *mat = dynamic_cast<TPZConservationLaw *>(fMaterial);
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
  REAL phistore[220],dphistore[660],dphixstore[660];
  TPZFMatrix phi(nshape,1,phistore,220);
  TPZFMatrix dphi(dim,nshape,dphistore,660),dphix(dim,nshape,dphixstore,660);
  TPZManVector<REAL> u(numdof);
  TPZFMatrix du(dim,numdof,0.);
  TPZFMatrix axes(3,3,0.);
  REAL jacstore[10],jacinvstore[10];
  TPZFMatrix jacobian(dim,dim,jacstore,10);
  TPZFMatrix jacinv(dim,dim,jacinvstore,10);
  TPZManVector<REAL> x(3);
  REAL detjac;
  fReference->Jacobian(qsi,jacobian,axes,detjac,jacinv);//(calcula axes)
  fReference->X(qsi,x);
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
	du(d,iv%numdof) += dphix(d,iv/numdof)*Sol(pos+jn,0);
      }
      iv++;
    }
  }
  mat->Solution(u,du,axes,var,sol);
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

void TPZCompElDisc::AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight){

  int el,i,npoints;
  TPZVec<REAL> pt(3),x(3,0.0);
  TPZFMatrix jacobian(3,3),jacinv(3,3),axes(3,3);
  REAL detjac,wt;

  TPZGeoEl *subgel = Reference();
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

