// _*_ c++ _*_ 
#include "pzelgc3d.h"
#include "pzelgt3d.h"
#include "pzelgpi3d.h"
#include "pzelgpr3d.h"
#include "pzelgt2d.h"
#include "pzelgq2d.h"
//#include "pzelcq2d.h"
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
//#include "TPZRefPattern.h"
#include "time.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include <math.h>
#include <stdio.h>

int TPZCompElDisc::gInterfaceDimension = 2;


TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,TPZGeoEl *ref,int &index) :
		TPZCompEl(mesh,index), fCenterPoint(3) {//dois pontos por eixo
  int i;
  fDegree = gOrder;
  fReference = ref;
  ref->SetReference(this);
  fMesh = &mesh;
  int materialid = ref->MaterialId();
  fMaterial = mesh.FindMaterial(materialid);
  fConnectIndex  = CreateMidSideConnect();
  mesh.ConnectVec()[fConnectIndex].IncrementElConnected();
  ref->CenterPoint(ref->NSides()-1,fCenterPoint);
  TPZVec<REAL> csi(fCenterPoint);
  ref->X(csi,fCenterPoint);
  fConstC = NormalizeConst();
  //criando os elementos interface
  CreateInterfaces();//argumento: número de faces
}

//INTERFACE
void TPZCompElDisc::CreateInterfaces(){
  
  int nsides = fReference->NSides();
  int side;
  nsides--;//last face
  for(side=nsides;side>=0;side--){
    if(fReference->SideDimension(side) != gInterfaceDimension) continue;
    TPZCompElSide thisside(this,side);
    if(ExistsInterface(thisside)) {
      cout << "TPZCompElDisc::CreateInterface inconsistent\n";
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
    TPZGeoEl *gel = fReference->CreateBCGeoEl(side,fMaterial->Id());
    //isto acertou as vizinhanas da interface geométrica com o atual
    int index;
    new TPZInterfaceElement(*fMesh,gel,index,*this);
    return;
  }
  //aqui não existe igual: só pode existir lower 
  //(pois isto foi verificado no CreateInterfaces())
  TPZCompElSide lower = thisside.LowerLevelElementList(0);
  if(lower.Exists()){
    //existem esquerdo e direito: this e lower
    TPZGeoEl *gel = fReference->CreateBCGeoEl(side,fMaterial->Id());
    int index;
    new TPZInterfaceElement(*fMesh,gel,index,*lower.Element());
    return;
  }
}

REAL TPZCompElDisc::NormalizeConst(){
  //maior distancia entre o ponto interior e os vértices do elemento
  int nnodes = fReference->NNodes(),i;
  REAL maxdist,dist;
  TPZGeoNode node = fReference->Mesh()->NodeVec()[0];
  maxdist = pow(node.Coord(0)-fCenterPoint[0],2.)+pow(node.Coord(1)-fCenterPoint[1],2.);
  maxdist += pow(node.Coord(2)-fCenterPoint[2],2.);
  maxdist = sqrt(maxdist);
  for(i=1;i<nnodes;i++){
    node = fReference->Mesh()->NodeVec()[i];
    dist = pow(node.Coord(0)-fCenterPoint[0],2.)+pow(node.Coord(1)-fCenterPoint[1],2.);
    dist += pow(node.Coord(2)-fCenterPoint[2],2.);
    dist = sqrt(maxdist);
    if(maxdist < dist) maxdist = dist;
  }
  return maxdist;  
}

void TPZCompElDisc::Shape(TPZVec<REAL> X, TPZFMatrix &phi, TPZFMatrix &dphi) {

  if(Dimension()==1)
    TPZShapeDisc::Shape1D(fConstC,fCenterPoint,X,fDegree,phi,dphi);
  if(Dimension()==2)
    TPZShapeDisc::Shape2D(fConstC,fCenterPoint,X,fDegree,phi,dphi);
  if(Dimension()==3)
    TPZShapeDisc::Shape3D(fConstC,fCenterPoint,X,fDegree,phi,dphi);  
}

void TPZCompElDisc::Print(ostream &out) {

  out << "\nDiscontinous element : \n";
  out << "\tGeometric reference id : " << fReference->Id() << endl
      << "\tMaterial id : " << fReference->MaterialId() << endl
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


int TPZCompElDisc::CreateMidSideConnect(){

  if(!Material())
    PZError << "\nTPZCompElDisc::CreateMidSideConnect Material nulo\n";
  int newnodeindex = Mesh()->AllocateNewConnect();
  return newnodeindex;
}

int TPZCompElDisc::NShapeF(){

  int dim = Dimension();
  int i,sum=0;
  switch(dim){  
    case 0:
      return 0;
    case 1:
      return (fDegree+1);
    case 2:
      return (fDegree+1)*(fDegree+2)/2;
    case 3:
      for(i=0;i<fDegree;i++) sum += (i+1)*(i+2)/2;;
      return sum;
    default:
      PZError << "TPZCompElDisc::NShapeF case not exists\n";
  }
}


void TPZCompElDisc::InternalPoint(TPZVec<REAL> &point){

  point.Resize(3,0.);
  point[0] = fCenterPoint[0];
  point[1] = fCenterPoint[1];
  point[2] = fCenterPoint[2];
}


void TPZCompElDisc::CalcStiffDisc(TPZFMatrix &ek, TPZFMatrix &ef){

  if(fMaterial == NULL){
    cout << "TPZCompElDisc::CalcStiff : no material for this element\n";
    return;
  }

  int dim = Dimension();
  int nstate = fMaterial->NStateVariables();
  int nshape = NShapeF();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix &MeshSol = Mesh()->Solution();
  int numeq = nshape * nstate;
  ek.Redim(numeq,numeq);
  ef.Redim(numeq,1);

  TPZFMatrix phix(nshape,1),dphix(dim,nshape);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL detjac,weight,C = ConstC();
  TPZIntPoints *intrule = Reference()->CreateSideIntegrationRule(Reference()->NSides()-1,Degree());
  int npoints = intrule->NPoints(),ip;
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
    TPZConnect *df = &Connect(0);
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
    fMaterial->Contribute(x,jacinv,sol,dsol,weight,axes,phix,dphix,ek,ef);
  }
}

void TPZCompElDisc::Divide(int index,TPZVec<int> &subindex,int degree){

  if (fMesh->ElementVec()[index] != this) {
    PZError << "TPZInterpolatedElement::Divide index error";
    subindex.Resize(0);
    return;
  }
  if(Type() == 17){
    PZError << "TPZInterpolatedElement::Divide element interface cannot be divided!\n";
    exit(-1);
  }

  RemoveInterfaces();

  if(0){//TESTE
    ofstream mesh("MALHADIV0.out");//TESTE
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

  cout << "\nTPZCompElDisc::Divide INTERPOLATE SOLUTION COM BUG\n";
  if(0){//COM PROBLEMAS
    Mesh()->ExpandSolution();
    for(i=0; i<nsubs; i++) {
      discel = (TPZCompElDisc *) fMesh->ElementVec()[subindex[i]];
      discel->InterpolateSolution(*this);
    }
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

int TPZCompElDisc::ExistsInterface(TPZCompElSide compsd){

  TPZGeoElSide geosd = compsd.Reference();
  TPZGeoElSide  neighside = geosd.Neighbour();
  while(neighside.Element() && neighside.Element() != geosd.Element()){
    TPZCompElSide neighcompside = neighside.Reference();
    neighside = neighside.Neighbour();
    if(!neighcompside.Element()) continue;
    if(neighcompside.Element()->Type() == 17) 
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
      if(list[i].Element()->Type() == 17) list[i] = TPZCompElSide();//tirando interface
    }
    for(i=0;i<size;i++){// percorre os elementos menores
      if(!list[i].Element()) continue;
      TPZGeoElSide geolist = list[i].Reference();//TESTE
      if(geolist.Dimension() != gInterfaceDimension) continue;
      equal.Resize(0);// para cada elemento menor e' preciso verificar a dimensao,
      list[i].EqualLevelElementList(equal,0,0);//montar a lista de elementos iguais (todos)
      equal.Push(list[i]);//não é incorporado no método anterior
      int neq = equal.NElements(),k=-1;
      while(++k < neq) if(equal[k].Element()->Type() != 17) break;//procurando elemento descontínuo cujo
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
  while(++i < size) if(list[i].Element()->Type() == 17) break;// procura aquele que e derivado de TPZInterfaceEl
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


/*
            //não apagar - NÃO APAGAR

    for(i=0;i<size;i++){
      if(list[i].Element()->Type() == 17){
	if(0){
	  TPZGeoElSide geoside = list[i].Reference();
	  TPZGeoElSide  neigh = geoside.Neighbour();
	  TPZCompElSide comp;
	  while(neigh.Element() && neigh.Element() != geoside.Element()){
	    comp = neigh.Reference();
	    neigh = neigh.Neighbour();
	    if(!comp.Element()) continue;
	    if(comp.Element()->Type() != 17) break;
	  }
	  if(neigh.Element() == geoside.Element()) list[i] = TPZCompElSide();
	  else list[i] = comp;//substitui-se o elem. interface por elem. de volume
	} else {
	  list[i] = TPZCompElSide();//tirando o elemento interface
	}
      }
    }
*/
