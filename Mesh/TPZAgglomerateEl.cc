//$Id: TPZAgglomerateEl.cc,v 1.13 2003-11-07 12:57:52 cedric Exp $

#include "TPZAgglomerateEl.h"
#include "TPZInterfaceEl.h"
#include "TPZEulerConsLaw.h"
#include "TPZConservationLaw.h"
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

TPZAgglomerateElement::TPZAgglomerateElement(int nummat,int &index,TPZCompMesh &aggcmesh,TPZCompMesh *finemesh) : 
  TPZCompElDisc(aggcmesh,index) {

  /** 
   * o algomerado aponta para nulo mais o elemento computacional 
   * que ele agrupa aponta para o geométrico original
   * a copia do material na nova malha (malha aglomerada) neste
   * ponto já devia existir
   */
  fMotherMesh = finemesh;
  TPZMaterial *mater = Mesh()->FindMaterial(nummat);
  if(mater){
    SetMaterial(mater);
  } else {
    //cria copia dos materiais da malha fina
    CreateMaterialCopy(aggcmesh);
    mater = Mesh()->FindMaterial(nummat);
    SetMaterial(mater);
    if(!mater) PZError << "TPZAgglomerateElement::TPZAgglomerateElement material not found\n";
  }
  CreateMidSideConnect();
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
  //depois a malha fina recupera as referências originais
  agg->SetReference2(agg->NIndexes()-1);
}

void TPZAgglomerateElement::InitializeElement() {

  int indsize = NIndexes();
  //verificar se os materiais dos sub-elementos são iguais
  int i,maxdeg = -1,mat,mat2;//bc degree
  mat = FineElement(0)->Reference()->MaterialId();
  for(i=1;i<indsize;i++){
    mat2 = FineElement(i)->Reference()->MaterialId();
    if(mat2 != mat){
      for(int k=0;k<10;k++)
	PZError << "TPZAgglomerateElement::TPZAgglomerateElement data error, distinct material\n";
      exit(-1);
    }
  }
  if(indsize < 1){
    PZError << "TPZAgglomerateElement::TPZAgglomerateElement empty list\n";
    return;
  }
  //tomando o grau como o máximo grau dos sub-elementos
  for(i=0;i<indsize;i++){
    int deg = dynamic_cast<TPZCompElDisc *>(SubElement(i))->Degree();
    if(deg > maxdeg) maxdeg = deg;
  }
  SetDegree(maxdeg);
  CenterPoint();
  REAL cons = NormalizeConst();
  SetConstC(cons);
  //caso o conjunto de elementos sendo aglomerado preenchen totalmente
  //um único elemento geométrico esse será a referencia
  TPZGeoEl *ref = CalculateReference();
  TPZCompElDisc::SetReference(ref);
  if(ref) ref->SetReference(this);
}

void TPZAgglomerateElement::AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight){
  //acumula as regras de integração dos elementos aglomerados (descontínuos)
  int nsubs = NIndexes(),i;
  for(i=0; i<nsubs; i++){
    TPZAgglomerateElement *agg = dynamic_cast<TPZAgglomerateElement *>(SubElement(i));
    if(agg){//recursive call: agglomerated of the agglomerated
      agg->AccumulateIntegrationRule(degree,point,weight);
    } else {
      TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(SubElement(i));
      if(!disc){
	PZError << "TPZAgglomerateElement::AccumulateIntegrationRule, null index element\n";
	exit(-1);
      }
      //acumule integration rule
      disc->AccumulateIntegrationRule(degree,point,weight);
    }
  }
}

void TPZAgglomerateElement::CenterPoint(){

  int indsize = NIndexes(),i;
  TPZVec<REAL> volumes(indsize);
  for(i=0;i<indsize;i++){
    TPZCompEl *cel = FineElement(i);
    if(cel->Type() == EAgglomerate){//o sub-elemento é aglomerado
      volumes[i] = dynamic_cast<TPZAgglomerateElement *>(cel)->VolumeOfEl();
    } else if(cel->Type() == EDiscontinuous){//o sub é descontínuo
      volumes[i] = cel->Reference()->Volume();
    }
  }
  REAL voltot = 0.0;
  for(i=0;i<indsize;i++) voltot += volumes[i];
  REAL centx=0.0,centy=0.0,centz=0.0;
  for(i=0;i<indsize;i++){
    //o decontínuo tem fCenterPoint
    TPZCompElDisc *cel = dynamic_cast<TPZCompElDisc *>(FineElement(i));
    centx += cel->CenterPoint(0) * volumes[i];
    centy += cel->CenterPoint(1) * volumes[i];
    centz += cel->CenterPoint(2) * volumes[i];
  }
  SetCenterPoint(0,centx/voltot);
  SetCenterPoint(1,centy/voltot);
  SetCenterPoint(2,centz/voltot);
}

REAL TPZAgglomerateElement::VolumeOfEl(){

  /** 
   * um elemento aglomerado é formado de grupos de sub-elementos aglomerados
   * ou de sub-elementos descontínuos não aglomerados
   */
  REAL volume = 0.0;
  int nindex = NIndexes(),i;
  for(i=0;i<nindex;i++){
    TPZCompEl *cel = FineElement(i);
    //caso comp é aglomerado: chamada recursiva
    if(cel->Type() == EAgglomerate){//aglomerado
      volume += dynamic_cast<TPZAgglomerateElement *>(cel)->VolumeOfEl();
    } else if(cel->Type() == EDiscontinuous){
      TPZGeoEl *geo = cel->Reference();
      volume += geo->Volume();
    }
  }
  return volume;
}

TPZCompEl *TPZAgglomerateElement::FineElement(int index){

  TPZCompEl *cel = fMotherMesh->ElementVec()[fIndexes[index]];
  if(!cel) PZError << "TPZAgglomerateElement::FineElement sub element index error\n";
  return cel;
}

// a regra de integra¢ão será criada no CalcStiff e CalcRes cada ves
void TPZAgglomerateElement::CalcResidual(TPZFMatrix &Rhs,TPZCompElDisc *el){

  PZError << "TPZAgglomerateElement::CalcResidual DEVE SER IMPLEMENTADO";
}

void TPZAgglomerateElement::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){

  if(Reference()) return TPZCompElDisc::CalcStiff(ek,ef);

  if(Material() == NULL){
    cout << "TPZCompElDisc::CalcStiff : no material for this element\n";
    return;
  }
  int ncon = NConnects();
  int dim = Dimension();
  int nstate = Material()->NStateVariables();
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
  if(ncon){//pode ser no máximo ncon = 1
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
  TPZFMatrix jacinv(dim,dim);
  TPZVec<REAL> x(3,0.);
  REAL weight;
  int integ = 2*Degree();
  TPZStack<REAL> points,weights;
  //acumula as regras dos obtidos por aglomeração
  AccumulateIntegrationRule(integ,points,weights);//integra fi*fj
  int ip,npoints = weights.NElements();
  TPZVec<REAL> sol(nstate,0.);
  TPZFMatrix dsol(dim,nstate,0.);

  for(ip=0;ip<npoints;ip++){
    x[0] = points[3*ip];
    x[1] = points[3*ip+1];
    x[2] = points[3*ip+2];
    weight = weights[ip];
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
    Material()->Contribute(x,jacinv,sol,dsol,weight,axes,phix,dphix,*ek.fMat,*ef.fMat);
  }
}


TPZCompEl *TPZAgglomerateElement::SubElement(int sub){

  int nsubs = NIndexes();
  if(sub < 0  || sub > nsubs){
    PZError << "TPZAgglomerateElement::SubElement sub-element out of range\n";
    return NULL;
  }
 //  return Mesh()->ElementVec()[fIndexes[sub]];
   return fMotherMesh->ElementVec()[fIndexes[sub]];
}

void TPZAgglomerateElement::SetReference2(int sub){

  TPZCompEl *cel = SubElement(sub);
  int type = cel->Type();
  //caso comp é aglomerado: chamada recursiva
  if(type == EAgglomerate){//aglomerado
    dynamic_cast<TPZAgglomerateElement *>(cel)->SetReference();
  } else if(type == EDiscontinuous){//descontínuo
    //o geométrico agrupado apontará para o atual computacional
    cel->Reference()->SetReference(this);
  }
}

void TPZAgglomerateElement::SetReference(){

  int nindex = NIndexes(),i;
  for(i=0;i<nindex;i++){
    TPZCompEl *cel = SubElement(i);
    int type = cel->Type();
    //caso comp é aglomerado: chamada recursiva
    if(type == EAgglomerate){//aglomerado
      SetReference();
    } else if(type == EDiscontinuous){//descontínuo
      //o geométrico agrupado apontará para o atual computacional
      cel->Reference()->SetReference(this);
    }
  }
}

REAL TPZAgglomerateElement::NormalizeConst(){
  //maior distancia entre o ponto interior e os vértices do elemento
  REAL maxsub,max = -1.0;
  int nindex = NIndexes(),i;
  for(i=0;i<nindex;i++){
    TPZCompEl *cel = SubElement(i);
    //caso comp é aglomerado: chamada recursiva
    if(cel->Type() == EAgglomerate){//aglomerado
      maxsub = NormalizeConst();
    } else if(cel->Type() == EDiscontinuous){//descontínuo
      //o geométrico agrupado apontará para o atual computacional
      maxsub = dynamic_cast<TPZCompElDisc *>(cel)->NormalizeConst();
    }
    if(max < maxsub) max = maxsub;
  }
  return max;  
}


int TPZAgglomerateElement::CreateMidSideConnect(){

  if(!Material())
    PZError << "\nTPZCompElDisc::CreateMidSideConnect Material nulo\n";

  TPZStack<TPZCompElSide> list;
  int dim = Dimension();
//   int existsconnect = 0;

  if(dim == gInterfaceDimension){
    // o atual é um elemento BC
    SetConnectIndex(0, -1);
    SetDegree(-1);//=> nshape = 0
    return ConnectIndex();
  }

//   if(!existsconnect){
    //o atual é um elemento de volume e
    //não achou-se um elemento superposto
    int nvar = Material()->NStateVariables();
    int newnodeindex = Mesh()->AllocateNewConnect();
    SetConnectIndex(0,newnodeindex);
    TPZConnect &newnod = Mesh()->ConnectVec()[newnodeindex];
    int seqnum = newnod.SequenceNumber();
    Mesh()->Block().Set(seqnum,nvar*NShapeF());
    Mesh()->ConnectVec()[ConnectIndex()].IncrementElConnected();
//   }

  return ConnectIndex();
}

int TPZAgglomerateElement::Dimension() const {

//  int nind = NIndexes();
  return (gInterfaceDimension + 1);//?!
//   if(!nind){
//     PZError << "TPZAgglomerateElement::Dimension() empty list of elements\n";
//     return -1;
//   }
//   return ( fMotherMesh->ElementVec()[fIndexes[0]]->Reference()->Dimension() );
}

void TPZAgglomerateElement::Print(ostream &out) {

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

  //se cada um dos elementos aglomerados tiver o mesmo
  //número de lados retorna esse número
  //caso contrário retorna -1
  int indsize = NIndexes(),i;
  TPZStack<int> nsides;
  for(i=0;i<indsize;i++){
    TPZCompEl *cel = FineElement(i);
    if(cel->Type() == EAgglomerate){
      dynamic_cast<TPZAgglomerateElement *>(cel)->NSides();
    } else if(cel->Type() == EDiscontinuous){
      nsides.Push(cel->Reference()->NSides());
    } else {
      PZError << "TPZAgglomerateElement::NSides unknow type element\n";
    }
  }
  int nel = nsides.NElements();
  int nsd = nsides[0];
  for(i=0;i<nel;i++){
    if(nsides[i] != nsd){
      PZError << "TPZAgglomerateElement::NSides incompatible number of sides\n";
      return -1;
    }
  }
  return nsd;
}

void TPZAgglomerateElement::ListOfDiscEl(TPZStack<TPZCompEl *> &elvec){
  //agrupa na lista todos os elementos discontinuos 
  //agrupados no elemento aglomerado atual
  int indsize = NIndexes(),i;
  for(i=0;i<indsize;i++){
    TPZCompEl *cel = FineElement(i);
    if(cel->Type() == EAgglomerate){
      dynamic_cast<TPZAgglomerateElement *>(cel)->ListOfDiscEl(elvec);
    } else if(cel->Type() == EDiscontinuous){
      elvec.Push(cel);//guarda todos os descontínuos
    } else {
      PZError << "TPZAgglomerateElement::NSides unknow type element\n";
    }
  }
}

void TPZAgglomerateElement::IndexesDiscSubEls(TPZStack<int> &elvec){
  //agrupa na lista todos os indexes dos elementos discontinuos 
  //agrupados no elemento aglomerado atual
  int indsize = NIndexes(),i;
  for(i=0;i<indsize;i++){
    TPZCompEl *cel = FineElement(i);
    if(cel->Type() == EAgglomerate){
      dynamic_cast<TPZAgglomerateElement *>(cel)->IndexesDiscSubEls(elvec);
    } else if(cel->Type() == EDiscontinuous){
      elvec.Push(cel->Index());//guarda todos os descontínuos
    } else {
      PZError << "TPZAgglomerateElement::NSides unknow type element\n";
    }
  }
}

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
      TPZConservationLaw * consmat = dynamic_cast<TPZConservationLaw *>(mat);
      if ( consmat ){
	mat = consmat->NewMaterial();//mat = new TPZEulerConsLaw(*euler);//copia
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

void TPZAgglomerateElement::Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol) {


  TPZCompElDisc::Solution(qsi,var,sol);
  return;


  if(Material() == NULL){
    PZError << "TPZAgglomerateElement::Solution : no Material for this element\n";
    return;
  }
  TPZStack<int> elvec;
  TPZVec<REAL> center(3),real(3),inter(3);
  InternalPoint(inter);
  IndexesDiscSubEls(elvec);
  int size = elvec.NElements(),i,imin;
  REAL dist,mindist = 10000.0;
  //procurando sub-elemento mais próximo do ponto centro do aglomerado fCenterPoint
  for(i=0;i<size;i++){
    TPZCompEl *comp = fMotherMesh->ElementVec()[elvec[i]];
    TPZGeoEl *geo = comp->Reference();
    if(!geo) PZError <<  "TPZAgglomerateElement::Solution null reference\n";
    geo->CenterPoint(geo->NSides()-1,center);
    geo->X(center,real);
    dist = TPZGeoEl::Distance(real,inter);
    if(dist < mindist){
      mindist = dist;
      imin = i;
    }
  }
  fMotherMesh->ElementVec()[elvec[imin]]->Solution(qsi,var,sol);
}

TPZGeoEl *TPZAgglomerateElement::CalculateReference(){

  //return TPZCompElDisc::Reference();
  //porenquanto interfaces não são aglomeradas
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

  //quantos sub-elementos father contém até nlevels níveis abaixo dele
  if(!father) return 0;
  if(nlevels == 1) return father->NSubElements();
  int numberels = 0,i,lev = 1;
  while(lev < nlevels){
    int nsubs = father->NSubElements();
    numberels = nsubs;
    for(i=0;i<nsubs;i++){
      TPZGeoEl *sub = father->SubElement(i);
      if(!sub->Reference()) NSubsOfLevels(sub,lev);//só computacionais
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
      numberels += NSubCompEl(sub);//só computacionais
    } else {
      numberels++;
    }
  }
  return numberels;
}

TPZGeoEl *TPZAgglomerateElement::FatherN(TPZGeoEl *sub,int n){

  //procura-se o ancestral n níveis acima de sub
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


int Level(TPZGeoEl *gel);
void TPZAgglomerateElement::ListOfGroupings(TPZCompMesh *finemesh,TPZVec<int> &accumlist,
					    int nivel,int &numaggl,int dim){
  //todo elemento de volume deve ser agrupado nem que for para um só elemento
  finemesh->SetDimModel(dim);
  cout << "\nTPZAgglomerateElement::ListOfGroupings para malha 2D\n\n";
  int nel = finemesh->NElements(),i;
  //não todo index é sub-elemento
  accumlist.Resize(nel,-1);
  int mdim = finemesh->Dimension();
  int sdim = mdim - 1;//dimensão superficial
  for(i=0;i<nel;i++){
    TPZCompEl *cel = finemesh->ElementVec()[i];
    if(!cel) continue;
    TPZGeoEl *gel = cel->Reference();
    int type = cel->Type(),eldim = cel->Dimension();
    //agrupando elementos computacionais de volume
    if(type == EInterface) continue;//interface será clonada
    if(eldim == sdim) continue;//discontinuo BC será clonado
    TPZGeoEl *father = gel->Father();
    if(!father) continue;
    while(father && Level(father) != nivel) father = father->Father();//nova
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
  //conta o número de elementos obtidos por aglomera¢ão
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
  if(newfat != numaggl) cout << "TPZAgglomerateElement::ListOfGroupings número de pais não confere\n";
}

int Level(TPZGeoEl *gel){
  //retorna o nível do elemento gel
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

// Void Tpzagglomerateelement::Projectresidual(){
void TPZAgglomerateElement::RestrictionOperator(){

  //projeta a solução dos elementos contidos na aglomeração
  //no elemento por eles aglomerado

  if(Reference()){
    //elemento grande com geometria definida
    RestrictionOperator2();
    cout << "TPZAgglomerateElement::RestrictionOperator CONTINUA\n";
    //return;
  }
  int dimension = Dimension();
  int aggmatsize = NShapeF();
  int nvar = Material()->NStateVariables();
  TPZFMatrix aggmat(aggmatsize,aggmatsize,0.);
  TPZFMatrix loadvec(aggmatsize,nvar,0.);
  //verificar que o grau do grande é <= que o grau de cada um dos pequenos ???
  TPZStack<int> elvec;
  IndexesDiscSubEls(elvec);
  int size = elvec.NElements(),i;
  int mindegree = 1000,coarsedeg = Degree(),maxdegree = 0;
  for(i=0;i<size;i++){
    TPZCompElDisc *disc = 
      dynamic_cast<TPZCompElDisc *>(MotherMesh()->ElementVec()[elvec[i]]);
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
  //fi base do grande, uh solução sobre os pequenos
  TPZGeoEl *ref = MotherMesh()->ElementVec()[elvec[0]]->Reference();
  TPZIntPoints *intrule = 
    ref->CreateSideIntegrationRule(ref->NSides()-1,maxdegree + coarsedeg);
  //eventualmente pode criar uma regra para cada sub-elemento dentro do laço
  //de integração para reduzir o número de pontos caso os grus são distintos
  TPZFMatrix aggphix(aggmatsize,1);
  TPZFMatrix aggdphix(dimension,aggmatsize);
  TPZVec<REAL> intpoint(dimension);
  TPZFMatrix jacobian(dimension,dimension),jacinv(dimension,dimension);
  TPZFMatrix axes(3,3,0.);
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> uh(nvar);
  int npoints = intrule->NPoints();
  REAL weight,detjac;
  int in,jn,kn,ip,ind;
//   TPZBlock &fineblock = MotherMesh()->Block();
//   TPZFMatrix &FineMeshSol = MotherMesh()->Solution();
  TPZCompElDisc *disc;

  for(ind=0;ind<size;ind++){
    disc = dynamic_cast<TPZCompElDisc *>(MotherMesh()->ElementVec()[elvec[ind]]);
    ref = disc->Reference();

    for(ip=0;ip<npoints;ip++){
      intrule->Point(ip,intpoint,weight);
      ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
      ref->X(intpoint, x);
      weight *= fabs(detjac);
      Shape(x,aggphix,aggdphix);
      FineSolution(disc,aggphix,uh);
      //projetando a solução fina no elemento aglomerado
      //a soma dos detjac dos sub-elementos dá o detjac do grande
      //a geometria do grande é a soma das geometrias dos pequenos
      for(in=0; in<aggmatsize; in++) {
	for(jn=0; jn<aggmatsize; jn++) {
	  //ordem do grande <= a menor ordem dos pequenos
	  // => a regra  dos pequenos integra ok
	  aggmat(in,jn) += weight*aggphix(in,0)*aggphix(jn,0);
	}
      }
      //a soma das regras dos pequenos cobre a geometria do grande
      for(kn=0; kn<nvar; kn++) {
	loadvec(in,kn) += weight*aggphix(in,0)*uh[kn];
      }
    }
  }
  //achando a solução restringida
  aggmat.SolveDirect(loadvec,ELU);
  //transferindo a solução obtida por restrição
  int iv=0;
  TPZBlock &block = Mesh()->Block();
  TPZConnect *df = &Connect(0);
  int dfseq = df->SequenceNumber();
  int dfvar = block.Size(dfseq);
  for(jn=0; jn<dfvar; jn++) {
    block(dfseq,0,jn,0) = loadvec(iv/nvar,iv%nvar);
    iv++;
  }
}

void TPZAgglomerateElement::FineSolution(TPZCompElDisc *disc,TPZFMatrix &aggphix,TPZVec<REAL> &uh){
   
  TPZBlock &fineblock = MotherMesh()->Block();
  int nstate = disc->Material()->NStateVariables();
  TPZFMatrix &FineMeshSol = MotherMesh()->Solution();

  TPZConnect *df = &disc->Connect(0);
  int dfseq = df->SequenceNumber();
  int dfvar = fineblock.Size(dfseq);
  int pos   = fineblock.Position(dfseq);
  int iv = 0,d;
  uh.Fill(0.);
  for(d=0; d<dfvar; d++) {
    uh[iv%nstate] += aggphix(iv/nstate,0)*FineMeshSol(pos+d,0);
    iv++;
  }
}

void TPZAgglomerateElement::RestrictionOperator2(){

  cout << "TPZAgglomerateElement::RestrictionOperator2 NOT IMPLEMENTED\n";
}
