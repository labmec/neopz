//$Id: TPZAgglomerateEl.cpp,v 1.7 2003-10-23 15:44:04 tiago Exp $

#include "TPZAgglomerateEl.h"
#include "TPZInterfaceEl.h"
#include "TPZEulerConsLaw.h"
#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzquad.h"
#include "pzelmat.h"

TPZAgglomerateElement::TPZAgglomerateElement(int nummat,int &index,TPZCompMesh &aggcmesh,TPZCompMesh *finemesh) : 
  TPZCompElDisc(aggcmesh,index) {

  /** 
   * o algomerado aponta para nulo mais o elemento computacional 
   * que ele agrupa aponta para o geométrico original
   */
  fMotherMesh = finemesh;
  TPZMaterial *mat = finemesh->FindMaterial(nummat);
  SetMaterial(mat);
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
  TPZCompEl *cel = FineElement(0);
  int matid =   cel->Reference()->MaterialId();
  TPZMaterial *mater = Mesh()->FindMaterial(matid),*material=0;
  if(mater){
    SetMaterial(mater);
  } else {
    //cria copia de material da malha fina
    //só entra algumas vezes por aqui
    mater = FineElement(0)->Material();
    if(matid < 0){
      mater = dynamic_cast<TPZBndCond *>(mater)->Material();
      material = Mesh()->FindMaterial(mater->Id());//existe a copia?
    }
    if( !strcmp(mater->Name(),"TPZEulerConsLaw") && !material ){
      TPZEulerConsLaw *euler = dynamic_cast<TPZEulerConsLaw *>(mater);
      material = new TPZEulerConsLaw(*euler);
      SetMaterial(material);
    }
    if(matid < 0){
      //aqui material é a copia do material
      TPZBndCond *cc = dynamic_cast<TPZBndCond *>( MotherMesh()->FindMaterial(matid) );
      TPZBndCond copy(*cc,material);
      SetMaterial(&copy);
    }
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
  //a partir daqui os sub-elementos conhecem o aglomerado
  //depois a malha fina recupera as referências originais
  for(i=0;i<indsize;i++) SetReference();
}

void TPZAgglomerateElement::AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight){

  cout << "TPZAgglomerateElement::AccumulateIntegrationRule NOT IMPLEMENTED\n\n";
  int nsubs = NIndexes(),i;
  for(i=0; i<nsubs; i++){
    TPZCompElDisc *agg = dynamic_cast<TPZCompElDisc *>(SubElement(i));
    agg->AccumulateIntegrationRule(degree,point,weight);
  }
}

void TPZAgglomerateElement::CenterPoint(){

  int indsize = NIndexes(),i;
  TPZVec<REAL> volumes(indsize);
  for(i=0;i<indsize;i++){
    TPZCompEl *cel = FineElement(i);
    if(cel->Type() == 17){//o sub-elemento é aglomerado
      volumes[i] = dynamic_cast<TPZAgglomerateElement *>(cel)->VolumeOfEl();
    } else if(cel->Type() == 15){//o sub é descontínuo
      volumes[i] = cel->Reference()->Volume();
    }
  }
  REAL voltot = 0.0;
  for(i=0;i<indsize;i++) voltot += volumes[i];
  REAL centx=0.0,centy=0.0,centz=0.0;
  for(i=0;i<indsize;i++){
    //o aglomerado tem fCenterPoint
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
    if(cel->Type() == 17){//aglomerado
      volume += dynamic_cast<TPZAgglomerateElement *>(cel)->VolumeOfEl();
    } else if(cel->Type() == 15){
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

// a regra de integra¢ão será criada no CalcStiff e CalcRes cada ves ??????????????????????
void TPZAgglomerateElement::CalcResidual(TPZFMatrix &Rhs,TPZCompElDisc *el){

  PZError << "TPZAgglomerateElement::CalcResidual DEVE SER IMPLEMENTADO";
}

void TPZAgglomerateElement::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){

  PZError << "TPZAgglomerateElement::CalcStiff DEVE SER IMPLEMENTADO";
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

void TPZAgglomerateElement::SetReference(){

  int nindex = NIndexes(),i;
  for(i=0;i<nindex;i++){
    TPZCompEl *cel = SubElement(i);
    //caso comp é aglomerado: chamada recursiva
    if(cel->Type() == 17){//aglomerado
      SetReference();
    } else if(cel->Type() == 15){//descontínuo
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
    if(cel->Type() == 17){//aglomerado
      maxsub = NormalizeConst();
    } else if(cel->Type() == 15){//descontínuo
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
    TPZConnect &newnod = Mesh()->ConnectVec()[newnodeindex];
    int seqnum = newnod.SequenceNumber();
    Mesh()->Block().Set(seqnum,nvar*NShapeF());
    SetConnectIndex(0,newnodeindex);
    Mesh()->ConnectVec()[ConnectIndex()].IncrementElConnected();
//   }

  return ConnectIndex();
}

int TPZAgglomerateElement::Dimension() const {

//  int nind = NIndexes();
  return (gInterfaceDimension + 1);
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
