
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

TPZAgglomerateElement::TPZAgglomerateElement(int &index,TPZCompMesh &cmesh,TPZCompMesh *finemesh) : 
  TPZCompElDisc(cmesh,index) {

  /** 
   * o agomerado aponta para nulo mais o elemento computacional 
   * que ele agrupa aponta para o geométrico original
   */

  fMotherMesh = finemesh;
  CreateMidSideConnect();
  //para conferir a atualaziza¢ão
}

void TPZAgglomerateElement::AddSubElementIndex(TPZCompMesh *cmesh,int subel,int father){

  TPZAgglomerateElement *agg = dynamic_cast<TPZAgglomerateElement *>(cmesh->ElementVec()[father]);
  if(!agg){
    for(int i=0;i<20;i++)
      PZError << "TPZAgglomerateElement::AddSubElementIndex null element\n";
    return;
  }
  agg->fIndexes.Push(subel);
}

void TPZAgglomerateElement::InitializeElement(int mat) {

  int indsize = NIndexes();
  //verificar se os materiais dos sub-elementos são iguais
  int mat2,i,maxdeg = 0;
  for(i=1;i<indsize;i++){
    mat2 = SubElement(i)->Material()->Id();
    if(mat2 != mat){
      for(int i=0;i<20;i++)
	PZError << "TPZAgglomerateElement::TPZAgglomerateElement data error, distinct material\n";
      exit(-1);
    }
  }
  if(indsize < 1){
    PZError << "TPZAgglomerateElement::TPZAgglomerateElement empty list\n";
    return;
  }
  TPZMaterial *mater = Mesh()->FindMaterial(mat);
  if(mater){
    SetMaterial(mater);
  } else {
    //cria copia de material da malha fina
    mater = FineElement(0)->Material();
    if( !strcmp(mater->Name(),"TPZEulerConsLaw") ){
      TPZEulerConsLaw *euler = dynamic_cast<TPZEulerConsLaw *>(mater);
      mater = new TPZEulerConsLaw(*euler);
      SetMaterial(mater);
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
  return Mesh()->ElementVec()[fIndexes[sub]];
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
      maxsub = TPZCompElDisc::NormalizeConst();
    }
    if(max < maxsub) max = maxsub;
  }
  return max;  
}
