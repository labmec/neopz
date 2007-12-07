// -*- c++ -*-

//$Id: TPZInterfaceEl.cpp,v 1.80 2007-12-07 18:37:22 cesar Exp $

#include "pzelmat.h"
#include "TPZInterfaceEl.h"
#include "TPZCompElDisc.h"
#include "pzgeoelside.h"
#include "pzquad.h"
#include "pzmaterial.h"
#include "TPZConservationLaw.h"
#include "pzconslaw.h"
#include "pzbndcond.h"
#include "pzintel.h"
#include "pzlog.h"
#include "pzinterpolationspace.h"
#include "pzmaterialdata.h"

using namespace std;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzinterfacelement"));
#endif


void TPZInterfaceElement::SetLeftRightElements(TPZCompElSide & left, TPZCompElSide & right){

  if(fLeftElSide.Element() && fRightElSide.Element()) this->DecreaseElConnected();

  TPZCompEl * cel = left.Element();
  if(cel){

    TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement *>(cel);
    TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc *>(cel);
    if (!intel && !disc){
      PZError << __PRETTY_FUNCTION__ << " - Left element is not a TPZInterpolatedElement or TPZCompElDisc.\n";
    }

    this->fLeftElSide.SetElement( left.Element() );
    this->fLeftElSide.SetSide( left.Side() );
  }
  else{
    PZError << __PRETTY_FUNCTION__ << " - Left element is null.\n";
  }

  cel = right.Element();
  if (cel){

    TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement *>(cel);
    TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc *>(cel);
    if (!intel && !disc){
      PZError << __PRETTY_FUNCTION__ << " - Right element is not a TPZInterpolatedElement or TPZCompElDisc.\n";
    }

    this->fRightElSide.SetElement( right.Element() );
    this->fRightElSide.SetSide( right.Side() );
  }
  else{
    PZError << __PRETTY_FUNCTION__ << " - Right element is null.\n";
  }
  this->ComputeNormal();

  this->IncrementElConnected();
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

TPZInterfaceElement::~TPZInterfaceElement(){
  if(Reference()){
    this->Reference()->DecrementNumInterfaces();
    this->Reference()->ResetReference();
  }
//  this->DecreaseElConnected();
};

TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,
                                         TPZCompElSide& left, TPZCompElSide& right)
   : TPZCompEl(mesh,geo,index){

  geo->SetReference(this);
  geo->IncrementNumInterfaces();

if (left.Side() == -1 || right.Side() == -1){
  PZError << "Error at " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " Side should not be -1\n";
}

  this->SetLeftRightElements(left, right);

  this->IncrementElConnected();
}

TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index)
   : TPZCompEl(mesh,geo,index), fLeftElSide(), fRightElSide(){
  geo->SetReference(this);
  geo->IncrementNumInterfaces();
  this->IncrementElConnected();
}

TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh, const TPZInterfaceElement &copy)
   : TPZCompEl(mesh,copy) {

   this->fLeftElSide.SetElement( mesh.ElementVec()[copy.fLeftElSide.Element()->Index()] );
   this->fLeftElSide.SetSide( copy.fLeftElSide.Side() );

   this->fRightElSide.SetElement( mesh.ElementVec()[copy.fRightElSide.Element()->Index()] );
   this->fRightElSide.SetSide( copy.fRightElSide.Side() );

#ifdef DEBUG
   if( !fLeftElSide.Element() || ! fRightElSide.Element() ) {
      cout << "Something wrong with clone of interface element\n";
      exit(-1);
   }
   if(fLeftElSide.Element()->Mesh() != &mesh || fRightElSide.Element()->Mesh() != &mesh) {
      cout << "The discontinuous elements should be cloned before the interface elements\n";
      exit(-1);
   }
#endif

   fNormal = copy.fNormal;

   TPZAutoPointer<TPZMaterial> mat = copy.Material();

   this->IncrementElConnected();

   if (this->Reference()){
    this->Reference()->IncrementNumInterfaces();
   }
   else{
    PZError << "ERROR at " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " - this->Reference() is NULL\n";
    exit(-1);
   }
}


TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,
                                         const TPZInterfaceElement &copy,
                                         std::map<int,int> &gl2lcConIdx,
                                         std::map<int,int> &gl2lcElIdx) : TPZCompEl(mesh,copy)
{

  int cplftIdx = copy.fLeftElSide.Element()->Index();
  int cprgtIdx = copy.fRightElSide.Element()->Index();
  if (gl2lcElIdx.find(cplftIdx) == gl2lcElIdx.end() || gl2lcElIdx.find(cprgtIdx) == gl2lcElIdx.end())
  {
    std::stringstream sout;
    sout << "ERROR in " << __PRETTY_FUNCTION__
         << " trying to generate an interface for discontinuous elements that are not cloned."
         << " Right idx = " << cprgtIdx << " left index = " << cplftIdx;
    LOGPZ_ERROR (logger,sout.str().c_str());
    exit(-1);
  }

  this->fLeftElSide.SetElement( mesh.ElementVec()[gl2lcElIdx[cplftIdx]] );
  this->fLeftElSide.SetSide( copy.fLeftElSide.Side() );

  this->fRightElSide.SetElement( mesh.ElementVec()[gl2lcElIdx[cprgtIdx]] );
  this->fRightElSide.SetSide( copy.fRightElSide.Side() );

#ifdef DEBUG
  if( !fLeftElSide.Element() || ! fRightElSide.Element() ) {
    cout << "Something wrong with clone of interface element\n";
    exit(-1);
  }
  if(fLeftElSide.Element()->Mesh() != &mesh || fRightElSide.Element()->Mesh() != &mesh) {
    cout << "The discontinuous elements should be cloned before the interface elements\n";
    exit(-1);
  }
#endif

  fNormal = copy.fNormal;
  TPZAutoPointer<TPZMaterial> mat = copy.Material();
  this->IncrementElConnected();

  if (this->Reference()){
    this->Reference()->IncrementNumInterfaces();
  }
  else{
    PZError << "ERROR at " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " - this->Reference() is NULL\n";
    //exit(-1); //the geometric elements are generated latter...
  }
}



TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,const TPZInterfaceElement &copy,int &index)
  : TPZCompEl(mesh,copy,index) {

  //ambos elementos esquerdo e direito j�foram clonados e moram na malha aglomerada
  //o geometrico da malha fina aponta para o computacional da malha aglomerada
  fNormal = copy.fNormal;

  this->fLeftElSide.SetElement( mesh.ElementVec()[copy.fLeftElSide.Element()->Index()] );
  this->fLeftElSide.SetSide( copy.fLeftElSide.Side() );

  this->fRightElSide.SetElement( mesh.ElementVec()[copy.fRightElSide.Element()->Index()] );
  this->fRightElSide.SetSide( copy.fRightElSide.Side() );

#ifdef DEBUG
  if( !fLeftElSide.Element() || ! fRightElSide.Element() ) {
    cout << "TPZInterfaceElement::TPZInterfaceElement Something wrong with clone of interface element\n";
    exit(-1);
  }
  if(fLeftElSide.Element()->Mesh() != &mesh || fRightElSide.Element()->Mesh() != &mesh) {
    cout << "TPZInterfaceElement::TPZInterfaceElement The discontinuous elements should be cloned "
	 << "before the interface elements\n";
    exit(-1);
  }
#endif

  this->IncrementElConnected();

  if (this->Reference()){
    this->Reference()->IncrementNumInterfaces();
  }
  else{
    PZError << "ERROR at " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " - this->Reference() is NULL\n";
  }
}

TPZInterfaceElement::TPZInterfaceElement() : TPZCompEl(), fLeftElSide(), fRightElSide(),
  fNormal(3,0.)
{
   //NOTHING TO BE DONE HERE
}

TPZCompEl * TPZInterfaceElement::CloneInterface(TPZCompMesh &aggmesh,int &index, /*TPZCompElDisc **/TPZCompElSide &left, /*TPZCompElDisc **/TPZCompElSide &right) const {
   return  new TPZInterfaceElement(aggmesh, this->Reference(), index, left, right);
}

void TPZInterfaceElement::CalcResidual(TPZElementMatrix &ef){
//  cout << "\nImplementar adequadamente: " << __PRETTY_FUNCTION__ << "\n";
  TPZElementMatrix fake_ek(this->Mesh(), TPZElementMatrix::EK);
  this->CalcStiff(fake_ek, ef);
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

int TPZInterfaceElement::ConnectIndex(int i) const {

   const int nleftcon = this->NLeftConnects();
   const int nrightcon = this->NRightConnects();
   const int ncon = nleftcon + nrightcon;

   if(i < 0 || i >= ncon){
      PZError << "TPZInterfaceElement::ConnectIndex wrong argument i, i = " << i << endl;
      return -1;
   }

   if(i < nleftcon){ //required connect is associated to left neighbour
      return fLeftElSide.Element()->ConnectIndex(i);
   }

   if(i < ncon){ //required connect is associated to right neighbour
      return fRightElSide.Element()->ConnectIndex(i-nleftcon);
   }
   return -1;
}

void TPZInterfaceElement::Print(std::ostream &out){

   TPZCompEl* fLeftEl  = this->LeftElement();
   TPZCompEl* fRightEl = this->RightElement();

  out << "\nInterface element : \n";
  //out << "\tId of the geometric reference : " << Reference()->Id() << endl;
  out << "\tGeometric reference of the left element of id : ";
  if(fLeftEl){
    if(fLeftEl->Type() == EAgglomerate) out << "EAgglomerate index " << LeftElement()->Index() << endl;
    else out <<  fLeftEl->Reference()->Id() << endl;
  } else {
    out << "Null" << endl;
    cout << "TPZInterfaceElement::Print null left element\n\n";
  }
  out << "\tGeometric reference of the right element of id : ";
  if(fRightEl){
    if(fRightEl->Type() == EAgglomerate) out << "EAgglomerate index " << RightElement()->Index() << endl;
    else out << fRightEl->Reference()->Id() << endl;
  } else {
    out << "Null" << endl;
    cout << "TPZInterfaceElement::Print null right element\n\n";
  }
  out << "\tMaterial id : " << Reference()->MaterialId() << endl;

  out << "\tNormal a interface : ";
  out << "(" << fNormal[0] << "," << fNormal[1] << "," << fNormal[2] << ")\n";

}

 void TPZInterfaceElement::SetConnectIndex(int node, int index) {
   cout << "TPZInterfaceElement::SetConnectIndex should never be called\n";
 }

int TPZInterfaceElement::main(TPZCompMesh &cmesh){
  // esta func� testa o correto desempenho do algoritmo que cria e
  // deleta elementos de interface numa malha sujeita a refinamento h

  // InterfaceDimension �a dimens� do elemento de interface
  // verifica-se para cada lado de dimens� InterfaceDimension do
  // elemento que existe um elemento interface e que este �nico

  int iel,iside,nel = cmesh.NElements();

  int InterfaceDimension;

  for(iel=0;iel<nel;iel++){
    TPZCompEl *cel = cmesh.ElementVec()[iel];
    if(!cel) continue;
    TPZGeoEl *geo = cel->Reference();
    InterfaceDimension = cel->Material()->Dimension() -1;
    if(!geo){
      PZError << "TPZInterfaceElement::main computational element with null reference\n";
      exit(-1);
    }
    int nsides = geo->NSides();;
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
  int cap = list.NElements();

  if(cap){
    //caso existem elementos pequenos n� deve existir
    //interface associada ao lado atual, o lado atual
    //deve apontar para elemento computacional nulo
    TPZGeoElSide geo = comp.Reference(),neigh;
    neigh = geo.Neighbour();
    while(neigh.Exists() && neigh != geo){
      if(neigh.Element()->Reference()){
	PZError << "TPZInterfaceElement::ExistInterfaces error of data structure\n";
	exit(-1);
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
    exit(-1);
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

  int iel,nel = cmesh.NElements();
  for(iel=0;iel<nel;iel++){
    TPZCompEl *cel = cmesh.ElementVec()[iel];
    if(!cel) continue;
    if(cel->Type() != EInterface) continue;//interessa s�interfaces
    TPZGeoEl *gel = cel->Reference();
    if(!gel){
      PZError << "TPZInterfaceElement::FreeInterface computational element with null reference\n";
      exit(-1);
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

void TPZInterfaceElement::ComputeNormal(){

   TPZCompEl * fLeftEl = this->LeftElement();
   TPZCompEl * fRightEl = this->RightElement();

  //  int dim = Reference()->Dimension();
  TPZGeoEl *ref = Reference();
  int face = ref->NSides()-1;
  //face: lado do elemento bidimensional ou aresta
  //do unidimensional ou canto do ponto
  fNormal.Resize(3,0.);
  fNormal.Fill(0.);
  int faceleft,faceright;

  TPZManVector<REAL, 10> param(3),centleft(3),centright(3),point(3,0.),result(3,0.),xint(3),xvolleft(3),xvolright(3),vec(3),rib(3);
  TPZFMatrix jacobian(3,3),jacinv(3,3),axes(3,3);
  REAL detjac,normalize;
  int i;

  faceleft = fLeftEl->Reference()->NSides()-1;//lado interior do elemento esquerdo
  faceright = fRightEl->Reference()->NSides()-1; // lado interior do element direito
  fLeftEl->Reference()->CenterPoint(faceleft,centleft);//ponto centro do elemento de volume
  fRightEl->Reference()->CenterPoint(faceright,centright);
  fLeftEl->Reference()->X(centleft,xvolleft);
  fRightEl->Reference()->X(centright,xvolright);
  for(i=0;i<3;i++) vec[i] = xvolright[i]-xvolleft[i];//n� deve ser nulo

  int InterfaceDimension =  fLeftEl->Material()->Dimension() - 1;

  switch(InterfaceDimension){
  case 0:
     fNormal[0] = 1.0;// a normal sempre apontar�na dire� positiva do eixo
     fNormal[1] = 0.;
     fNormal[2] = 0.;
   break;
  case 1:
    ref->Jacobian(param,jacobian,axes,detjac,jacinv);
    for(i=0;i<3;i++) rib[i] = axes(0,i);//dire� da aresta
    this->VetorialProd(rib,vec,result);
    this->VetorialProd(result,rib,fNormal);
    //normalizando a normal
    normalize = 0.;
    for(i=0;i<3;i++) normalize += fNormal[i]*fNormal[i];
    if(normalize == 0.0)
      PZError << "TPZInterfaceElement::NormalToFace null normal vetor\n";
    normalize = sqrt(normalize);
    for(i=0;i<3;i++) fNormal[i] = fNormal[i]/normalize;
    break;
  case 2:
    ref->CenterPoint(face,param);//ponto da face
    ref->Jacobian(param,jacobian,axes,detjac,jacinv);
    for(i=0;i<3;i++) fNormal[i] = axes(2,i);
    break;
  default:
    PZError << "TPZInterfaceElement::NormalToFace in case that not treated\n";
    fNormal.Resize(0);
    return;
  }

  //to guarantee the normal points from left to right neighbours:
  normalize = 0.;
  for(i=0; i<3; i++) normalize += fNormal[i]*vec[i];
  if(normalize < 0.) {
    for(i=0; i<3; i++) fNormal[i] = -fNormal[i];
  }
}

void TPZInterfaceElement::VetorialProd(TPZVec<REAL> &ivet,TPZVec<REAL> &jvet,TPZVec<REAL> &kvet){

  kvet.Resize(3);
  kvet[0] =  ivet[1]*jvet[2] - ivet[2]*jvet[1];
  kvet[1] = -ivet[0]*jvet[2] + ivet[2]*jvet[0];
  kvet[2] =  ivet[0]*jvet[1] - ivet[1]*jvet[0];
}

void TPZInterfaceElement::Normal(TPZVec<REAL> &normal) {
  normal.Resize(3);
  for(int i=0;i<3;i++) normal[i] = fNormal[i];
}

void TPZInterfaceElement::EvaluateError(void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
  TPZVec<REAL> &errors, TPZBlock * /*flux */) {
   errors.Fill(0.0);
}


  /**
  * returns the unique identifier for reading/writing objects to streams
  */
int TPZInterfaceElement::ClassId() const
{
  return TPZINTERFACEELEMENTID;
}

template class
    TPZRestoreClass< TPZInterfaceElement, TPZINTERFACEELEMENTID>;

  /**
  Save the element data to a stream
  */
void TPZInterfaceElement::Write(TPZStream &buf, int withclassid)
{
  TPZCompEl::Write(buf,withclassid);
  int leftelindex = fLeftElSide.Element()->Index();
  int rightelindex = fRightElSide.Element()->Index();
  if ( (this->Index() < leftelindex) || (this->Index() < rightelindex) ){
     PZError << __PRETTY_FUNCTION__ << endl
	     << "Indices of neighbours are less than interface index:" << endl
	     << "Left: " << leftelindex << ", Right: " << rightelindex << ", this: " << this->Index() << endl;
  }

  int leftside = fLeftElSide.Side();
  int rightside = fRightElSide.Side();

  buf.Write(&leftelindex,1);
  buf.Write(&leftside,1);
  buf.Write(&rightelindex,1);
  buf.Write(&rightside,1);
  WriteObjects(buf,fNormal);
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
   }
  int leftelindex;
  int rightelindex;
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

  ReadObjects(buf,fNormal);
}

void TPZInterfaceElement::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){

  TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(Material().operator ->());
  if(!mat || mat->Name() != "no_name"){
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

  TPZMaterialData data;
  const int dim = this->Dimension();
  const int diml = left->Dimension();
  const int dimr = right->Dimension();
  int nshapel = left ->NShapeF();
  int nshaper = right->NShapeF();
  const int nstatel = left->Material()->NStateVariables();
  const int nstater = right->Material()->NStateVariables();
  this->InitMaterialData(data,left,right);

   TPZManVector<TPZConnect*> ConnectL, ConnectR;
   TPZManVector<int> ConnectIndexL, ConnectIndexR;

   this->GetConnects( this->LeftElementSide(),  ConnectL, ConnectIndexL );
   this->GetConnects( this->RightElementSide(), ConnectR, ConnectIndexR );
   const int ncon = ConnectL.NElements() + ConnectR.NElements();
   const int neql = nshapel * nstatel;
   const int neqr = nshaper * nstater;
   const int neq = neql + neqr;
   ek.fMat.Redim(neq,neq);
   ef.fMat.Redim(neq,1);
   ek.fBlock.SetNBlocks(ncon);
   ef.fBlock.SetNBlocks(ncon);
   ek.fConnect.Resize(ncon);
   ef.fConnect.Resize(ncon);

   int ic = 0;
   int n = ConnectL.NElements();
   for(int i = 0; i < n; i++) {
    const int nshape = left->NConnectShapeF(i);
    const int con_neq = nstatel * nshape;
    ek.fBlock.Set(ic,con_neq );
    ef.fBlock.Set(ic,con_neq);
    (ef.fConnect)[ic] = ConnectIndexL[i];
    (ek.fConnect)[ic] = ConnectIndexL[i];
    ic++;
   }
   n = ConnectR.NElements();
   for(int i = 0; i < n; i++) {
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

   //LOOKING FOR MAX INTERPOLATION ORDER
   data.leftp = left->MaxOrder();
   data.rightp = right->MaxOrder();
   //Max interpolation order
   const int p = (data.leftp > data.rightp) ? data.leftp : data.rightp;

   TPZGeoEl *ref = Reference();
   TPZIntPoints *intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, 2*(p+1) );
   if(mat->HasForcingFunction()){
      TPZManVector<int,10> order(3);
      intrule->GetOrder(order);
      int maxorder = intrule->GetMaxOrder();
      order.Fill(maxorder);
      intrule->SetOrder(order);
   }
   const int npoints = intrule->NPoints();

   //integration points in left and right elements: making transformations to neighbour elements
   TPZTransform TransfLeft, TransfRight;
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

#ifdef DEBUG
      this->CheckConsistencyOfMappedQsi(this->LeftElementSide(), intpoint, LeftIntPoint);
      this->CheckConsistencyOfMappedQsi(this->RightElementSide(), intpoint, RightIntPoint);
#endif

      left->ComputeShape(LeftIntPoint, data.x, data.leftjac, data.axesleft, data.leftdetjac, data.leftjacinv, data.phil, data.dphixl);
      right->ComputeShape(RightIntPoint, data.x, data.rightjac, data.axesright, data.rightdetjac, data.rightjacinv, data.phir, data.dphixr);

      this->ComputeRequiredData(data, left, right, intpoint, LeftIntPoint, RightIntPoint);
      mat->ContributeInterface(data, weight, ek.fMat, ef.fMat);

   }//loop over integration points

   delete intrule;
}

void TPZInterfaceElement::GetConnects(TPZCompElSide &elside, TPZVec<TPZConnect*> &connects, TPZVec<int> &connectindex){

   TPZCompEl * el = elside.Element();

   if(el){
      int ncon = el->NConnects();
      connects.Resize(ncon);
      connects.Fill(NULL);
      connectindex.Resize(ncon);
      connectindex.Fill(-1);
      int i, index;
      for(i = 0; i < ncon; i++){
	 index = el->ConnectIndex(i);
	 connectindex[i] = index;
	 connects[i] = &(fMesh->ConnectVec()[ index ]);
      }//for

   }
   else{   //if (!el)
      connects.Resize(0);
      connectindex.Resize(0);
   }

}//end of method

void TPZInterfaceElement::EvaluateInterfaceJumps(TPZVec<REAL> &errors){
   TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(Material().operator ->());
   if(!mat || mat->Name() != "no_name"){
      PZError << "TPZInterfaceElement::CalcStiff interface material null, do nothing\n";
      return;
   }

   int NErrors = this->Material()->NEvalErrors();
   errors.Resize(NErrors);
   errors.Fill(0.0);

   TPZMaterialData data;
   int nstatel = this->LeftElement()->Material()->NStateVariables();
   int nstater = this->RightElement()->Material()->NStateVariables();
   int diml = this->LeftElement()->Dimension();
   int dimr = this->RightElement()->Dimension();

   //LOOKING FOR MAX INTERPOLATION ORDER
   TPZGeoEl *ref = Reference();
   TPZIntPoints *intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, 2 );
   TPZManVector<int> order(3);
   intrule->GetOrder(order);
   int maxorder = intrule->GetMaxOrder();
   order.Fill(maxorder);
   intrule->SetOrder(order);
   const int npoints = intrule->NPoints();
   TPZManVector<REAL> intpoint(3), x(3);
   REAL weight;
   //LOOP OVER INTEGRATION POINTS
   for(int ip = 0; ip < npoints; ip++){
      intrule->Point(ip,intpoint,weight);
      ref->Jacobian( intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
      weight *= fabs(data.detjac);
      ref->X(intpoint, x);

      //method NeighbourSolution will compute the transformation in this->MapQsi every time it is called
      //(which means for all integration point). Instead of calling NeighbourSolution the whole method
      //may be written here but keeping the transformation computed at first integration point (as done in CalcStiff).
      this->NeighbourSolution(this->LeftElementSide(), intpoint, data.soll, data.dsoll, data.axesleft);
      this->NeighbourSolution(this->RightElementSide(), intpoint, data.solr, data.dsolr, data.axesright);

      TPZManVector<REAL> leftNormalDeriv(nstatel), rightNormalDeriv(nstater), normal;
      this->Normal(normal);

      if (data.soll.NElements()){
        for(int iv = 0; iv < nstatel; iv++){
          leftNormalDeriv[iv] = 0.;
          for(int d = 0; d < diml; d++){
            leftNormalDeriv[iv]  += data.dsoll(d,iv)* normal[d];
          }//for d
        }//for iv
      }//if

      if (data.solr.NElements()){
        for(int iv = 0; iv < nstater; iv++){
          rightNormalDeriv[iv] = 0.;
          for(int d = 0; d < dimr; d++){
            rightNormalDeriv[iv] += data.dsolr(d,iv)* normal[d];
          }//for d
        }  //for iv
      }//if

      TPZManVector<REAL> localerror(NErrors);
      mat->InterfaceJumps(data.x, data.soll, leftNormalDeriv, data.solr, rightNormalDeriv, localerror);

      for(int ier = 0; ier < NErrors; ier++){
        errors[ier] += localerror[ier]*weight;
      }
   }//loop over integration points
   delete intrule;
   //Norma sobre o elemento
   for(int ier = 0; ier < NErrors; ier++)
   errors[ier] = sqrt(errors[ier]);

}//method

void TPZInterfaceElement::ComputeError(int errorid,
                                       TPZVec<REAL> &errorL,
                                       TPZVec<REAL> &errorR){

   TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(Material().operator ->());
   if(!mat){
      PZError << "TPZInterfaceElement::CalcStiff interface material null, do nothing\n";
      return;
   }


   TPZInterpolationSpace * left = dynamic_cast<TPZInterpolationSpace*>(this->LeftElement());
   TPZInterpolationSpace * right = dynamic_cast<TPZInterpolationSpace*>(this->RightElement());

   if (!left || !right){
     PZError << "\nError at TPZInterfaceElement::CalcStiff null neighbour\n";
     return;
   }
   if(!left->Material() || !right->Material()){
      PZError << "\n Error at TPZInterfaceElement::CalcStiff null material\n";
      return;
   }

  TPZMaterialData data;
  data.SetAllRequirements(true);
  this->InitMaterialData(data,left,right);

//   data.fPrimalExactSol = fp;
//   data.fDualExactSol = fd;


  const int dim = this->Dimension();
  const int diml = left->Dimension();
  const int dimr = right->Dimension();

   //LOOKING FOR MAX INTERPOLATION ORDER
   data.leftp = left->MaxOrder();
   data.rightp = right->MaxOrder();
   //Max interpolation order
   const int p = (data.leftp > data.rightp) ? data.leftp : data.rightp;

   TPZGeoEl *ref = Reference();
   TPZIntPoints *intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, 2*(p+1) );
   if(mat->HasForcingFunction()){
      TPZManVector<int,10> order(3);
      intrule->GetOrder(order);
      int maxorder = intrule->GetMaxOrder();
      order.Fill(maxorder);
      intrule->SetOrder(order);
   }
   const int npoints = intrule->NPoints();

   //integration points in left and right elements: making transformations to neighbour elements
   TPZTransform TransfLeft, TransfRight;
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

#ifdef DEBUG
      this->CheckConsistencyOfMappedQsi(this->LeftElementSide(), intpoint, LeftIntPoint);
      this->CheckConsistencyOfMappedQsi(this->RightElementSide(), intpoint, RightIntPoint);
#endif

      left->ComputeShape(LeftIntPoint, data.x, data.leftjac, data.axesleft, data.leftdetjac, data.leftjacinv, data.phil, data.dphixl);
      right->ComputeShape(RightIntPoint, data.x, data.rightjac, data.axesright, data.rightdetjac, data.rightjacinv, data.phir, data.dphixr);

      this->ComputeRequiredData(data, left, right, intpoint, LeftIntPoint, RightIntPoint);

      mat->ContributeInterfaceErrors(data, weight,errorL,errorR,errorid);

   }//loop over integration points

   delete intrule;
}

void TPZInterfaceElement::Integrate(int variable, TPZVec<REAL> & value){
  const int varsize = this->Material()->NSolutionVariables(variable);
  value.Resize(varsize);
  value.Fill(0.);
}

void TPZInterfaceElement::IntegrateInterface(int variable, TPZVec<REAL> & value){
  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " : no material for this element\n";
    return;
  }
  if (!this->Reference()){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " : no reference element\n";
    return;
  }
  const int dim = this->Dimension();
  TPZInterpolationSpace *left = dynamic_cast<TPZInterpolationSpace*>(this->LeftElement());
  TPZInterpolationSpace *right = dynamic_cast<TPZInterpolationSpace*>(this->RightElement());
  if (!left || !right){
    PZError << "\nError at TPZInterfaceElement::CalcStiff null neighbour\n";
    return;
  }
  if(!left->Material() || !right->Material()){
    PZError << "\n Error at TPZInterfaceElement::CalcStiff null material\n";
    return;
  }

  ///local variables
  REAL weight;
  TPZMaterialData data;
  this->InitMaterialData(data,left,right);
  TPZGeoEl *ref = Reference();
  data.leftp = left->MaxOrder();
  data.rightp = right->MaxOrder();
  TPZManVector<REAL, 3> intpoint(dim,0.);
  const int varsize = material->NSolutionVariables(variable);
  ///Max interpolation order
  const int p = (data.leftp > data.rightp) ? data.leftp : data.rightp;

  ///Integration rule
  TPZIntPoints *intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, 2*(p+1));
  if(material->HasForcingFunction()){
    TPZManVector<int,10> order(3);
    intrule->GetOrder(order);
    int maxorder = intrule->GetMaxOrder();
    order.Fill(maxorder);
    intrule->SetOrder(order);
  }

  ///loop over integration points
  const int npoints = intrule->NPoints();
  int ip, iv;
  value.Resize(varsize);
  value.Fill(0.);
  TPZManVector<REAL> locval(varsize);
  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    ref->Jacobian(intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
    weight *= fabs(data.detjac);
    this->ComputeSolution(intpoint, data.phi, data.dphix, data.axes, data.sol, data.dsol);
    this->NeighbourSolution(this->LeftElementSide(), intpoint, data.soll, data.dsoll, data.axesleft);
    this->NeighbourSolution(this->RightElementSide(), intpoint, data.solr, data.dsolr, data.axesright);
    material->Solution(data, variable, locval);
    for(iv = 0; iv < varsize; iv++){
      value[iv] += locval[iv]*weight;
    }///for iv
  }///for ip
  delete intrule;
}///method

void TPZInterfaceElement::ComputeSideTransform(TPZCompElSide &Neighbor, TPZTransform &transf){
  TPZGeoEl * neighel = Neighbor.Element()->Reference();
  const int dim = this->Dimension();
  TPZTransform LocalTransf(dim);
  TPZGeoElSide thisgeoside(this->Reference(), this->Reference()->NSides()-1);
  TPZGeoElSide neighgeoside(neighel, Neighbor.Side());
  thisgeoside.SideTransform3(neighgeoside, LocalTransf);

  TPZGeoElSide highdim(neighel, neighel->NSides()-1);
  transf = neighgeoside.SideToSideTransform(highdim).Multiply(LocalTransf);
}//ComputeSideTransform

void TPZInterfaceElement::MapQsi(TPZCompElSide &Neighbor, TPZVec<REAL> &qsi, TPZVec<REAL> &NeighIntPoint){
  TPZTransform Transf;
  this->ComputeSideTransform(Neighbor, Transf);
  Transf.Apply( qsi, NeighIntPoint );
#ifdef DEBUG
  this->CheckConsistencyOfMappedQsi(Neighbor, qsi, NeighIntPoint);
#endif
}//MapQsi

bool TPZInterfaceElement::CheckConsistencyOfMappedQsi(TPZCompElSide &Neighbor, TPZVec<REAL> &qsi, TPZVec<REAL>&NeighIntPoint){
  const REAL tol = 1.e-10;
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
    return false;
  }
  return true;
}//void

void TPZInterfaceElement::NeighbourSolution(TPZCompElSide & Neighbor, TPZVec<REAL> & qsi,
                                            TPZVec<REAL> &sol, TPZFMatrix &dsol,
                                            TPZFMatrix &NeighborAxes){
  TPZGeoEl * neighel = Neighbor.Element()->Reference();
  const int neighdim = neighel->Dimension();
  TPZManVector<REAL,3> NeighIntPoint(neighdim);
  this->MapQsi(Neighbor, qsi, NeighIntPoint);
  Neighbor.Element()->ComputeSolution(NeighIntPoint, sol, dsol, NeighborAxes);
}

void TPZInterfaceElement::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
                                          const TPZFMatrix &axes, TPZVec<REAL> &sol, TPZFMatrix &dsol){
  sol.Resize(0);
  dsol.Resize(0,0);
}

void TPZInterfaceElement::ComputeSolution(TPZVec<REAL> &qsi,
                                          TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix &axes){
  sol.Resize(0);
  dsol.Resize(0,0);
  axes.Zero();
}

void TPZInterfaceElement::ComputeSolution(TPZVec<REAL> &qsi,
                                 TPZVec<REAL> &normal,
                                 TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
                                 TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes){
  normal = this->fNormal;
  this->NeighbourSolution(this->fLeftElSide, qsi, leftsol, dleftsol, leftaxes);
  this->NeighbourSolution(this->fRightElSide, qsi, rightsol, drightsol, rightaxes);
}//method

void TPZInterfaceElement::InitMaterialData(TPZMaterialData &data, TPZInterpolationSpace *left, TPZInterpolationSpace *right){
  TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(Material().operator ->());
  if (!mat){
    PZError << "FATAL ERROR AT "  << __PRETTY_FUNCTION__ << "\n";
  }
  mat->FillDataRequirementsInterface(data);
  int nshapel = left->NShapeF();
  int nshaper = right->NShapeF();
  const int nstatel = left->Material()->NStateVariables();
  const int nstater = right->Material()->NStateVariables();
  const int dim = this->Dimension();
  const int diml = left->Dimension();
  const int dimr = right->Dimension();
  data.phil.Redim(nshapel,1);
  data.dphixl.Redim(diml,nshapel);
  data.phir.Redim(nshaper,1);
  data.dphixr.Redim(dimr,nshaper);
  data.axes.Redim(3,3);
  data.axesleft.Redim(3,3);
  data.axesright.Redim(3,3);
  data.jacobian.Redim(dim,dim);
  data.leftjac.Redim(diml,diml);
  data.rightjac.Redim(dimr,dimr);
  data.jacinv.Redim(dim,dim);
  data.leftjacinv.Redim(diml,diml);
  data.rightjacinv.Redim(dimr,dimr);
  data.x.Resize(3);
  if (data.fNeedsNeighborSol){
    data.soll.Resize(nstatel);
    data.solr.Resize(nstater);
    data.dsoll.Redim(diml,nstatel);
    data.dsolr.Redim(dimr,nstater);
  }
  data.normal = this->fNormal;
}//void

void TPZInterfaceElement::ComputeRequiredData(TPZMaterialData &data,
                                              TPZInterpolationSpace *left, TPZInterpolationSpace *right,
                                              TPZVec<REAL> &qsi,
                                              TPZVec<REAL> &LeftIntPoint, TPZVec<REAL> &RightIntPoint){
  if (data.fNeedsNeighborSol){
    left->ComputeSolution(  LeftIntPoint, data.phil, data.dphixl, data.axesleft, data.soll, data.dsoll );
    right->ComputeSolution( RightIntPoint,data.phir, data.dphixr, data.axesright, data.solr, data.dsolr );
  }

  if (data.fNeedsSol){
    this->ComputeSolution(qsi, data.phi, data.dphix, data.axes, data.sol, data.dsol);
    //this->ComputeSolution(qsi, data.sol, data.dsol, data.axes);//chamando acima porque senao nao chega
                                                                 //a TPZReferredCompEl<TPZInterfaceElement>
  }

  if (data.fNeedsHSize){
    const int dim = this->Dimension();
    REAL faceSize;
    if (dim == 0){//it means I am a point
      //2*(a+b)/2
      faceSize = left->InnerRadius() + right->InnerRadius();
    }
    else{
      faceSize = 2.*this->Reference()->ElementRadius();//Ivo Mozolevski's suggestion. It works well for elements with small aspect ratio
    }
    data.HSize = faceSize;
  }
}//void

