// -*- c++ -*-

//$Id: TPZInterfaceEl.cpp,v 1.47 2005-12-19 12:00:27 tiago Exp $

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

int TPZInterfaceElement::gCalcStiff = 1;
using namespace std;

void TPZInterfaceElement::SetLeftRightElements(TPZCompElSide & left, TPZCompElSide & right){

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
}//method

void TPZInterfaceElement::IncrementElConnected(){
   const int ncon = this->NConnects();
   for(int i = 0; i < ncon; i++){
      int index = this->ConnectIndex(i);
      fMesh->ConnectVec()[index].IncrementElConnected();
   }
}

/**
 * Para CloneInterface.
 * A normal �clonada, e n� recalculada.
 */
TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompEl *left,TPZCompEl *right, const TPZVec<REAL> & normal)
   : TPZCompEl(mesh,geo,index), fNormal(3,0.) {

  geo->SetReference(this);
  int materialid = geo->MaterialId();

  //poderia eliminar esta vari�el e carrega-la do elemento de volume associado
  fMaterial = mesh.FindMaterial(materialid);
  
  TPZCompElSide leftside(left, -1);
  TPZCompElSide rightside(right, -1);
  this->SetLeftRightElements(leftside, rightside);

  for (int i = 0; i < fNormal.NElements(); i++)
    fNormal[i] = normal[i];

  this->IncrementElConnected();

}

//construtor para o elemento descontinuo
TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,
 					 TPZCompEl *left,TPZCompEl *right)
   : TPZCompEl(mesh,geo,index){
//    TPZManVector<REAL,3> NORMAL(3, 0.0);
//    TPZInterfaceElement::TPZInterfaceElement(mesh,geo,index,left,right,NORMAL);

  geo->SetReference(this);
  int materialid = geo->MaterialId();

  //poderia eliminar esta vari�el e carrega-la do elemento de volume associado
  fMaterial = mesh.FindMaterial(materialid);
  
  TPZCompElSide leftside(left, -1);
  TPZCompElSide rightside(right, -1);
  this->SetLeftRightElements(leftside, rightside);

  this->NormalToFace(fNormal);

  this->IncrementElConnected();
}

//construtor para o elemento continuo/descontinuo
TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,
 					 TPZCompEl *left,TPZCompEl *right, int leftside, int rightside)
   : TPZCompEl(mesh,geo,index){

  geo->SetReference(this);
  int materialid = geo->MaterialId();

  //poderia eliminar esta vari�el e carrega-la do elemento de volume associado
  fMaterial = mesh.FindMaterial(materialid);
  TPZCompElSide leftCompElside(left, leftside);
  TPZCompElSide rightCompElside(right, rightside);
  this->SetLeftRightElements(leftCompElside, rightCompElside);

  this->NormalToFace(fNormal);

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

   TPZMaterial *mat = copy.Material();
   if(mat) {
      int materialid = mat->Id();
      fMaterial = mesh.FindMaterial(materialid);
   } else {
      fMaterial = 0;
   }

   this->IncrementElConnected();
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

  TPZMaterial *mat = copy.Material();
  if(mat) {
    int materialid = mat->Id();
    fMaterial = mesh.FindMaterial(materialid);
  } else {
    fMaterial = 0;
  }

  this->IncrementElConnected();
}

TPZInterfaceElement::TPZInterfaceElement() : TPZCompEl(), fLeftElSide(), fRightElSide(),
  fNormal(3,0.), fMaterial(0)
{
   //NOTHING TO BE DONE HERE
}

TPZCompEl * TPZInterfaceElement::CloneInterface(TPZCompMesh &aggmesh,int &index, TPZCompElDisc *left, TPZCompElDisc *right) const {
   return  new TPZInterfaceElement(aggmesh, this->Reference(), index, left, right, this->fNormal);
}

void TPZInterfaceElement::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){

   switch (TPZInterfaceElement::gCalcStiff)
   {
      case 1 :
	 this->CalcStiffStandard(ek, ef);
	 break;

      case 2 :
	 this->CalcStiffPenalty(ek, ef);
	 break;

      case 3 :
	 this->CalcStiffContDisc(ek, ef);
	 break;

      default:
	 PZError << "TPZInterfaceElement::CalcStiff - CalcStiff method not implemented." << endl;
   }

}

void TPZInterfaceElement::CalcResidual(TPZElementMatrix &ef){

   switch (TPZInterfaceElement::gCalcStiff)
   {
      case 1 :
	 this->CalcResidualStandard(ef);
	 break;

	 PZError << "TPZInterfaceElement::CalcStiff - CalcStiff method not implemented." << endl;
   }

}

void TPZInterfaceElement::CalcStiffStandard(TPZElementMatrix &ek, TPZElementMatrix &ef){


  TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial);
#ifndef NODEBUG
  if(!mat || !strcmp("no_name",mat->Name())){
    PZError << "TPZInterfaceElement::CalcStiff interface material null, do nothing\n";
    return;
  }
#endif

  TPZCompElDisc  *left = dynamic_cast<TPZCompElDisc*>( LeftElement() );
  TPZCompElDisc *right = dynamic_cast<TPZCompElDisc*>( RightElement() );

#ifndef NODEBUG
  if(!left->Material() || !right->Material()){
    PZError << "TPZInterfaceElement::CalcStiff null material\n";
    return;
  }
#endif

  int fConnectIndexL = -1;
  TPZConnect* fConnectL = NULL;
  if (left->NConnects()){
     fConnectIndexL = left->ConnectIndex();
     fConnectL = &( left->Connect(0) );//Discontinuous elements has only 1 connect
  }

  int fConnectIndexR = -1;
  TPZConnect* fConnectR = NULL;
  if (right->NConnects()){
     fConnectIndexR = right->ConnectIndex();
     fConnectR = &( right->Connect(0) );//Discontinuous elements has only 1 connect
  }

  

  int nshapel = left->NShapeF();
  int nshaper = right->NShapeF();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix &MeshSol = Mesh()->Solution();
  int nstatel = left->Material()->NStateVariables();
  int nstater = right->Material()->NStateVariables();
  int neql = nshapel * nstatel;
  int neqr = nshaper * nstater;
  int dim = Dimension();
  int diml = left->Dimension();
  int dimr = right->Dimension();
  int ncon = 0;
  if(fConnectL) ncon++;
  if(fConnectR) ncon++;


  int neq = neql + neqr;
  ek.fMat.Redim(neq,neq);
  ef.fMat.Redim(neq,1);
  int ic = 0;
  ek.fBlock.SetNBlocks(ncon);
  ef.fBlock.SetNBlocks(ncon);
  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);
  if(fConnectL) {
      ek.fBlock.Set(ic,neql);
      ef.fBlock.Set(ic,neql);
    (ef.fConnect)[ic] = fConnectIndexL;
    (ek.fConnect)[ic] = fConnectIndexL;
      ic++;
  }
  if(fConnectR) {
    ek.fBlock.Set(ic,neqr);
    ef.fBlock.Set(ic,neqr);
    (ef.fConnect)[ic] = fConnectIndexR;
    (ek.fConnect)[ic] = fConnectIndexR;
  }
  ek.fBlock.Resequence();
  ef.fBlock.Resequence();

  TPZFNMatrix<100> phixl(nshapel,1),dphixl(diml,nshapel);
  TPZFNMatrix<100> phixr(nshaper,1),dphixr(dimr,nshaper);
  TPZFNMatrix<9> axes(3,3);
  TPZFNMatrix<9> jacobian(dim,dim);
  TPZFNMatrix<9> jacinv(dim,dim);
  TPZManVector<REAL,3> x(3);
  TPZManVector<REAL,3> intpoint(dim);
  REAL detjac,weight;
  TPZManVector<REAL,5> soll(nstatel),solr(nstater);
  TPZFNMatrix<15> dsoll(diml,nstatel),dsolr(dimr,nstater);
  int pl = left->Degree();
  int pr = right->Degree();
  int p = (pl > pr) ? pl : pr;
  TPZGeoEl *ref = Reference();
  int face = ref->NSides()-1;
  TPZIntPoints *intrule = ref->CreateSideIntegrationRule(face,2*p);//integra u(n)*fi
  if(fMaterial->HasForcingFunction())
  {
  	TPZManVector<int> order(3);
	intrule->GetOrder(order);
	int maxorder = intrule->GetMaxOrder();
	order.Fill(maxorder);
	intrule->SetOrder(order);
  }
  int npoints = intrule->NPoints();
  int ip;
  //  TPZManVector<REAL,3> point(3);
//  TPZBndCond *bcleft = 0,*bcright=0;

  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    weight *= fabs(detjac);
    ref->X(intpoint, x);
    //solu� da itera� anterior
    if(fConnectL){
      left->Shape(x,phixl,dphixl);
      soll.Fill(0.);
      dsoll.Zero();
      TPZConnect *df = fConnectL;
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      int pos = block.Position(dfseq);
      int iv = 0,d;
      for(int jn=0; jn<dfvar; jn++) {
	soll[iv%nstatel] += phixl(iv/nstatel,0)*MeshSol(pos+jn,0);
	for(d=0; d<diml; d++)
	  dsoll(d,iv%nstatel) += dphixl(d,iv/nstatel)*MeshSol(pos+jn,0);
	iv++;
      }
    }
    //solu� da itera� anterior
    if(fConnectR){
      right->Shape(x,phixr,dphixr);
      solr.Fill(0.);
      dsolr.Zero();
      TPZConnect *df = fConnectR;
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      int pos = block.Position(dfseq);
      int iv = 0,d;
      for(int jn=0; jn<dfvar; jn++) {
	solr[iv%nstater] += phixr(iv/nstater,0)*MeshSol(pos+jn,0);
	for(d=0; d<dimr; d++)
	  dsolr(d,iv%nstater) += dphixr(d,iv/nstater)*MeshSol(pos+jn,0);
	iv++;
      }
    }
    mat->ContributeInterface(x,soll,solr,dsoll,dsolr,weight,fNormal,phixl,phixr,dphixl,dphixr,ek.fMat,ef.fMat);
  }
  delete intrule;
}


void TPZInterfaceElement::CalcResidualStandard(TPZElementMatrix &ef){


  TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial);
#ifndef NODEBUG
  if(!mat || !strcmp("no_name",mat->Name())){
    PZError << "TPZInterfaceElement::CalcResidual interface material null, do nothing\n";
    return;
  }
#endif
  TPZCompElDisc  *left = dynamic_cast<TPZCompElDisc*>( LeftElement() );
  TPZCompElDisc *right = dynamic_cast<TPZCompElDisc*>( RightElement() );
#ifndef NODEBUG
  if(!left->Material() || !right->Material()){
    PZError << "TPZInterfaceElement::CalcResidual null material\n";
    return;
  }
#endif
  //  cout << "TPZInterfaceElement::CalcStiff normal" << fNormal << "left " << left->Reference()->Id() << " right " << right->Reference()->Id() << endl;

  int fConnectIndexL = -1;
  TPZConnect* fConnectL = NULL;
  if (left->NConnects()){
     fConnectIndexL = left->ConnectIndex();
     fConnectL = &( left->Connect(0) );//Discontinuous elements has only 1 connect
  }

  int fConnectIndexR = -1;
  TPZConnect* fConnectR = NULL;
  if (right->NConnects()){
     fConnectIndexR = right->ConnectIndex();
     fConnectR = &( right->Connect(0) );//Discontinuous elements has only 1 connect
  }

  int nshapel = left->NShapeF();
  int nshaper = right->NShapeF();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix &MeshSol = Mesh()->Solution();
  int nstatel = left->Material()->NStateVariables();
  int nstater = right->Material()->NStateVariables();
  int neql = nshapel * nstatel;
  int neqr = nshaper * nstater;
  int dim = Dimension();
  int diml = left->Dimension();
  int dimr = right->Dimension();
  int ncon = 0;
  if(fConnectL) ncon++;
  if(fConnectR) ncon++;


  int neq = neql + neqr;
  //ek.fMat.Redim(neq,neq);
  ef.fMat.Redim(neq,1);
  int ic = 0;
  //ek.fBlock.SetNBlocks(ncon);
  ef.fBlock.SetNBlocks(ncon);
  //ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);
  if(fConnectL) {
      //ek.fBlock.Set(ic,neql);
      ef.fBlock.Set(ic,neql);
    (ef.fConnect)[ic] = fConnectIndexL;
    //(ek.fConnect)[ic] = fConnectIndexL;
      ic++;
  }
  if(fConnectR) {
    //ek.fBlock.Set(ic,neqr);
    ef.fBlock.Set(ic,neqr);
    (ef.fConnect)[ic] = fConnectIndexR;
    //(ek.fConnect)[ic] = fConnectIndexR;
  }
  //ek.fBlock.Resequence();
  ef.fBlock.Resequence();

  TPZFNMatrix<100> phixl(nshapel,1),dphixl(diml,nshapel);
  TPZFNMatrix<100> phixr(nshaper,1),dphixr(dimr,nshaper);
  TPZFNMatrix<9> axes(3,3);
  TPZFNMatrix<9> jacobian(dim,dim);
  TPZFNMatrix<9> jacinv(dim,dim);
  TPZManVector<REAL,3> x(3);
  TPZManVector<REAL,3> intpoint(dim);
  REAL detjac,weight;
  TPZManVector<REAL,5> soll(nstatel),solr(nstater);
  TPZFNMatrix<15> dsoll(diml,nstatel),dsolr(dimr,nstater);
  int pl = left->Degree();
  int pr = right->Degree();
  int p = (pl > pr) ? pl : pr;
  TPZGeoEl *ref = Reference();
  int face = ref->NSides()-1;
  TPZIntPoints *intrule = ref->CreateSideIntegrationRule(face,2*p);//integra u(n)*fi
  int npoints = intrule->NPoints();
  int ip;
  //  TPZManVector<REAL,3> point(3);
//  TPZBndCond *bcleft = 0,*bcright=0;

  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    weight *= fabs(detjac);
    ref->X(intpoint, x);
    //solu� da itera� anterior
    if(fConnectL){
      left->Shape(x,phixl,dphixl);
      soll.Fill(0.);
      dsoll.Zero();
      TPZConnect *df = fConnectL;
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      int pos = block.Position(dfseq);
      int iv = 0,d;
      for(int jn=0; jn<dfvar; jn++) {
	soll[iv%nstatel] += phixl(iv/nstatel,0)*MeshSol(pos+jn,0);
	for(d=0; d<diml; d++)
	  dsoll(d,iv%nstatel) += dphixl(d,iv/nstatel)*MeshSol(pos+jn,0);
	iv++;
      }
    }
    //solu� da itera� anterior
    if(fConnectR){
      right->Shape(x,phixr,dphixr);
      solr.Fill(0.);
      dsolr.Zero();
      TPZConnect *df = fConnectR;
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      int pos = block.Position(dfseq);
      int iv = 0,d;
      for(int jn=0; jn<dfvar; jn++) {
	solr[iv%nstater] += phixr(iv/nstater,0)*MeshSol(pos+jn,0);
	for(d=0; d<dimr; d++)
	  dsolr(d,iv%nstater) += dphixr(d,iv/nstater)*MeshSol(pos+jn,0);
	iv++;
      }
    }

    mat->ContributeInterface(x,soll,solr,dsoll,dsolr,weight,fNormal,phixl,phixr,dphixl,dphixr,/*ek.fMat,*/ef.fMat);
  }

  delete intrule;
}


void TPZInterfaceElement::CalcStiffPenalty(TPZElementMatrix &ek, TPZElementMatrix &ef){

  TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial);
  if(!mat || !strcmp("no_name",mat->Name())){
    PZError << "TPZInterfaceElement::CalcStiff interface material null, do nothing\n";
    return;
  }
  TPZCompElDisc  *left = dynamic_cast<TPZCompElDisc*>( LeftElement() );
  TPZCompElDisc *right = dynamic_cast<TPZCompElDisc*>( RightElement() );
  if(!left->Material() || !right->Material()){
    PZError << "TPZInterfaceElement::CalcStiff null material\n";
    return;
  }

  int nshapel = left->NShapeF();
  int nshaper = right->NShapeF();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix &MeshSol = Mesh()->Solution();
  int nstatel = left->Material()->NStateVariables();
  int nstater = right->Material()->NStateVariables();
  int neql = nshapel * nstatel;
  int neqr = nshaper * nstater;
  int dim = Dimension();
  int diml = left->Dimension();
  int dimr = right->Dimension();
  int ncon = NConnects(),i;


  int neq = neql + neqr;
  ek.fMat.Redim(neq,neq);
  ef.fMat.Redim(neq,1);
  if(ncon){//no m�imo ncon = 1
    int ic = 0;
    ek.fBlock.SetNBlocks(ncon);
    ef.fBlock.SetNBlocks(ncon);
    if(left->NConnects()) {
      ek.fBlock.Set(ic,left->NShapeF()*nstatel);
      ef.fBlock.Set(ic,left->NShapeF()*nstatel);
      ic++;
    }
    if(right->NConnects()) {
      ek.fBlock.Set(ic,right->NShapeF()*nstater);
      ef.fBlock.Set(ic,right->NShapeF()*nstater);
    }
    ek.fBlock.Resequence();
    ef.fBlock.Resequence();
  }
  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);
  for(i=0;i<ncon;i++){//pelo menos 2: 1 para cada elemento right e left
    (ef.fConnect)[i] = ConnectIndex(i);
    (ek.fConnect)[i] = ConnectIndex(i);
  }

  TPZFMatrix phixl(neql,1,0.),dphixl(diml,neql);
  TPZFMatrix phixr(neqr,1,0.),dphixr(dimr,neqr);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL detjac,weight;
  TPZVec<REAL> soll(nstatel,0.),solr(nstater,0.);
  TPZFMatrix dsoll(diml,nstatel,0.),dsolr(dimr,nstater,0.);
  int pl = left->Degree();
  int pr = right->Degree();
  int p = (pl > pr) ? pl : pr;
  TPZGeoEl *ref = Reference();
  int face = ref->NSides()-1;
  TPZIntPoints *intrule = ref->CreateSideIntegrationRule(face,2*p);//integra u(n)*fi
  int npoints = intrule->NPoints();
  int ip;
  TPZVec<REAL> point(3,0.),normal(3,0.);
//  TPZBndCond *bcleft = 0,*bcright=0;

  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    weight *= fabs(detjac);
    ref->X(intpoint, x);
    left->Shape(x,phixl,dphixl);
    right->Shape(x,phixr,dphixr);
    int nderivr = dphixr.Rows();
    int nderivl = dphixl.Rows();
    dsoll.Redim(nderivl,nstatel);
    dsolr.Redim(nderivr,nstater);
    //solu� da itera� anterior
    soll.Fill(0.);
    dsoll.Zero();
    if(left->NConnects()){
      TPZConnect *df = &left->Connect(0);
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      int pos = block.Position(dfseq);
      int iv = 0,d;
      for(int jn=0; jn<dfvar; jn++) {
	soll[iv%nstatel] += phixl(iv/nstatel,0)*MeshSol(pos+jn,0);
	for(d=0; d<nderivl; d++)
	  dsoll(d,iv%nstatel) += dphixl(d,iv/nstatel)*MeshSol(pos+jn,0);
	iv++;
      }
    } 

    //solu� da itera� anterior
    solr.Fill(0.);
    dsolr.Zero();
    if(right->NConnects()){
      TPZConnect *df = &right->Connect(0);
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      int pos = block.Position(dfseq);
      int iv = 0,d;
      for(int jn=0; jn<dfvar; jn++) {
	solr[iv%nstater] += phixr(iv/nstater,0)*MeshSol(pos+jn,0);
	for(d=0; d<nderivr; d++)
	  dsolr(d,iv%nstater) += dphixr(d,iv/nstater)*MeshSol(pos+jn,0);
	iv++;
      }
    } 

    REAL faceSize;
    if (this->Reference()->Dimension() == 0){ //it means point
       //2*(a+b)/2
       faceSize = left->InnerRadius() + right->InnerRadius();//we cannot use right->Reference()->ElementRadius because of right may be an agglomerate
    }
    else{
       faceSize = 2. * this->Reference()->ElementRadius();
    }

    mat->ContributeInterface(x,soll,solr,dsoll,dsolr,weight,fNormal,phixl,phixr,dphixl,dphixr,ek.fMat,ef.fMat, left->Degree(), right->Degree(), faceSize);

  }
delete intrule;
}





void TPZInterfaceElement::GetTransformsLeftAndRight(TPZTransform &tl,TPZTransform &tr){

   TPZCompEl* fLeftEl  = fLeftElSide.Element();
   TPZCompEl* fRightEl = fRightElSide.Element();

  TPZGeoEl *refl = fLeftEl->Reference();
  TPZGeoEl *refr = fRightEl->Reference();
  TPZGeoElSide leftside;
  TPZGeoElSide rightside;
  TPZGeoEl *ref = Reference();
  int i,face = ref->NSides()-1;

  TPZStack<TPZCompElSide> list;
  list.Resize(0);
  TPZCompElSide thisface(this,face);
  //a interface sempre pertence a face menor
  thisface.EqualLevelElementList(list,0,0);
  //procurando left ou right vizinho da atual interface
  int cap = list.NElements();
  for(i=0;i<cap;i++){
    TPZGeoElSide geoside = list[i].Reference();
    if(geoside.Element() == refl) leftside = geoside;
    if(geoside.Element() == refr) rightside = geoside;
  }
  //o elemento interface n� tem pai
  TPZGeoElSide thisgeoside(ref,face);
  TPZCompElSide lower;
  if(!rightside.Exists() && leftside.Exists()){
    lower =  leftside.Reference().LowerLevelElementList(0);
    if(lower.Exists()){
      rightside = lower.Reference();
      TPZTransform acum(leftside.Dimension());
      thisgeoside.SideTransform3(leftside,acum);
      tl = acum;
      leftside.SideTransform3(rightside,acum);//acumula?
      tr = acum;
      return;
    }
  } else
  if(!leftside.Exists() && rightside.Exists()){
    lower = rightside.Reference().LowerLevelElementList(0);
    if(lower.Exists()){
      leftside = lower.Reference();
      TPZTransform acum(rightside.Dimension());
      thisgeoside.SideTransform3(rightside,acum);
      tr = acum;
      rightside.SideTransform3(leftside,acum);//acumula?
      tl = acum;
      return;
    }
  }
  if(!leftside.Exists() || !rightside.Exists()){
    PZError << "TPZInterfaceElement::CalcStiff element does not exist, is impossible\n";
    int dummy;
    cin >> dummy;
  }
  //aqui left e right s� vizinhos
  TPZTransform t2l(leftside.Dimension()),t2r(rightside.Dimension());
  thisgeoside.SideTransform3(leftside,t2l);
  thisgeoside.SideTransform3(rightside,t2r);
  tl = t2l;
  tr = t2r;
}

int TPZInterfaceElement::NConnects() {
   return this->NLeftConnects() + this->NRightConnects();
}

int TPZInterfaceElement::NLeftConnects(){
   TPZCompEl * LeftEl  = fLeftElSide.Element();
   return LeftEl->NConnects();
}

int TPZInterfaceElement::NRightConnects(){
   TPZCompEl * RightEl = fRightElSide.Element();
   return RightEl->NConnects();
}

int TPZInterfaceElement::ConnectIndex(int i) {

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
}

void TPZInterfaceElement::Print(ostream &out){

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

void VetorialProd(TPZVec<REAL> &ivet,TPZVec<REAL> &jvet,TPZVec<REAL> &kvet);

void TPZInterfaceElement::NormalToFace(TPZVec<REAL> &normal /*,int leftside*/){

   TPZCompEl * fLeftEl = this->LeftElement();
   TPZCompEl * fRightEl = this->RightElement();

  //  int dim = Reference()->Dimension();
  TPZGeoEl *ref = Reference();
  int face = ref->NSides()-1;
  //face: lado do elemento bidimensional ou aresta 
  //do unidimensional ou canto do ponto
  normal.Resize(3,0.);
  normal.Fill(0.);
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
     normal[0] = 1.0;// a normal sempre apontar�na dire� positiva do eixo
     normal[1] = 0.;
     normal[2] = 0.;
   break;
  case 1:
    ref->Jacobian(param,jacobian,axes,detjac,jacinv);
    for(i=0;i<3;i++) rib[i] = axes(0,i);//dire� da aresta
    VetorialProd(rib,vec,result);
    VetorialProd(result,rib,normal);
    //normalizando a normal
    normalize = 0.;
    for(i=0;i<3;i++) normalize += normal[i]*normal[i];
    if(normalize == 0.0)
      PZError << "TPZInterfaceElement::NormalToFace null normal vetor\n";
    normalize = sqrt(normalize);
    for(i=0;i<3;i++) normal[i] = normal[i]/normalize;
    
    break;
  case 2: 
    ref->CenterPoint(face,param);//ponto da face
    ref->Jacobian(param,jacobian,axes,detjac,jacinv);
    for(i=0;i<3;i++) normal[i] = axes(2,i);
    break;
  default:
    PZError << "TPZInterfaceElement::NormalToFace in case that not treated\n";
    normal.Resize(0);
    return;
  }

  //to guarantee the normal points from left to right neighbours:
  normalize = 0.;
  for(i=0; i<3; i++) normalize += normal[i]*vec[i];
  if(normalize < 0.) {
    for(i=0; i<3; i++) normal[i] = -normal[i];
  }
}

void VetorialProd(TPZVec<REAL> &ivet,TPZVec<REAL> &jvet,TPZVec<REAL> &kvet){

  kvet.Resize(3);
  kvet[0] =  ivet[1]*jvet[2] - ivet[2]*jvet[1];
  kvet[1] = -ivet[0]*jvet[2] + ivet[2]*jvet[0];
  kvet[2] =  ivet[0]*jvet[1] - ivet[1]*jvet[0];
}

void TPZInterfaceElement::Normal(TPZVec<REAL> &normal) {

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
  /**
  Save the element data to a stream
  */
void TPZInterfaceElement::Write(TPZStream &buf, int withclassid)
{
  TPZCompEl::Write(buf,withclassid);
  int leftelindex = fLeftElSide.Element()->Index();
  int rightelindex = fRightElSide.Element()->Index();
  if ( (this->Index() > leftelindex) || (this->Index() > rightelindex) ){
     PZError << __PRETTY_FUNCTION__ << endl
	     << "Indices of neighbours are less than interface index:" << endl
	     << "Left: " << leftelindex << ", Right: " << rightelindex << ", this: " << this->Index() << endl;
  }
  int matid = fMaterial->Id();

  int leftside = fLeftElSide.Side();
  int rightside = fRightElSide.Side();

  buf.Write(&leftelindex,1);
  buf.Write(&leftside,1);
  buf.Write(&rightelindex,1);
  buf.Write(&rightside,1);
  buf.Write(&matid,1);
  WriteObjects(buf,fNormal);
}
  
  /**
  Read the element data from a stream
  */
void TPZInterfaceElement::Read(TPZStream &buf, void *context)
{
  TPZCompEl::Read(buf,context);
  int leftelindex;
  int rightelindex;
  int leftside, rightside;
  int matid;
  buf.Read(&leftelindex,1);
  buf.Read(&leftside,1);
  buf.Read(&rightelindex,1);
  buf.Read(&rightside,1);
  buf.Read(&matid,1);
  fMaterial = Mesh()->FindMaterial(matid);
  if(!fMaterial) PZError << "TPZInterfaceElement::Read - null material\n";

  this->fLeftElSide.SetElement ( Mesh()->ElementVec()[leftelindex]  );
  this->fRightElSide.SetElement( Mesh()->ElementVec()[rightelindex] );
  this->fLeftElSide.SetSide( leftside );
  this->fRightElSide.SetSide( rightside );

  ReadObjects(buf,fNormal);
}

void TPZInterfaceElement::CalcStiffContDisc(TPZElementMatrix &ek, TPZElementMatrix &ef){


   TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial);
#ifndef NODEBUG
   if(!mat || !strcmp("no_name",mat->Name())){
      PZError << "TPZInterfaceElement::CalcStiff interface material null, do nothing\n";
      return;
   }
#endif

   TPZCompElDisc  *discL = NULL;
   TPZCompElDisc  *discR = NULL;
   TPZInterpolatedElement* intelL = NULL;
   TPZInterpolatedElement* intelR = NULL;
   const int leftside = fLeftElSide.Side();
   const int rightside = fRightElSide.Side();
   TPZCompEl * left = this->LeftElement();
   TPZCompEl * right = this->RightElement();
  
   if(leftside == -1){ //i.e. TPZCompElDisc
      discL = dynamic_cast<TPZCompElDisc*>( left );  
   } 
   else{//i.e. TPZInterpolatedElement
      intelL = dynamic_cast<TPZInterpolatedElement*>( left );
   }

   if (rightside == -1){//i.e. TPZCompElDisc
      discR = dynamic_cast<TPZCompElDisc*>( right );
   }
   else{//i.e. TPZInterpolatedElement
      intelR = dynamic_cast<TPZInterpolatedElement*>( right );
   }

#ifndef NODEBUG
   if(!left->Material() || !right->Material()){
      PZError << "TPZInterfaceElement::CalcStiff null material\n";
      return;
   }
#endif

   TPZManVector<TPZConnect*> ConnectL, ConnectR;
   TPZManVector<int> ConnectIndexL, ConnectIndexR;

   this->GetConnects( fLeftElSide,  ConnectL, ConnectIndexL );
   this->GetConnects( fRightElSide, ConnectR, ConnectIndexR );

   int nshapel = -1;
   if (discL)  nshapel = discL ->NShapeF();
   if (intelL) nshapel = intelL->NShapeF();
   int nshaper = -1;
   if (discR)  nshaper = discR->NShapeF();
   if (intelR) nshaper = intelR->NShapeF();

   TPZBlock &block = Mesh()->Block();
   TPZFMatrix &MeshSol = Mesh()->Solution();
   const int nstatel = left->Material()->NStateVariables();
   const int nstater = right->Material()->NStateVariables();
   const int neql = nshapel * nstatel;
   const int neqr = nshaper * nstater;
   const int dim = this->Dimension();
   const int diml = left->Dimension();
   const int dimr = right->Dimension();
   const int ncon = ConnectL.NElements() + ConnectR.NElements();

   const int neq = neql + neqr;
   ek.fMat.Redim(neq,neq);
   ef.fMat.Redim(neq,1);
   ek.fBlock.SetNBlocks(ncon);
   ef.fBlock.SetNBlocks(ncon);
   ek.fConnect.Resize(ncon);
   ef.fConnect.Resize(ncon);

   int ic = 0;
   {
      int n = ConnectL.NElements();
      for(int i = 0; i < n; i++) {
	 int nshape = 0;
	 if(discL) nshape = discL->NShapeF();
	 if(intelL) nshape = intelL->NConnectShapeF(i);
	 int con_neq = nstatel * nshape;
	 ek.fBlock.Set(ic,con_neq );
	 ef.fBlock.Set(ic,con_neq);
	 (ef.fConnect)[ic] = ConnectIndexL[i];
	 (ek.fConnect)[ic] = ConnectIndexL[i];
	 ic++;
      }
   }

   {
      int n = ConnectR.NElements();
      for(int i = 0; i < n; i++) {
	 int nshape = 0;
	 if(discR) nshape = discR->NShapeF();
	 if(intelR) nshape = intelR->NConnectShapeF(i);
	 int con_neq = nstater * nshape;
	 ek.fBlock.Set(ic,con_neq );
	 ef.fBlock.Set(ic,con_neq);
	 (ef.fConnect)[ic] = ConnectIndexR[i];
	 (ek.fConnect)[ic] = ConnectIndexR[i];
	 ic++;
      }
   }

   ek.fBlock.Resequence();
   ef.fBlock.Resequence();

   TPZFNMatrix<100> phixl(nshapel,1),dphixl(diml,nshapel);
   TPZFNMatrix<100> phixr(nshaper,1),dphixr(dimr,nshaper);
   TPZFNMatrix<9> axes(3,3);
   TPZFNMatrix<9> jacobian(dim,dim);
   TPZFNMatrix<9> jacinv(dim,dim);
   TPZManVector<REAL,3> x(3);
   TPZManVector<REAL,3> intpoint(dim), LeftIntPoint(diml), RightIntPoint(dimr);
   REAL detjac,weight;
   TPZManVector<REAL,5> soll(nstatel),solr(nstater);
   TPZFNMatrix<15> dsoll(diml,nstatel),dsolr(dimr,nstater);

   //LOOKING FOR MAX INTERPOLATION ORDER
   int pl, pr;

   //Left element
   if (discL)  pl = discL->Degree();
   if (intelL){
      int is, nsides = intelL->Reference()->NSides();
      int order = 0;
      pl = 0;
      for(is = 0; is < nsides; is++){
	 order = intelL->SideOrder( is );
	 if (order > pl) pl = order;
      }
   }

   //Right element
   if (discR)  pr = discR->Degree();
   if (intelR) {
      int is, nsides = intelR->Reference()->NSides();
      int order = 0;
      pr = 0;
      for(is = 0; is < nsides; is++){
	 order = intelR->SideOrder( is );
	 if (order > pr) pr = order;
      }
   }
   
   //Max interpolation order
   const int p = (pl > pr) ? pl : pr;

   TPZGeoEl *ref = Reference();
   TPZIntPoints *intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, 2*(p+1) );
   if(fMaterial->HasForcingFunction()) {
      TPZManVector<int> order(3);
      intrule->GetOrder(order);
      int maxorder = intrule->GetMaxOrder();
      order.Fill(maxorder);
      intrule->SetOrder(order);
   }
   const int npoints = intrule->NPoints();


   //integration points in left and right elements: making transformations to interpolated elements
   TPZTransform TransfLeft(dim/*l*/), TransfRight(dim/*r*/);

   if (intelL){

      TPZGeoElSide thisgeoside(this->Reference(), this->Reference()->NSides()-1);
      TPZGeoElSide leftgeoside(left->Reference(), leftside);
      thisgeoside.SideTransform3( leftgeoside, TransfLeft );

      TPZGeoElSide highdim(left->Reference(), left->Reference()->NSides()-1);
      TransfLeft = leftgeoside.SideToSideTransform(highdim).Multiply(TransfLeft);
   }


   if(intelR){

      TPZGeoElSide thisgeoside(this->Reference(), this->Reference()->NSides()-1);
      TPZGeoElSide rightgeoside(right->Reference(), rightside);
      thisgeoside.SideTransform3( rightgeoside,TransfRight );

      TPZGeoElSide highdim( right->Reference(), right->Reference()->NSides()-1 );
      TransfRight = rightgeoside.SideToSideTransform(highdim).Multiply(TransfRight);
   }
   
   
   //LOOP OVER INTEGRATION POINTS
   for(int ip = 0; ip < npoints; ip++){

      intrule->Point(ip,intpoint,weight);
      ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
      weight *= fabs(detjac);
      ref->X(intpoint, x);

      if (intelL) TransfLeft.Apply( intpoint, LeftIntPoint );
      if (intelR) TransfRight.Apply( intpoint, RightIntPoint );

#ifdef DEBUG
      {
	 const REAL tol = 1.e-10;

	 if (intelL){
	    TPZManVector<REAL> FaceXPoint(3), LeftXPoint(3);
	    this->Reference()->X( intpoint, FaceXPoint);
	    intelL->Reference()->X( LeftIntPoint, LeftXPoint);
	    int i, n = FaceXPoint.NElements();
	    if (n != LeftXPoint.NElements() ){
	       PZError << __PRETTY_FUNCTION__ << endl
		       << "Face X point and LeftElement X point have not same dimension." << endl;
	    }
	    REAL erro = 0.;
	    for(i = 0; i < n; i++){
	       erro += (LeftXPoint[i] - FaceXPoint[i])*(LeftXPoint[i] - FaceXPoint[i]);
	    }
	    erro = sqrt(erro);
	    if (erro > tol) PZError << __PRETTY_FUNCTION__ << endl 
				    << "Face X point and LeftElement X point are not same." << endl;
	 }

	 if (intelR){
	    TPZManVector<REAL> FaceXPoint(3), RightXPoint(3);
	    this->Reference()->X( intpoint, FaceXPoint);
	    intelR->Reference()->X( LeftIntPoint, RightXPoint);
	    int i, n = FaceXPoint.NElements();
	    if (n != RightXPoint.NElements() ){
	       PZError << __PRETTY_FUNCTION__ << endl
		       << "Face X point and RightElement X point have not same dimension." << endl;
	    }
	    REAL erro = 0.;
	    for(i = 0; i < n; i++){
	       erro += (RightXPoint[i] - FaceXPoint[i])*(RightXPoint[i] - FaceXPoint[i]);
	    }
	    erro = sqrt(erro);
	    if (erro > tol) PZError << __PRETTY_FUNCTION__ << endl 
				    << "Face X point and RightElement X point are not same." << endl;
	 }
	 
      }
#endif


      //COMPUTING SHAPE FUNCTIONS - LEFT
      if (discL) discL->Shape(x,phixl,dphixl);

      if (intelL){  
	 this->ComputeShape(intelL, phixl, dphixl, LeftIntPoint );
      }//intelL


      //COMPUTING SHAPE FUNCTIONS - RIGHT
      if (discR) discR->Shape(x,phixr,dphixr);

      if (intelR){
	 this->ComputeShape(intelR, phixr, dphixr, RightIntPoint );
      }//intelR




      //solu� da iteracao anterior - Left
      int n = ConnectL.NElements();
      for( int i = 0; i < n; i++ ) {
	 soll.Fill(0.);
	 dsoll.Zero();
	 TPZConnect *df = ConnectL[i];
	 int dfseq = df->SequenceNumber();
	 int dfvar = block.Size(dfseq);
	 int pos = block.Position(dfseq);
	 int iv = 0,d;
	 for(int jn=0; jn<dfvar; jn++) {
	    soll[iv%nstatel] += phixl(iv/nstatel,0)*MeshSol(pos+jn,0);
	    for(d=0; d<diml; d++)
	       dsoll(d,iv%nstatel) += dphixl(d,iv/nstatel)*MeshSol(pos+jn,0);
	    iv++;
	 }
      }

      //solu� da iteracao anterior - Right
      n = ConnectR.NElements();
      for( int i = 0; i < n; i++ ){
	 solr.Fill(0.);
	 dsolr.Zero();
	 TPZConnect *df = ConnectR[i];
	 int dfseq = df->SequenceNumber();
	 int dfvar = block.Size(dfseq);
	 int pos = block.Position(dfseq);
	 int iv = 0,d;
	 for(int jn=0; jn<dfvar; jn++) {
	    solr[iv%nstater] += phixr(iv/nstater,0)*MeshSol(pos+jn,0);
	    for(d=0; d<dimr; d++)
	       dsolr(d,iv%nstater) += dphixr(d,iv/nstater)*MeshSol(pos+jn,0);
	    iv++;
	 }
      }

      //CONTRIBUTING TO STIFF MATRIX AND LOAD VECTOR
      mat->ContributeInterface(x,soll,solr,dsoll,dsolr,weight,fNormal,phixl,phixr,dphixl,dphixr,ek.fMat,ef.fMat);

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

void TPZInterfaceElement::ComputeShape(TPZInterpolatedElement* intel, TPZFMatrix &phix, TPZFMatrix &dphix, TPZVec<REAL> &IntPoint ) {

   //COMPUTING SHAPE FUNCTIONS TO INTERPOLATED NEIGHBOUR
   
   {//teste
     TPZVec< TPZVec<REAL> > copia(3);
     TPZVec<REAL> *local;
     local = &copia[0];
     local->Resize(2);
     (*local)[0] = IntPoint[0];
     (*local)[1] = 0.;
     
     local = &copia[1];
     local->Resize(2);
     (*local)[0] = -1. * IntPoint[0];
     (*local)[1] = 0.;
     
     local = &copia[2];
     local->Resize(2);
     (*local)[0] = 0.;
     (*local)[1] = 0.;        
     
     for(int ivec = 0; ivec < copia.NElements(); ivec++){
       
       TPZVec<REAL> &Ponto = copia[ivec];
        
       TPZFNMatrix<100> dphi( dphix.Rows(), dphix.Cols() );
      
      intel->Shape( Ponto, phix, dphi);

      const int dim    = intel->Dimension();
      const int nshape = intel->NShapeF();
      TPZFNMatrix<9> axes(3,3);
      TPZFNMatrix<9> jacobian(dim,dim), jacinv(dim,dim);
      REAL detjac;

      intel->Reference()->Jacobian( Ponto, jacobian, axes, detjac, jacinv );
	    
      int ieq;
      switch(dim) {
	 case 0:
	    break;
	 case 1:
	    dphix = dphi;
	    dphix *= (1./detjac);
	    break;
	 case 2:
	    for(ieq = 0; ieq < nshape; ieq++) {
	       dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
	       dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
	    }
	    break;
	 case 3:
	    for(ieq = 0; ieq < nshape; ieq++) {
	       dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
	       dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
	       dphix(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
	    }
	    break;
	 default:
	    PZError << "TPZInterface please implement the " << dim << "d Jacobian and inverse\n";
	    PZError.flush();
      }   
      
      }
   
   }
   
   

   if (intel){  
      
      TPZFNMatrix<100> dphi( dphix.Rows(), dphix.Cols() );
      
      intel->Shape( IntPoint, phix, dphi);

      const int dim    = intel->Dimension();
      const int nshape = intel->NShapeF();
      TPZFNMatrix<9> axes(3,3);
      TPZFNMatrix<9> jacobian(dim,dim), jacinv(dim,dim);
      REAL detjac;

      intel->Reference()->Jacobian( IntPoint, jacobian, axes, detjac, jacinv );
	    
      int ieq;
      switch(dim) {
	 case 0:
	    break;
	 case 1:
	    dphix = dphi;
	    dphix *= (1./detjac);
	    break;
	 case 2:
	    for(ieq = 0; ieq < nshape; ieq++) {
	       dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
	       dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
	    }
	    break;
	 case 3:
	    for(ieq = 0; ieq < nshape; ieq++) {
	       dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
	       dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
	       dphix(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
	    }
	    break;
	 default:
	    PZError << "TPZInterface please implement the " << dim << "d Jacobian and inverse\n";
	    PZError.flush();
      }

#warning Falta aplicar AXES
   }//intelL

}


