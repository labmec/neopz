//$Id: TPZInterfaceEl.cpp,v 1.34 2004-04-06 14:55:43 erick Exp $

#include "pzelmat.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzquad.h"
#include "pzmaterial.h"
#include "TPZConservationLaw.h"
#include "pzconslaw.h"
#include "pzbndcond.h"

int TPZInterfaceElement::gCalcStiff = 1;

/**
 * Para CloneInterface.
 * A normal é clonada, e não recalculada.
 */
TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElDisc *left,TPZCompElDisc *right, TPZVec<REAL> normal)
  : TPZCompEl(mesh,index), fNormal(3,0.), fConnectL(0), fConnectR(0), fConnectIndexL(-1), fConnectIndexR(-1) {

  fReference = geo;
  geo->SetReference(this);
  int materialid = fReference->MaterialId();
  //  if(materialid<0){
  //    PZError << "TPZInterfaceElement::TPZInterfaceElement the interface element is not a BC condition\n";
  //  }
  //poderia eliminar esta variável e carrega-la do elemento de volume associado
  fMaterial = mesh.FindMaterial(materialid);
  fLeftEl = left;
  fRightEl = right;
  int ic = 0;
  int icon;
  for(icon = 0; icon < left->NConnects(); icon++){
    left ->Connect(icon).IncrementElConnected();
    fConnectL = &left->Connect(icon);
    fConnectIndexL = left->ConnectIndex(icon);
    ic++;
  }

  for(icon = 0; icon < right->NConnects(); icon++){
    right ->Connect(icon).IncrementElConnected();
    fConnectR = &right->Connect(icon);
    fConnectIndexR = right->ConnectIndex(icon);
    ic++;
  }

  for (int i = 0; i < fNormal.NElements(); i++)
    fNormal[i] = normal[i];
  //NormalToFace(fNormal/*,leftside*/);
}




//construtor para o elemento descontinuo
TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,
					 TPZCompElDisc *left,TPZCompElDisc *right/*,int leftside*/) 
  : TPZCompEl(mesh,index), fNormal(3,0.) , fConnectL(0), fConnectR(0), fConnectIndexL(-1), fConnectIndexR(-1){

  fReference = geo;
  geo->SetReference(this);
  int materialid = fReference->MaterialId();
  //  if(materialid<0){
  //    PZError << "TPZInterfaceElement::TPZInterfaceElement the interface element is not a BC condition\n";
  //  }
  //poderia eliminar esta variável e carrega-la do elemento de volume associado
  fMaterial = mesh.FindMaterial(materialid);
  fLeftEl = left;
  fRightEl = right;
  int icon;
  int ic=0;
  for(icon = 0; icon < left->NConnects(); icon++){
    left ->Connect(icon).IncrementElConnected();
    fConnectL = &left->Connect(icon);
    fConnectIndexL = left->ConnectIndex(icon);
    ic++;
  }

  for(icon = 0; icon < right->NConnects(); icon++){
    right ->Connect(icon).IncrementElConnected();
    fConnectR = &right->Connect(icon);
    fConnectIndexR = right->ConnectIndex(icon);
    ic++;
  }

  NormalToFace(fNormal/*,leftside*/);
}

TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh, const TPZInterfaceElement &copy)
  : TPZCompEl(mesh,copy), fNormal(copy.fNormal), fConnectL(0), fConnectR(0), fConnectIndexL(-1), fConnectIndexR(-1) {
  fLeftEl = dynamic_cast<TPZCompElDisc *>(mesh.ElementVec()[copy.fLeftEl->Index()]);
  fRightEl = dynamic_cast<TPZCompElDisc *>(mesh.ElementVec()[copy.fRightEl->Index()]);

#ifdef DEBUG
  if(!fLeftEl || ! fRightEl) {
    cout << "Something wrong with clone of interface element\n";
    exit(-1);
  }
  if(fLeftEl->Mesh() != &mesh || fRightEl->Mesh() != &mesh) {
    cout << "The discontinuous elements should be cloned before the interface elements\n";
    exit(-1);
  }
#endif

  fReference = copy.fReference;
  TPZMaterial *mat = copy.Material();
  if(mat) {
    int materialid = mat->Id();
    fMaterial = mesh.FindMaterial(materialid);
  } else {
    fMaterial = 0;
  }
  //  fMaterial = copy.fMaterial;

  int ic=0;
  int icon;
  for(icon = 0; icon < fLeftEl->NConnects(); icon++){
    fLeftEl ->Connect(icon).IncrementElConnected();
    fConnectL = &fLeftEl->Connect(icon);
    fConnectIndexR = fLeftEl->ConnectIndex(icon);
    ic++;
  }

  for(icon = 0; icon < fRightEl->NConnects(); icon++){
    fRightEl ->Connect(icon).IncrementElConnected();
    fConnectR = &fRightEl->Connect(icon);
    fConnectIndexR = fRightEl->ConnectIndex(icon);
    ic++;
  }
}

TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,const TPZInterfaceElement &copy,int &index) 
  : TPZCompEl(mesh,index), fNormal(copy.fNormal), fConnectL(0), fConnectR(0), fConnectIndexL(-1), fConnectIndexR(-1) {

  //ambos elementos esquerdo e direito já foram clonados e moram na malha aglomerada
  //o geometrico da malha fina aponta para o computacional da malha aglomerada
  TPZCompEl *left = copy.fLeftEl->Reference()->Reference();
  TPZCompEl *right = copy.fRightEl->Reference()->Reference();
  fLeftEl = dynamic_cast<TPZCompElDisc *>(left);
  fRightEl = dynamic_cast<TPZCompElDisc *>(right);
  
#ifdef DEBUG
  if(!fLeftEl || ! fRightEl) {
    cout << "TPZInterfaceElement::TPZInterfaceElement Something wrong with clone of interface element\n";
    exit(-1);
  } 
  if(fLeftEl->Mesh() != &mesh || fRightEl->Mesh() != &mesh) {
    cout << "TPZInterfaceElement::TPZInterfaceElement The discontinuous elements should be cloned "
	 << "before the interface elements\n";
    exit(-1);
  }
#endif 

  fReference = copy.fReference;
  TPZMaterial *mat = copy.Material();
  if(mat) {
    int materialid = mat->Id();
    fMaterial = mesh.FindMaterial(materialid);
  } else {
    fMaterial = 0;
  }
  //  fMaterial = copy.fMaterial;
  int ic=0;
  int icon;
  for(icon = 0; icon < fLeftEl->NConnects(); icon++){
    fLeftEl ->Connect(icon).IncrementElConnected();
    fConnectL = &fLeftEl->Connect(icon);
    fConnectIndexL = fLeftEl->ConnectIndex(icon);
    ic++;
  }

  for(icon = 0; icon < fRightEl->NConnects(); icon++){
    fRightEl ->Connect(icon).IncrementElConnected();
    fConnectR = &fRightEl->Connect(icon);
    fConnectIndexR = fRightEl->ConnectIndex(icon);
    ic++;
  }

}

TPZCompEl * TPZInterfaceElement::CloneInterface(TPZCompMesh &aggmesh,int &index, TPZCompElDisc *left, TPZCompElDisc *right) const {
  //  return new TPZInterfaceElement(aggmesh, *this, index);
  //<!>O certo eh esse, mas ate o Phil voltar:  return new TPZInterfaceElement(aggmesh, this->Reference(), index, left, right);
  return new TPZInterfaceElement(aggmesh, this->Reference(), index, left, right, this->fNormal);
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
  TPZCompElDisc *left = LeftElement();
  TPZCompElDisc *right = RightElement();
#ifndef NODEBUG
  if(!left->Material() || !right->Material()){
    PZError << "TPZInterfaceElement::CalcStiff null material\n";
    return;
  }
#endif
  //  cout << "TPZInterfaceElement::CalcStiff normal" << fNormal << "left " << left->Reference()->Id() << " right " << right->Reference()->Id() << endl;

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
  int face = fReference->NSides()-1;
  TPZIntPoints *intrule = fReference->CreateSideIntegrationRule(face,2*p);//integra u(n)*fi
  int npoints = intrule->NPoints();
  int ip;
  //  TPZManVector<REAL,3> point(3);
//  TPZBndCond *bcleft = 0,*bcright=0;

  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    weight *= fabs(detjac);
    fReference->X(intpoint, x);
    //solu¢ão da itera¢ão anterior
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
    //solu¢ão da itera¢ão anterior
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
  delete intrule;
}
    mat->ContributeInterface(x,soll,solr,dsoll,dsolr,weight,fNormal,phixl,phixr,dphixl,dphixr,ek.fMat,ef.fMat);
  }


void TPZInterfaceElement::CalcResidualStandard(TPZElementMatrix &ef){


  TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial);
#ifndef NODEBUG
  if(!mat || !strcmp("no_name",mat->Name())){
    PZError << "TPZInterfaceElement::CalcResidual interface material null, do nothing\n";
    return;
  }
#endif
  TPZCompElDisc *left = LeftElement();
  TPZCompElDisc *right = RightElement();
#ifndef NODEBUG
  if(!left->Material() || !right->Material()){
    PZError << "TPZInterfaceElement::CalcResidual null material\n";
    return;
  }
#endif
  //  cout << "TPZInterfaceElement::CalcStiff normal" << fNormal << "left " << left->Reference()->Id() << " right " << right->Reference()->Id() << endl;

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
  int face = fReference->NSides()-1;
  TPZIntPoints *intrule = fReference->CreateSideIntegrationRule(face,2*p);//integra u(n)*fi
  int npoints = intrule->NPoints();
  int ip;
  //  TPZManVector<REAL,3> point(3);
//  TPZBndCond *bcleft = 0,*bcright=0;

  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    weight *= fabs(detjac);
    fReference->X(intpoint, x);
    //solu¢ão da itera¢ão anterior
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
    //solu¢ão da itera¢ão anterior
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
  TPZCompElDisc *left = LeftElement();
  TPZCompElDisc *right = RightElement();
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
  if(ncon){//no máximo ncon = 1
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
  int face = fReference->NSides()-1;
  TPZIntPoints *intrule = fReference->CreateSideIntegrationRule(face,2*p);//integra u(n)*fi
  int npoints = intrule->NPoints();
  int ip;
  TPZVec<REAL> point(3,0.),normal(3,0.);
//  TPZBndCond *bcleft = 0,*bcright=0;

  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    weight *= fabs(detjac);
    fReference->X(intpoint, x);
    left->Shape(x,phixl,dphixl);
    right->Shape(x,phixr,dphixr);
    //solu¢ão da itera¢ão anterior
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
	for(d=0; d<diml; d++)
	  dsoll(d,iv%nstatel) += dphixl(d,iv/nstatel)*MeshSol(pos+jn,0);
	iv++;
      }
    } 

    //solu¢ão da itera¢ão anterior
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
	for(d=0; d<dimr; d++)
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

  TPZGeoEl *refl = fLeftEl->Reference();
  TPZGeoEl *refr = fRightEl->Reference();
  TPZGeoElSide leftside;
  TPZGeoElSide rightside;
  int i,face = fReference->NSides()-1;

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
  //o elemento interface não tem pai
  TPZGeoElSide thisgeoside(fReference,face);
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
  //aqui left e right são vizinhos
  TPZTransform t2l(leftside.Dimension()),t2r(rightside.Dimension());
  thisgeoside.SideTransform3(leftside,t2l);
  thisgeoside.SideTransform3(rightside,t2r);
  tl = t2l;
  tr = t2r;
}

int TPZInterfaceElement::NConnects() {

  if(!fLeftEl && !fRightEl) return 0;
  if(!fLeftEl || !fRightEl) return 1;
  int nl = fLeftEl->NConnects();
  int nr = fRightEl->NConnects();
  return (nl+nr);//that's the right and left element
}

int TPZInterfaceElement::ConnectIndex(int i) {


  if(i<0 || i>1)
    PZError << "TPZInterfaceElement::ConnectIndex wrong argument i, i = " << i << endl;
  int conl = fLeftEl->ConnectIndex(0);
  int conr = fRightEl->ConnectIndex(0);
  //ambos devem ser não nulos: left e right
  if(i == 0 && fLeftEl){
    if(conl > -1) return conl;//de preferência o esquerdo
    if(conr > -1) return conr;
  }
  if(i == 1 && fRightEl){
    if(conr > -1) return conr;//de preferência o direito
    if(conl > -1) return conl;
  }
  return -1;//ambos nulos, esquerdo e direito
}

void TPZInterfaceElement::Print(ostream &out){

  out << "\nInterface element : \n";
  //out << "\tId of the geometric reference : " << fReference->Id() << endl;
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
  out << "\tMaterial id : " << fReference->MaterialId() << endl;
  
  out << "\tNormal a interface : ";
  out << "(" << fNormal[0] << "," << fNormal[1] << "," << fNormal[2] << ")\n";

}

void TPZInterfaceElement::SetConnectIndex(int node, int index) {
  cout << "TPZInterfaceElement::SetConnectIndex should never be called\n";
}

int TPZInterfaceElement::main(TPZCompMesh &cmesh){
  // esta funcão testa o correto desempenho do algoritmo que cria e 
  // deleta elementos de interface numa malha sujeita a refinamento h

  // InterfaceDimension é a dimensão do elemento de interface
  // verifica-se para cada lado de dimensão InterfaceDimension do
  // elemento que existe um elemento interface e que este é único

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
    //caso existem elementos pequenos não deve existir
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
    //caso o vizinho não existe todo bem
    //caso existe não pode ter referência computacional
    return 1;//sem problemas
  }
  
  //neste estagio o lado atual enxerga um elemento vizinho ou
  //está comtido no lado de um elemento maior, portanto deve 
  //ter associado um elemento interface
  TPZGeoElSide geo = comp.Reference();
  if(!geo.Exists()){
    PZError << "TPZInterfaceElement::ExistInterfaces error of data structure\n";
    exit(-1);
  }
  TPZGeoElSide  neigh = geo.Neighbour();
  int exists = 0;
  if(comp.Element()->Type() == EInterface) exists++;//o próprio é interface
  
  while(neigh.Element() && neigh.Element() != geo.Element()){
    TPZCompElSide comp = neigh.Reference();
    neigh = neigh.Neighbour();
    if(!comp.Element()) continue;
    if(comp.Element()->Type() == EInterface) exists++;
  }
  if(exists != 1) return 0;
  return 1;//existe uma única interface
}

int TPZInterfaceElement::FreeInterface(TPZCompMesh &cmesh){

  int iel,nel = cmesh.NElements();
  for(iel=0;iel<nel;iel++){
    TPZCompEl *cel = cmesh.ElementVec()[iel];
    if(!cel) continue;
    if(cel->Type() != EInterface) continue;//interessa só interfaces
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
    //só pode haver 1 ou 2 elementos de volume associados a um el. interface
    if(exists < 1 || exists > 2) return 0;
  }
  return 1;
}

void VetorialProd(TPZVec<REAL> &ivet,TPZVec<REAL> &jvet,TPZVec<REAL> &kvet);

void TPZInterfaceElement::NormalToFace(TPZVec<REAL> &normal /*,int leftside*/){

  //  int dim = fReference->Dimension();
  int face = fReference->NSides()-1;
  //face: lado do elemento bidimensional ou aresta 
  //do unidimensional ou canto do ponto
  normal.Resize(3,0.);
  normal.Fill(0.);
  int faceleft,faceright;

  TPZVec<REAL> param(3),centleft(3),centright(3),point(3,0.),result(3,0.),xint(3),xvolleft(3),xvolright(3),vec(3),rib(3);
  TPZFMatrix jacobian(3,3),jacinv(3,3),axes(3,3);
  REAL detjac,normalize;
  int i;

  faceleft = fLeftEl->Reference()->NSides()-1;//lado interior do elemento esquerdo
  faceright = fRightEl->Reference()->NSides()-1; // lado interior do element direito
  fLeftEl->Reference()->CenterPoint(faceleft,centleft);//ponto centro do elemento de volume
  fRightEl->Reference()->CenterPoint(faceright,centright);
  fLeftEl->Reference()->X(centleft,xvolleft);
  fRightEl->Reference()->X(centright,xvolright);
  for(i=0;i<3;i++) vec[i] = xvolright[i]-xvolleft[i];//não deve ser nulo

  int InterfaceDimension =  fLeftEl->Material()->Dimension() - 1;

  switch(InterfaceDimension){
  case 0:
    normal[0] = 1.0;// a normal sempre apontará na dire¢ão positiva do eixo
    normal[1] = 0.;
    normal[2] = 0.;
    break;
  case 1:
    fReference->Jacobian(param,jacobian,axes,detjac,jacinv);
    for(i=0;i<3;i++) rib[i] = axes(0,i);//dire¢ão da aresta
    VetorialProd(rib,vec,result);
    VetorialProd(result,rib,normal);
    //normalizando a normal
    normalize = 0.;
    for(i=0;i<3;i++) normalize += normal[i]*normal[i];
    if(!normalize)
      PZError << "TPZInterfaceElement::NormalToFace null normal vetor\n";
    normalize = sqrt(normalize);
    for(i=0;i<3;i++) normal[i] = normal[i]/normalize;
    
    break;
  case 2:
    fReference->CenterPoint(face,param);//ponto da face
    fReference->Jacobian(param,jacobian,axes,detjac,jacinv);
    for(i=0;i<3;i++) normal[i] = axes(2,i);
    break;
  default:
    PZError << "TPZInterfaceElement::NormalToFace in case that not treated\n";
    normal.Resize(0);
    return;
  }
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
  TPZVec<REAL> &errors, TPZBlock * /*flux */) {}

