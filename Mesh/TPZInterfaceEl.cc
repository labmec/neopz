// $Id: TPZInterfaceEl.cc,v 1.14 2003-10-01 21:07:15 tiago Exp $

#include "pzelmat.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzquad.h"
#include "pzmaterial.h"
#include "TPZConservationLaw.h"
#include "pzbndcond.h"

//construtor para o elemento descontinuo
TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,
					 TPZCompElDisc *left,TPZCompElDisc *right,int leftside)
  : TPZCompEl(mesh,index), fNormal(3,0.) {

  fReference = geo;
  geo->SetReference(this);
  int materialid = fReference->MaterialId();
  if(materialid<0){
    PZError << "TPZInterfaceElement::TPZInterfaceElement the interface element is not a BC condition\n";
  }
  //poderia eliminar esta variável e carrega-la do elemento de volume associado
  fMaterial = mesh.FindMaterial(materialid);
  if(!fMaterial) PZError << "TPZInterfaceElement::TPZInterfaceElement material not found\n";
  fLeftEl = left;
  fRightEl = right;
  if(leftside > -1) NormalToFace(fNormal,leftside);
}

void TPZInterfaceElement::CloneInterface(TPZCompMesh *aggmesh,int left,int right) {

  TPZCompElDisc *leftcel = dynamic_cast<TPZCompElDisc *>(aggmesh->ElementVec()[left]);
  TPZCompElDisc *rightcel = dynamic_cast<TPZCompElDisc *>(aggmesh->ElementVec()[right]);
  //copiando ou transferindo o material
  int matid = Material()->Id();
  TPZMaterial *mater = aggmesh->FindMaterial(matid);
  if(!mater){
    //cria copia de material da malha fina    
    mater = Mesh()->FindMaterial(matid);
    mater->Clone(aggmesh->MaterialVec());
//    if( !strcmp(mater->Name(),"TPZEulerConsLaw") ){
//      TPZEulerConsLaw *euler = dynamic_cast<TPZEulerConsLaw *>(mater);
//      mater = new TPZEulerConsLaw(*euler);
//    }
//    aggmesh->InsertMaterialObject(mater);
  }
  int leftside = -1,index;
  TPZGeoEl *gel = Reference();
  //criando a interface na nova malha aglomerada
  TPZInterfaceElement *interf = new TPZInterfaceElement(*aggmesh,gel,index,leftcel,rightcel,leftside);
  //copiando a normal
  TPZVec<REAL> normal(3);
  Normal(normal);
  interf->SetNormal(normal);
}

void TPZInterfaceElement::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){


  TPZConservationLaw *mat = dynamic_cast<TPZConservationLaw *>(fMaterial);
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

  // clean ek and ef
  if(!ek.fMat) ek.fMat = new TPZFMatrix();
  if(!ef.fMat) ef.fMat = new TPZFMatrix();
  if(!ek.fBlock) ek.fBlock = new TPZBlock(ek.fMat);
  if(!ef.fBlock) ef.fBlock = new TPZBlock(ef.fMat);

  int neq = neql + neqr;
  ek.fMat->Redim(neq,neq);
  ef.fMat->Redim(neq,1);
  if(ncon){//no máximo ncon = 1
    ek.fBlock->SetNBlocks(ncon);
    ef.fBlock->SetNBlocks(ncon);
    ek.fBlock->Set(0,left->NShapeF()*nstatel,right->NShapeF()*nstater);
    ef.fBlock->Set(0,left->NShapeF()*nstatel);
  }
  if( !ek.fMat || !ef.fMat || !ek.fBlock || !ef.fBlock){
    cout << "TPZInterfaceElement::CalcStiff : not enough storage for local stifness"
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
  TPZBndCond *bcleft = 0,*bcright=0;

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
    } else {
      bcleft = dynamic_cast<TPZBndCond *> (left->Material());
      if(!bcleft) PZError << "TPZInterfaceElement::CalcStiff material does not exists\n";
      mat->ComputeSolLeft(solr,soll,fNormal,bcleft);
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
    } else {
      bcright = dynamic_cast<TPZBndCond *> (right->Material());
      if(!bcright) PZError << "TPZInterfaceElement::CalcStiff material does not exists\n";
      mat->ComputeSolRight(solr,soll,fNormal,bcright);
    }
    mat->ContributeInterface(x,soll,solr,dsoll,dsolr,weight,fNormal,phixl,phixr,dphixl,dphixr,*ek.fMat,*ef.fMat);
  }
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

int TPZInterfaceElement::NConnects(){ 

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
  out << "\tId of the geometric reference : " << fReference->Id() << endl;
  if(fLeftEl){
    if(fLeftEl->Type() == 17) out << "\tElement index of the left element : "<< fLeftEl->Index() << endl;
    else {
      out << "\tGeometric reference of the left element of id : ";
      out <<  fLeftEl->Reference()->Id() << endl;
    }
  } else {
    out << "Null" << endl;
    cout << "TPZInterfaceElement::Print null left element\n\n";
  }
  if(fRightEl){
    if(fRightEl->Type() == 17) out << "\tElement index of the right element : "<< fRightEl->Index() << endl;
    else {
      out << "\tGeometric reference of the right element of id : ";
      out << fRightEl->Reference()->Id() << endl;
    }
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

  // gInterfaceDimension é a dimensão do elemento de interface
  // verifica-se para cada lado de dimensão gInterfaceDimension do
  // elemento que existe um elemento interface e que este é único

  int iel,iside,nel = cmesh.NElements();

  for(iel=0;iel<nel;iel++){
    TPZCompEl *cel = cmesh.ElementVec()[iel];
    if(!cel) continue;
    TPZGeoEl *geo = cel->Reference();
    if(!geo){
      PZError << "TPZInterfaceElement::main computational element with null reference\n";
      exit(-1);
    }
    int nsides = geo->NSides();;
    for(iside=0;iside<nsides;iside++){
      if(geo->SideDimension(iside) != TPZCompElDisc::gInterfaceDimension) continue;
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
  if(comp.Element()->Type() == 16) exists++;//o próprio é interface
  
  while(neigh.Element() && neigh.Element() != geo.Element()){
    TPZCompElSide comp = neigh.Reference();
    neigh = neigh.Neighbour();
    if(!comp.Element()) continue;
    if(comp.Element()->Type() == 16) exists++;
  }
  if(exists != 1) return 0;
  return 1;//existe uma única interface
}

int TPZInterfaceElement::FreeInterface(TPZCompMesh &cmesh){

  int iel,nel = cmesh.NElements();
  for(iel=0;iel<nel;iel++){
    TPZCompEl *cel = cmesh.ElementVec()[iel];
    if(!cel) continue;
    if(cel->Type() != 16) continue;//interessa só interfaces
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
      if(comp.Element()->Type() != 16) exists++;
    }
    //só pode haver 1 ou 2 elementos de volume associados a um el. interface
    if(exists < 1 || exists > 2) return 0;
  }
  return 1;
}

void VetorialProd(TPZVec<REAL> &ivet,TPZVec<REAL> &jvet,TPZVec<REAL> &kvet);

void TPZInterfaceElement::NormalToFace(TPZVec<REAL> &normal,int leftside){

  //  int dim = fReference->Dimension();
  int face = fReference->NSides()-1;
  //face: lado do elemento bidimensional ou aresta
  //do unidimensional ou canto do ponto
  normal.Resize(3,0.);
  TPZCompElSide neigh(fLeftEl,leftside);
  TPZGeoElSide neighside = neigh.Reference();
  TPZGeoEl *geoneigh = neighside.Element();

  TPZVec<REAL> param(3),cent(3),point(3,0.),result(3,0.),xint(3),xvol(3),vec(3),rib(3);
  TPZFMatrix jacobian(3,3),jacinv(3,3),axes(3,3);
  REAL detjac,normalize;
  int i;

  switch(TPZCompElDisc::gInterfaceDimension){
  case 0:
    normal[0] = 1.0;// a normal sempre apontará na dire¢ão positiva do eixo
    break;
  case 1:
    fReference->CenterPoint(face,param);//ponto interior da interface
    fReference->X(param,xint);
    fReference->Jacobian(param,jacobian,axes,detjac,jacinv);
    face = geoneigh->NSides()-1;//lado interior
    geoneigh->CenterPoint(face,cent);//ponto centro do elemento de volume
    geoneigh->X(cent,xvol);
    for(i=0;i<3;i++) vec[i] = xvol[i]-xint[i];//não deve ser nulo
    for(i=0;i<3;i++) rib[i] = axes(0,i);//dire¢ão da aresta
    VetorialProd(rib,vec,result);
    VetorialProd(result,rib,normal);
    //normalizando a normal
    for(i=0;i<3;i++) vec[i] = 0.;
    normalize = fReference->Distance(normal,vec);
    if(!normalize)
      PZError << "TPZInterfaceElement::NormalToFace null normal vetor\n";
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
  int leftdim = fLeftEl->Dimension();
  int lefts = fLeftEl->Reference()->NSides();
  int rightdim = fRightEl->Dimension();
  int rights = fRightEl->Reference()->NSides();
  TPZVec<REAL> leftp(leftdim), rightp(rightdim), leftx(3),rightx(3);
  fLeftEl->Reference()->CenterPoint(lefts-1,leftp);
  fLeftEl->Reference()->X(leftp,leftx);
  fRightEl->Reference()->CenterPoint(rights-1,rightp);
  fRightEl->Reference()->X(rightp,rightx);
  REAL prod = normal[0]*(rightx[0]-leftx[0])+normal[1]*(rightx[1]-leftx[1])+normal[2]*(rightx[2]-leftx[2]);
  if(prod < 0.) {
  	normal[0] *= -1.;
  	normal[1] *= -1.;
  	normal[2] *= -1.;
  }

//  cout << "Normal : " << normal[0] << " , " << normal[1] << " , " << normal[2] << endl;
}

void VetorialProd(TPZVec<REAL> &ivet,TPZVec<REAL> &jvet,TPZVec<REAL> &kvet){

  kvet.Resize(3);
  kvet[0] =  ivet[1]*jvet[2] - ivet[2]*jvet[1];
  kvet[1] = -ivet[0]*jvet[2] + ivet[2]*jvet[0];
  kvet[2] =  ivet[0]*jvet[1] - ivet[1]*jvet[0];
}

void TPZInterfaceElement::Normal(TPZVec<REAL> &normal){

  for(int i=0;i<3;i++) normal[i] = fNormal[i];
}

void TPZInterfaceElement::SetNormal(TPZVec<REAL> &normal){

  for(int i=0;i<3;i++) fNormal[i] = normal[i];
}

/**
 *  ESTA ALTERA¢ÃO É SÓ PARA ATUALIZAR O PZREPOSITORY (CVS)
 */
  //achando o lado do vizinho esquerdo
//   int side = fReference->NSides()-1;
//   TPZGeoElSide refside(fReference,side);
//   TPZGeoElSide neigh = refside.Neighbour();
//   TPZGeoEl *gel = neigh.Element();
//   while(gel && neigh != refside){
//     TPZCompEl *cel = gel->Reference();
//     if(cel && cel == LeftElement()) break;
//     neigh = neigh.Neighbour();
//     gel = neigh.Element();
//   }
//   if(!gel || neigh == refside){
//     PZError << "TPZInterfaceElement::CloneInterface left element not found\n";
//     exit(-1);
//   }
//  int leftside = neigh.Side(),index;
