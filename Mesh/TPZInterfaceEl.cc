
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzquad.h"
#include "pzmaterial.h"
#include "TPZConservationLaw.h"

TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompEl &thirdel) 
  : TPZCompEl(mesh,index){

  fReference = geo;
  geo->SetReference(this);
  fLeftEl = NULL;
  fRightEl = NULL;
  VolumeEls(thirdel);//identifica elementos esquerdo e direito conectados
}

TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index) 
  : TPZCompEl(mesh,index){

  fReference = geo;
  geo->SetReference(this);
  fLeftEl = NULL;
  fRightEl = NULL;
  //VolumeEls();//identifica elementos esquerdo e direito conectados
}

void TPZInterfaceElement::VolumeEls(TPZCompEl &thirdel){

  /**A ORDEM DE DEFINICÃO DOS VÉRTICES DO ELEMENTO INTERFACE 
     DETERMINA QUEM SÃO OS ELEMENTOS ESQUERDO E DIREITO ASSOCIADOS*/
  REAL detjac;
  TPZVec<REAL> param(3),normal(3);
  TPZFMatrix jacobian(3,3),jacinv(3,3),axes(3,3);
  int face = fReference->NSides()-1;//face: lado do elemento bidimensional ou aresta do unidimensional
  //cada face tem no máximo 3 elementos ligados, um de interface + esquerdo e direito
  //se a face é de fronteira um dos elementos é de contorno
  //pelos menos 1 é de volume
  fReference->CenterPoint(face,param);//ponto da face  
  fReference->Jacobian(param,jacobian,axes,detjac,jacinv);//normal: 3a linha de axes
  TPZStack<TPZCompElSide> list;
  list.Resize(0);
  //a próxima linha não retorna um el. interface
  TPZCompElSide(this,face).EqualLevelElementList(list,0,1);
  TPZGeoElSide gs = list[0].Reference();
  if(!gs.Exists())
    PZError << "TPZInterfaceElement::VolumeEls neighbor does not exist, inconsistency of data\n";
  TPZGeoEl *neigh = gs.Element(),*neigh2;
  int dim = fReference->Dimension();//dimensão do atual que é 2D
  neigh2 = thirdel.Reference();
  if(neigh == fReference || neigh2 == fReference){
    cout << "TPZInterfaceElement::VolumeEls error data (nao acha comp. de volume no ciclo: impossivel)\n";
    exit(-1);
  }
  TPZVec<REAL> x0(3);
  fReference->X(param,x0);//ponto da interface
  param.Resize(2);
  //elemento de volume associado
  TPZVec<REAL> x1;
  TPZCompElDisc *neighdisc = (TPZCompElDisc *) neigh->Reference();
  //aqui neighdisc for¢osamente existe
  if (neighdisc){
    neighdisc->InternalPoint(x1);//ponto interior ao volume
    TPZVec<REAL> vec(3);
    int i;
    for(i=0;i<3;i++) vec[i] = x1[i]-x0[i];//não deve ser nulo
    REAL prod = vec[0]*axes(2,0)+vec[1]*axes(2,1)+vec[2]*axes(2,2);//se pord = 0 os vetores são paralelos (superpostos)
    if(prod < 0)//ângulo maior que 90
      fLeftEl = neighdisc;
    else
      fRightEl = neighdisc;//a normal a interface aponta para o elemento de volume
    //segundo elemento de volume associado
    neighdisc = (TPZCompElDisc *) neigh2->Reference();
    if(!fLeftEl) fLeftEl = neighdisc;
    else fRightEl = neighdisc;
  }
  if(!fLeftEl || !fRightEl)
    PZError << "TPZInterfaceElement::VolumeEls not identified left or right element\n";
}

void TPZInterfaceElement::CalcStiff(TPZFMatrix &ek,TPZFMatrix &ef){

  TPZCompElDisc *left = LeftElement();
  TPZCompElDisc *right = RightElement();
  if(!left->Material() || right->Material()){
    PZError << "TPZInterfaceElement::CalcStiffInterf null material\n";
    return;
  }
  int nshapel = left->NShapeF();
  int nshaper = right->NShapeF();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix &MeshSol = Mesh()->Solution();
  int nstate = left->Material()->NStateVariables();
  int neql = nshapel * nstate;
  int neqr = nshaper * nstate;
  int neqsum = neql+neqr;
  int dim = Dimension();
  ek.Redim(neql,neqr);
  ef.Redim(neql,1);
  TPZFMatrix phixl(neql,1,0.),dphixl(dim,neql);
  TPZFMatrix phixr(neqr,1,0.),dphixr(dim,neqr);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL detjac,weight;
  TPZVec<REAL> soll(nstate,0.),solr(nstate,0.);
  TPZFMatrix dsoll(dim,nstate,0.),dsolr(dim,nstate,0.);
  int pl = left->Degree();
  int pr = right->Degree();
  int p = (pl > pr) ? pl : pr;
  int face = fReference->NSides()-1;
  TPZIntPoints *intrule = fReference->CreateSideIntegrationRule(face,p);
  int npoints = intrule->NPoints();
  //calcular transformaões entre face e volume
  TPZGeoEl *refl = left->Reference();
  TPZGeoEl *refr = right->Reference();
  TPZGeoElSide leftside(refl,refl->NSides()-1);
  TPZGeoElSide rightside(refr,refr->NSides()-1);
  TPZTransform tl(leftside.Dimension()),tr(rightside.Dimension());
  TPZGeoElSide thisgeoside(fReference,face);
  thisgeoside.SideTransform3(leftside,tl);
  thisgeoside.SideTransform3(rightside,tr);
  int ip,k,j;
  TPZVec<REAL> point(3,0.);
  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    weight *= fabs(detjac);
    tl.Apply(intpoint,point);
    fReference->X(point, x);
    left->Shape(x,phixl,dphixl);
    tr.Apply(intpoint,point);
    fReference->X(point, x);
    right->Shape(x,phixr,dphixr);
    //solu¢ão da itera¢ão anterior
    soll.Fill(0.);
    dsoll.Zero();
    TPZConnect *df = &left->Connect(0);
    int dfseq = df->SequenceNumber();
    int dfvar = block.Size(dfseq);
    int pos = block.Position(dfseq);
    int iv = 0,d;
    for(int jn=0; jn<dfvar; jn++) {
      soll[iv%nstate] += phixl(iv/nstate,0)*MeshSol(pos+jn,0);
      for(d=0; d<dim; d++)
	dsoll(d,iv%nstate) += dphixl(d,iv/nstate)*MeshSol(pos+jn,0);
      iv++;
    }
    //solu¢ão da itera¢ão anterior
    solr.Fill(0.);
    dsolr.Zero();
    df = &right->Connect(0);
    dfseq = df->SequenceNumber();
    dfvar = block.Size(dfseq);
    pos = block.Position(dfseq);
    iv = 0;
    for(int jn=0; jn<dfvar; jn++) {
      solr[iv%nstate] += phixr(iv/nstate,0)*MeshSol(pos+jn,0);
      for(d=0; d<dim; d++)
	dsolr(d,iv%nstate) += dphixr(d,iv/nstate)*MeshSol(pos+jn,0);
      iv++;
    }
    TPZConservationLaw *mat = (TPZConservationLaw *) fMaterial;
    if(!mat){
      PZError << "TPZInterfaceElement::CalcStiff interface material null, do nothing\n";
      return;
    }
    mat->ContributeInterface(x,soll,solr,dsoll,dsolr,weight,axes,phixl,phixr,dphixl,dphixr,ek,ef);
  }
}

/*
void TPZCompElDisc::CalcStiffDisc(TPZFMatrix &ek, TPZFMatrix &ef){

  if(fMaterial == NULL){
    cout << "TPZCompElDisc::CalcStiff : no material for this element\n";
    return;
  }

  int dim = Dimension();
  int nstate = fMaterial->NStateVariables();
  int nshape = NShapeF();
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
    TPZConnect *df = &Connect(in);
    int dfseq = df->SequenceNumber();
    int dfvar = block.Size(dfseq);
    int pos = block.Position(dfseq);
    for(int jn=0; jn<dfvar; jn++) {
      sol[iv%nstate] += phi(iv/nstate,0)*MeshSol(pos+jn,0);
      for(d=0; d<dim; d++)
	dsol(d,iv%nstate) += dphix(d,iv/nstate)*MeshSol(pos+jn,0);
      iv++;
    }
    fMaterial->Contribute(x,jacinv,sol,dsol,weight,axes,phix,dphix,ek,ef);
  }
}
*/

int TPZInterfaceElement::NConnects(){ 

  if(!fLeftEl && !fRightEl) return 0;
  if(!fLeftEl || !fRightEl) return 1;
  return 2;//that's the right and left element
}

int TPZInterfaceElement::ConnectIndex(int i) {


  if(i<0 || i>1)
    PZError << "TPZInterfaceElement::ConnectIndex argument i error, i= " << i << endl;

  if(i == 0 && fLeftEl){
    return fLeftEl->ConnectIndex(0);
  }
  if(i == 1 && fRightEl)
   return fRightEl->ConnectIndex(0);

  return -1;//ambos nulos, esquerdo e direito
}

void TPZInterfaceElement::Print(ostream &out){

  out << "\nInterface element : \n";
  out << "\tGeometric reference id : " << fReference->Id() << endl;
  out << "\tGemotric reference id of the left  element : ";
  if(fLeftEl){
    out <<  fLeftEl->Reference()->Id() << endl;
  } else {
    out << "Null" << endl;
    cout << "TPZInterfaceElement::Print null left element\n\n";
  }
  out << "\tGemotric reference id of the right element : ";
  if(fRightEl){
    out << fRightEl->Reference()->Id() << endl;
  } else {
    out << "Null" << endl;
    cout << "TPZInterfaceElement::Print null right element\n\n";
  }
  out << "\tMaterial id : " << fReference->MaterialId() << endl;

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
	PZError << "number of interfaces existing: " << nint << endl;
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
  if(comp.Element()->Type() == 17) exists++;//o próprio é interface
  
  while(neigh.Element() && neigh.Element() != geo.Element()){
    TPZCompElSide comp = neigh.Reference();
    neigh = neigh.Neighbour();
    if(!comp.Element()) continue;
    if(comp.Element()->Type() == 17) exists++;
  }
  if(exists != 1) return 0;
  return 1;//existe uma única interface
}

int TPZInterfaceElement::FreeInterface(TPZCompMesh &cmesh){

  int iel,nel = cmesh.NElements();
  for(iel=0;iel<nel;iel++){
    TPZCompEl *cel = cmesh.ElementVec()[iel];
    if(!cel) continue;
    if(cel->Type() != 17) continue;//interessa só interfaces
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
      if(comp.Element()->Type() != 17) exists++;
    }
    //só pode haver 1 ou 2 elementos de volume associados a um el. interface
    if(exists < 1 || exists > 2) return 0;
  }
  return 1;
}

void TPZInterfaceElement::NormalToFace(TPZVec<REAL> &normal){

  REAL detjac;
  TPZVec<REAL> param(3);
  TPZFMatrix jacobian(3,3),jacinv(3,3),axes(3,3);
  int face = fReference->NSides()-1;//face: lado do elemento bidimensional ou aresta do unidimensional
  //cada face tem no máximo 3 elementos ligados, um de interface + esquerdo e direito
  //se a face é de fronteira um dos elementos é de contorno
  //pelos menos 1 é de volume
  fReference->CenterPoint(face,param);//ponto da face  
  fReference->Jacobian(param,jacobian,axes,detjac,jacinv);
  int i;
  for(i=0;i<3;i++) normal[i] = axes(2,i);
  //vetor que sai do elemento esquerdo e entra no 
  //elemento direito (veja TPZInterfaceElement::VolumeEls)

}
