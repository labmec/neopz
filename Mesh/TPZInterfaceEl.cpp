
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzquad.h"


TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index) 
  : TPZCompEl(mesh,index){

  fReference = geo;
  geo->SetReference(this);
  fLeftEl = NULL;
  fRightEl = NULL;
  VolumeEls();//identifica elementos esquerdo e direito conectados
}

/* TPZInterfaceElement::TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoElT2d *geo,int &index)  */
/*   : TPZCompEl(mesh,index){ */

/*   fReference = geo; */
/*   geo->SetReference(this); */
/*   fLeftEl = NULL; */
/*   fRightEl = NULL; */
/*   VolumeEls();//identifica elementos esquerdo e direito conectados */
/* } */

void TPZInterfaceElement::VolumeEls(){

  /**A ORDEM DE DEFINICÃO DOS VÉRTICES DO ELEMENTO INTERFACE 
     DETERMINA QUEM SÃO OS ELEMENTOS ESQUERDO E DIREITO ASSOCIADOS*/
  REAL detjac;
  TPZVec<REAL> param(3),normal(3);
  TPZFMatrix jacobian(3,3),jacinv(3,3),axes(3,3);
  int face = fReference->NSides()-1;//face: lado do elemento bidimensional
  //cada face tem no máximo 3 elementos ligados um de interface + esquerdo e direito
  //se a face é de fronteira um dos elementos é de contorno
  //pelos menos 1 é de volume, o cíclo é:  this->disc1->disc2->this
  fReference->CenterPoint(face,param);//ponto da face  
  fReference->Jacobian(param,jacobian,axes,detjac,jacinv);//normal: 3a linha de axes
  TPZGeoElSide gs = fReference->Neighbour(face);
  TPZGeoEl *neigh = gs.Element(),*neigh2;
  int dim = fReference->Dimension();//dimensão do atual que é 2D
  gs = gs.Neighbour();
  neigh2 = gs.Element();
/*   if(neigh->Dimension() == dim){ */
/*     TPZGeoELSide *neighaux = neigh; */
/*     neigh = neigh2; */
/*     neigh2 = neighaux; */
/*   } */
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
}

void TPZInterfaceElement::CalcStiffInterf(TPZFMatrix &jacob,TPZFMatrix &res){

  TPZCompElDisc *left = LeftElement();
  TPZCompElDisc *right = RightElement();
  if(!left->Material() || right->Material()){
    PZError << "TPZInterfaceElement::CalcStiffInterf null material\n";
    return;
  }
  int nshapel = left->NShapeF();
  int nshaper = right->NShapeF();
  int nstate = left->Material()->NStateVariables();
  int neql = nshapel * nstate;
  int neqr = nshaper * nstate;
  int neqsum = neql+neqr;
  int dim = Dimension();
  jacob.Resize(neql,neqr);
  res.Resize(neql,1);
  TPZFMatrix phix(neql+neqr,1,0.),dphix(dim,neql+neqr);
  TPZFMatrix phixr(neqr,1,0.),dphixr(dim,neqr);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL detjac,weight;
  TPZVec<REAL> sol(nstate,0.);
  TPZFMatrix dsol(dim,nstate,0.);
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
  TPZVec<REAL> volpoint(3,0.);
  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    weight *= fabs(detjac);
    tl.Apply(intpoint,volpoint);
    fReference->X(volpoint, x);
    left->Shape(x,phix,dphix);
    tr.Apply(intpoint,volpoint);
    fReference->X(volpoint, x);
    right->Shape(x,phixr,dphixr);
    for(k=0;k<neqr;k++){
      j = k + neql;
      phix(j,0) = phixr(k,0);
    }
    fMaterial->Contribute(x,jacinv,sol,dsol,weight,axes,phix,dphix,jacob,res);
    //o número de linhas de ek indica o número de shapes do elemento esquerdo
    //o número de coluas de ek indica o número de shapes do elemento direito
  }
}

int TPZInterfaceElement::NConnects(){ 

  if(!fLeftEl && !fRightEl) return 0;
  if(!fLeftEl || !fRightEl) return 1;
  return 2;//that's the right and left element
}

int TPZInterfaceElement::ConnectIndex(int i) {


  if(i<0 || i>1)
    PZError << "TPZInterfaceElement::ConnectIndex argument i error, i= " << i << endl;

  if(i == 0 && fLeftEl){
    int k = fLeftEl->ConnectIndex(0);
    return k;
  }
  if(i == 1 && fRightEl)
   return fRightEl->ConnectIndex(0);

  return -1;//ambos nulos, esquerdo e direito
}

void TPZInterfaceElement::Print(ostream &out){

  out << "\nInterface element : \n";
  out << "\tGeometric reference id : " << fReference->Id() << endl
      << "\tGemotric reference id of the left  element : " <<  fLeftEl->Reference()->Id() << endl
      << "\tGemotric reference id of the right element : " << fRightEl->Reference()->Id() << endl
      << "\tMaterial id : " << fReference->MaterialId() << endl;

}



  //TPZGeoElSide::SideTransform2(TPZGeoElSide neighbour,TPZTransform &t);
  //TPZTransform t(thisgeoside.Dimension());
  //thisgeoside.SideTransform2(largeside,t);

/*   gs = gs.Neighbour();//próximo vizinho */
/*   neigh = gs.Element();//geométrico de neighdisc */
/*   TPZGeoEl *ref = neighdisc->Reference(); */
/*   while(neigh && neigh != ref){ */
/*     if(neigh->Reference() && neigh->Reference() != this) break;//computacional de volume */
/*     gs = gs.Neighbour(); */
/*     neigh = gs.Element(); */
/*   } */
  //si um dos elementos esquerdo ou direito não foi criado computacionalmente ainda 
  //quando seja a vez dele o elemento que falta será atualizado
/*   if(neigh == ref){ */
/*     cout << "TPZInterfaceElement::VolumeEls error data (nao acha comp. de volume no ciclo: impossivel)\n"; */
/*     exit(-1); */
/*   } */

void TPZInterfaceElement::SetConnectIndex(int node, int index) {
  cout << "TPZInterfaceElement::SetConnectIndex should never be called\n";
}
