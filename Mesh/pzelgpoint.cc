// -*- c++ -*-
/**METHODS DEFINITION FOR CLASS ELEM1D*/

#include "pzelgpoint.h"
#include "pzgnode.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzelcpoint.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include "pzrefpoint.h"
#include <math.h>


static TPZCompEl *CreateEl(TPZGeoElPoint *gel, TPZCompMesh &mesh, int &index) {
  return new TPZCompElPoint(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElPoint::fp)(TPZGeoElPoint *,TPZCompMesh &,int &) = CreateEl;

TPZGeoElPoint::TPZGeoElPoint(int id,TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh)
  : TPZGeoEl(id,matind,mesh) {

  int nnod = nodeindices.NElements();
  if(nnod != 1) {
    PZError << "TPZGeoElPoint Constuctor, number of nodes : " << nnod << endl;
    return;
  }
  fNodeIndexes = nodeindices[0];
  fSubEl = 0;

}

TPZGeoElPoint::TPZGeoElPoint( TPZVec<int> &nodeindices, int matind, TPZGeoMesh &mesh)
  : TPZGeoEl(matind,mesh) {

  int nnod = nodeindices.NElements();
  if(nnod != 1) {
    PZError << "TPZGeoElPoint Constuctor, number of nodes : " << nnod << endl;
    return;
  }
  fNodeIndexes = nodeindices[0];
  fSubEl = 0;
}

TPZGeoElPoint *TPZGeoElPoint::CreateGeoEl(TPZVec<int> &np, int matind, TPZGeoMesh &mesh) {
  return new TPZGeoElPoint(np,matind,mesh);
}

TPZGeoElPoint::~TPZGeoElPoint() {
}

int TPZGeoElPoint::NSubElements() {
	return 0;
}

void TPZGeoElPoint::Jacobian(TPZVec<REAL> &/*fl*/,TPZFMatrix &result,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){

   result.Redim(1,1);
   jacinv.Redim(1,1);
   result(0,0) = 1.;
   detjac = 1.;
   jacinv(0,0) = 1.;
   axes.Zero();
   axes.Redim(3,3);
   axes.Zero();
   axes(0,0) = 1.;
   axes(1,1) = 1.;
   axes(2,2) = 1.;
}


void TPZGeoElPoint::X(TPZVec<REAL> & /*par*/, TPZVec<REAL> &result){

  int j;
  for(j=0; j<3; j++)  result[j] = NodePtr(0)->Coord(j);
}

int TPZGeoElPoint::NSideNodes(int side) {
  if(side == 0) return 1;
  return 0;
}

int TPZGeoElPoint::SideNodeIndex(int side, int /*node*/){
  switch(side) {
     case 0:
       return fNodeIndexes;
     default:
       PZError << "TPZGeoElPoint::SideNodeIndex. Bad parameter side.\n";
     }
  return -1;
}

void TPZGeoElPoint::MidSideNodeIndex(int side,int &index) {
  switch(side) {
  case 0:
    index = fNodeIndexes;
    return;
  default:
    PZError << "TPZGeoElPoint::MidSideNode called for side " << side << endl;
    PZError.flush();
    index = -1;
  }
}

int TPZGeoElPoint::SideDimension(int side) {
  switch(side) {
  case 0:
    return 0;
  default:
    PZError << "TPZGeoElPoint::SideDimension. Bad parameter side.\n";
    return -1;
  }
}

int TPZGeoElPoint::NodeIndex(int node) {
  if(node == 0) return fNodeIndexes;
  return -1;
}

void TPZGeoElPoint::Divide(TPZVec<TPZGeoEl *> &pv) {
  fSubEl = this;
  pv.Resize(1);
}

void TPZGeoElPoint::CreateNewNodes(int *gnodindex) {

  //int numnod = 1;
  int midsidenodeindex;
  MidSideNodeIndex(2,midsidenodeindex);
  TPZGeoNode *midsidenode;
  if(midsidenodeindex == -1) {
    //first look over all neighbours the middle node        // falta pesquisar nos vizinhos (1-Jorge)
    TPZGeoElSide neigh = Neighbour(2);
    while(neigh.Exists() && this!=neigh.Element()) {
      neigh.Element()->MidSideNodeIndex(neigh.Side(),midsidenodeindex);
      if(midsidenodeindex!=-1) break;
      neigh = neigh.Neighbour();
    }
    if(midsidenodeindex==-1) {
      TPZVec<REAL> gco(3);
      TPZVec<REAL> par(1);
      par[0] = 0.;
      X(par,gco);
      midsidenodeindex = Mesh()->NodeVec().AllocateNewElement();
      midsidenode = &Mesh()->NodeVec()[midsidenodeindex];
      midsidenode->Initialize(gco,*Mesh());
    }
  }
  //gnodindex[0] = NodeIndex(0);
  //gnodindex[2] = midsidenodeindex;
  gnodindex[0] = midsidenodeindex;
  //gnodindex[4] = NodeIndex(numnod-1);

}

void TPZGeoElPoint::NormalVector(int /*side*/, TPZVec<REAL> &/*loc*/, TPZVec<REAL> &normal,
			      TPZFMatrix &axes, TPZFMatrix &/*jac*/){
    axes(0,0) = 1.;
    axes(1,1) = 1.;
    axes(2,2) = 1.;
    normal[0] = 0.;
    normal[1] = 0.;
    normal[2] = 0.;
}

/**Accumulates the transformation of the jacobian which maps the current
master element space into the space of the master element of the father*/
//void TPZGeoElPoint::BuildTransform(int /*side*/, TPZGeoEl * /*father*/, TPZTransform &/*t*/) {
//	cout <<  "\nTPZGeoElPoint::BuildTransform called\n";
//}

TPZTransform TPZGeoElPoint::SideToSideTransform(int sidefrom, int sideto) {
  if(sideto != 0 && sidefrom !=0) {
    PZError << "TPZGeoElPoint:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto << endl;
  }
  TPZTransform t(1,0);                                                          //   !!!   ???
  t.Sum()(0,0) = .0;
  return t;
}

int TPZGeoElPoint::NSideSubElements(int side) {
  if(side==0) return 0;
  return -1;
}



TPZGeoElSide TPZGeoElPoint::SideSubElement(int side,int /*subel*/) {
  if(side != 0) {
    PZError << "TGeoElQ2d::SideSubElement called for side " << side << endl;
    return TPZGeoElSide(0,0);
  }
  if(!fSubEl) return TPZGeoElSide(0,0);
  return TPZGeoElSide(fSubEl,side);
}

/**Neigh é o lado de um elemento vizinho ao elemento atual El que esta sendo dividido rfnd guarda
   os dois nós globais deste lado no sentido antihorario do elemento El, posições 0 e 1 de rfnd*/
void TPZGeoElPoint::GetSubElement(int side,TPZVec<int> &/*rfndindex*/,TPZVec<TPZGeoElSide> &sub) {
  if(!HasSubElement()) {
    sub.Resize(0);
    return;
  }
  sub[0] = TPZGeoElSide();
}



TPZGeoElSide TPZGeoElPoint::Father(int /*side*/) {
  if(!fFather) return TPZGeoElSide();
  return TPZGeoElSide();
}

void TPZGeoElPoint::LowerDimensionSides(int side,TPZStack<int> &smallsides) {

}

TPZGeoElSide TPZGeoElPoint::HigherDimensionSides(int /*side*/,int /*targetdimension*/) {
   return TPZGeoElSide();
}
void TPZGeoElPoint::AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides){
	return;
}

/**Inicializa os coeficientes do par de nós do lado I do elemento de referencia*/
//void TPZGeoElPoint::SideMasterCo(int /*side*/,TPZVec<REAL> &IVec,TPZVec<REAL> &JVec) {

//  IVec.Fill(0.,0);
//  JVec.Fill(0.,0);
//}

//void TPZGeoElPoint::SideMasterCo(int side,TPZFMatrix &coord) {
//  if(side == 0) coord(0,0) = 1.;

//}

TPZGeoEl *TPZGeoElPoint::CreateBCGeoEl(int side, int bc) {

  if(side==0) {
    TPZManVector<int> nodes(1);
    nodes[0] = fNodeIndexes;
    TPZGeoElPoint *gel = CreateGeoEl(nodes,bc,*Mesh());
    TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,0));
    return gel;
  }
  else PZError << "TPZGeoElPoint::CreateBCCompEl. Side = " << side << endl;
  return 0;
}

TPZGeoElSide TPZGeoElPoint::Father2(int /*side*/){//Augusto:09/01/01
  //  cout << "TPZGeoElPoint::Father2 to be implemented\n";
  return TPZGeoElSide();
}

int TPZGeoElPoint::FatherSide(int side, int son){
	return TPZRefPoint::FatherSide(side,son);
}
 
TPZTransform TPZGeoElPoint::BuildTransform2(int /*side*/, TPZGeoEl * /*father*/, TPZTransform & /*t*/){//Augusto:09/01/01
  cout << "TPZGeoElPoint::BuildTransform2 makes no sense\n";
  return TPZTransform(0,0);
}


void TPZGeoElPoint::GetSubElements2(int /*side*/, TPZStack<TPZGeoElSide> &/*subel*/){//Augusto:09/01/01
  cout << "TPZGeoElPoint::GetSubElements2 to be implemented\n";
}

int TPZGeoElPoint::NSideSubElements2(int side) {
    cout << "TPZGeoElPoint::NSideSubElements2 to be given a meaning\n";
	if(side == 0) return 1;
	else return -1;
}

void TPZGeoElPoint::SetSubElement(int id, TPZGeoEl *el){
  if (id != 0){
    PZError << "TPZGeoElPoint::Trying do define subelement :" << id << endl;
    return;
  }
  fSubEl=el;
  return;
}
  

TPZIntPoints * TPZGeoElPoint::CreateSideIntegrationRule(int /*side*/, int /*order*/)
{
   return new TPZInt1Point();

}


TPZTransform TPZGeoElPoint::GetTransform(int side,int son) {
	return TPZRefPoint::GetTransform(side,son);
}

void TPZGeoElPoint::CenterPoint(int side, TPZVec<REAL> &masscent){

  if(side < 0 || side > NSides()-1){
    PZError << "TPZGeoElPoint::CenterPoint error side = " << side << endl;
    return;
  }

  masscent[0] = 0;//NodePtr(0)->Coord(0);
  masscent[1] = 0;//NodePtr(0)->Coord(1);
  masscent[2] = 0;//NodePtr(0)->Coord(2);
}
