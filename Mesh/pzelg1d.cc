//METHODS DEFINITION FOR CLASS ELEM1D

// -*- c++ -*-
#include "pzelg1d.h"
#include "pzelgpoint.h"
#include "pzgnode.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzelc1d.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include <math.h>
#include "TPZRefLinear.h"
#include "TPZGeoLinear.h"


static TPZCompEl *CreateEl(TPZGeoEl1d *gel, TPZCompMesh &mesh, int &index) {
  return new TPZCompEl1d(mesh,gel,index);
}

TPZCompEl *(*TPZGeoEl1d::fp)(TPZGeoEl1d *,TPZCompMesh &,int &) = CreateEl;

TPZGeoEl1d::TPZGeoEl1d(int id,TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,int refind)
  : TPZGeoEl(id,matind,mesh) {

  fYAxisIndex = refind;
  int nnod = nodeindices.NElements();
  if(nnod <2 || nnod > 2) {
    PZError << "TPZGeoEl1d Constuctor, number of nodes : " << nnod << endl;
    return;
  }
  fNodeIndexes[0] = nodeindices[0];
  fNodeIndexes[1] = nodeindices[1];
  //if(nnod == 3) fNodeIndexes[2] = nodeindices[2];
  fNodeIndexes[2] = -1;
  fSubEl[0] = 0;
  fSubEl[1] = 0;
}

TPZGeoEl1d::TPZGeoEl1d( TPZVec<int> &nodeindices, int matind, TPZGeoMesh &mesh, int refind)
  : TPZGeoEl(matind,mesh) {

  fYAxisIndex = refind;
  int nnod = nodeindices.NElements();
  if(nnod <2 || nnod > 2) {
    PZError << "TPZGeoEl1d Constuctor, number of nodes : " << nnod << endl;
    return;
  }
  fNodeIndexes[0] = nodeindices[0];
  fNodeIndexes[1] = nodeindices[1];
  //if(nnod == 3) fNodeIndexes[2] = nodeindices[2];
  fNodeIndexes[2] = -1;
  fSubEl[0] = 0;
  fSubEl[1] = 0;
}

TPZGeoEl1d *TPZGeoEl1d::CreateGeoEl(TPZVec<int> &np, int matind, TPZGeoMesh &mesh, int refind) {
  return new TPZGeoEl1d(np,matind,mesh,refind);
}

TPZGeoEl1d::~TPZGeoEl1d() {
}

int TPZGeoEl1d::NSubElements() {
	return TPZRefLinear::NSubEl;
}

void TPZGeoEl1d::Jacobian(TPZVec<REAL> &fl,TPZFMatrix &result,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){

	TPZFMatrix xco(3,TPZShapeLinear::NNodes);
	int i,j;
	for(i=0; i<TPZShapeLinear::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoLinear::Jacobian(xco,fl,result,axes,detjac,jacinv);
}


void TPZGeoEl1d::X(TPZVec<REAL> & par, TPZVec<REAL> &result){

	TPZFMatrix xco(3,TPZShapeLinear::NNodes);
	int i,j;
	for(i=0; i<TPZShapeLinear::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoLinear::X(xco,par,result);
}

int TPZGeoEl1d::NSideNodes(int side) {
	return TPZShapeLinear::NSideNodes(side);
}

int TPZGeoEl1d::SideNodeIndex(int side, int node){
	int locind = TPZShapeLinear::SideNodeLocId(side,node);
	return NodeIndex(locind);
}

void TPZGeoEl1d::MidSideNodeIndex(int side,int &index) {
	TPZRefLinear::MidSideNodeIndex(this,side,index);
}

int TPZGeoEl1d::SideDimension(int side) {
	return TPZShapeLinear::SideDimension(side);
}

int TPZGeoEl1d::NodeIndex(int node) {
  if(0 <= node && node< 2) return fNodeIndexes[node];
  return -1;
}

void TPZGeoEl1d::Divide(TPZVec<TPZGeoEl *> &pv) {
	TPZRefLinear::Divide(this,pv);
}

/*
void TPZGeoEl1d::CreateNewNodes(int *gnodindex) {
  
  int numnod = NNodes();
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
  gnodindex[0] = NodeIndex(0);
  gnodindex[2] = midsidenodeindex;
  gnodindex[4] = NodeIndex(numnod-1);
  if(numnod == 3) {
    // look for the midsidenodes of the sons of the neighbouring elements
    TPZVec<int> nodid(2),subelside(2);
    int side = 2;
    TPZVec<TPZGeoElSide> subelements(2);
    nodid[0] = NodeIndex(0);
    nodid[1] = NodeIndex(1);
    TPZGeoElSide neighbour;
    neighbour = Neighbour(side);//this e o elemento atual que esta sendo dividido
    if(neighbour.Exists()) {
      // try to find the midsidenodes of the neighbouring elements
      while(this != neighbour.Element()) {
	if(neighbour.HasSubElement()) {
	  neighbour.GetSubElement(nodid, subelements);
	  TPZGeoEl *gel = subelements[0].Element();
	  if(gnodindex[1]==-1 && gel) gel->MidSideNodeIndex(subelements[0].Side(),gnodindex[1]);
	  gel = subelements[1].Element();
	  if(gnodindex[3]==-1 && gel) gel->MidSideNodeIndex(subelements[1].Side(),gnodindex[3]);
	}
	neighbour = neighbour.Neighbour();
      }
    }
    if(gnodindex[1] != -1 && gnodindex[3] != -1) return;
    TPZVec<REAL> gco(3);
    TPZVec<REAL> par(1);
    if(gnodindex[1] == -1) {
      par[0] = -0.5;
      X(par,gco);
      gnodindex[1] = Mesh()->NodeVec().AllocateNewElement();
      TPZGeoNode *gnodptr = &Mesh()->NodeVec()[gnodindex[1]];
      gnodptr->Initialize(gco,*Mesh());
    }
    if(gnodindex[3] == -1) {
      par[0] = 0.5;
      X(par,gco);
      gnodindex[3] = Mesh()->NodeVec().AllocateNewElement();
      TPZGeoNode *gnodptr = &Mesh()->NodeVec()[gnodindex[3]];
      gnodptr->Initialize(gco,*Mesh());
    }
  }
}
*/
/*
void TPZGeoEl1d::NormalVector(int side, TPZVec<REAL> &loc, TPZVec<REAL> &normal,
			      TPZFMatrix &axes, TPZFMatrix &jac){
  
  REAL spacephi[9], spacedphi[18];
  int numnodes = NNodes();
  TPZFMatrix phi(numnodes,1,spacephi,9),dphi(1,numnodes,spacedphi,18);
  Shape1d(loc[0],numnodes,phi,dphi);                                     //Cuidado com Shape1d (orientacao) (3-Jorge)
  TPZGeoNode *np;
  int ic;
  TPZVec<REAL> yaxis(3,0.),v1(3,0.), vi(3,0.), v2(3,0.), v3(3,0.);
  REAL mod1 = 0;
  REAL modi = 0;
  REAL mod2 = 0;
  if(fYAxisIndex != -1) {
    TPZGeoNode &node = Mesh()->NodeVec()[fYAxisIndex];
    yaxis[0] = node.Coord(0);
    yaxis[1] = node.Coord(0);
    yaxis[2] = node.Coord(0);
  }
  
  for(int i=0; i < numnodes; i++) {
    np = NodePtr(i);
    for(ic = 0; ic < 3; ic++) {
      v1[ic] += np->Coord(ic)*dphi(0,i);
      vi[ic] = yaxis[i];
    }
  }

  for(ic=0; ic<3; ic++) {
    mod1 += v1[ic]*v1[ic];
    modi += vi[ic]*vi[ic];
  }
  
  mod1 = sqrt(mod1);
  modi = sqrt(modi);
  
  for(ic=0; ic < 3; ic++) {
    v1[ic] = v1[ic]/mod1;
    if(modi) vi[ic] = vi[ic]/modi;
  }

  REAL kk = 0;
  for(ic=0; ic<3; ic++) kk += vi[ic]*v1[ic];

  for(ic=0; ic<3; ic++) {
    v2[ic] = vi[ic]-kk*v1[ic];
  }
  
  for(ic=0; ic<3; ic++) {
    mod2 += v2[ic]*v2[ic];
  }
  
  mod2 = sqrt(mod2);

  if(mod2) {
    for(ic=0; ic<3; ic++) {
      v2[ic] = v2[ic]/mod2;
    }
  }
  
  v3[0] = v1[1]*v2[2]-v2[1]*v1[2];
  v3[1] = v1[2]*v2[0]-v1[0]*v2[2];
  v3[2] = v1[0]*v2[1]-v1[1]*v2[0];
  
  for(ic=0; ic<3; ic++) {
    axes(0, ic) = v1[ic];
    axes(1, ic) = v2[ic];
    axes(2, ic) = v3[ic];
    normal[ic] = v1[ic];
  }
  switch(side) {
  case 1:
    normal[0] *= -1.;
    normal[1] *= -1.;
    normal[2] *= -1.;
    break;
  case 0:
  case 2:
  default:
    break;
  }
}
*/
/*
void TPZGeoEl1d::SideSubElements(int side, VoidPtrVec &sub) {
        if(!fSubEl) {
        cout << "Error : TPZGeoEl1d::SideSubElements called for an element without subelements\n";
      Print();
        return;
   }
   if(side < 0 || side > 2) {
        cout << "Error : TPZGeoEl1d::SideSubElements called for side " << side << endl;
      Print();
        return;
   }
   if(side < 2) {
   	sub.resize(1);
	   sub[0] = (*fSubEl)[side];
   } else {
   	sub.resize(2);
	   sub[0] = (*fSubEl)[0];
	   sub[1] = (*fSubEl)[1];
	}
}
*/

/**Accumulates the transformation of the jacobian which maps the current
master element space into the space of the master element of the father*/
/*
void TPZGeoEl1d::BuildTransform(int side, TPZGeoEl *father, TPZTransform &t) {
  if(side != 2) return;
  int isub = 0;
  TPZGeoEl *locfather = fFather;
  if(!locfather) {
    cout << "TPZGeoEl1d::BuildTransform could not identify the father element\n";
    return;
  }
  for(; isub<2; isub++) {
    if(locfather->SubElement(isub) == this) break;
  }
  TPZTransform tloc(1);
  REAL store[2];
  TPZFMatrix mult(1,1,store,1);
  TPZFMatrix sum(1,1,store+1,1);
  switch(isub) {
  case 0:
    mult(0,0) = 0.5;
    sum(0,0) = -0.5;
    break;
  case 1:
    mult(0,0) = 0.5;
    sum(0,0) = 0.5;
    break;
  default:
    mult(0,0) = 1.;
    sum(0,0) = 0.;
    cout << "TPZGeoEl1d::BuildTransform subelement not detected within the"
      " father element\n";
  }
  tloc.SetMatrix(mult,sum);
  t = tloc.Multiply(t);
  if(locfather != father) locfather->BuildTransform(side,father,t);
}
*/
TPZTransform TPZGeoEl1d::SideToSideTransform(int sidefrom, int sideto) {
	return TPZShapeLinear::SideToSideTransform(sidefrom,sideto);
}

/*
TPZGeoElSide TPZGeoEl1d::SideSubElement(int side,int subel) {
  if(side < 0 || side > 2) {
    PZError << "TGeoElQ2d::SideSubElement called for side " << side << endl;
    return TPZGeoElSide(0,0);
  }
  if(!fSubEl[0]) return TPZGeoElSide(0,0);

  if(side==2)
    return TPZGeoElSide(fSubEl[subel],2);
  return TPZGeoElSide(fSubEl[side],side);
}
*/
/**Neigh é o lado de um elemento vizinho ao elemento atual El que esta sendo dividido rfnd guarda
   os dois nós globais deste lado no sentido antihorario do elemento El, posições 0 e 1 de rfnd*/
/*
void TPZGeoEl1d::GetSubElement(int side,TPZVec<int> &rfndindex,TPZVec<TPZGeoElSide> &sub) {
  if(!HasSubElement()) {
    sub.Resize(0);
    return;
  }
  switch(side) {
  case 0:
  case 1:
    sub.Resize(1);
    sub[0] = TPZGeoElSide(SubElement(side),side);
    break;
  case 2:
    sub.Resize(2);
    if(SubElement(0)->NodeIndex(0) == rfndindex[0]) {
      sub[0] = TPZGeoElSide(SubElement(0),2);
      sub[1] = TPZGeoElSide(SubElement(1),2);
    } else {
      sub[0] = TPZGeoElSide(SubElement(1),2);
      sub[1] = TPZGeoElSide(SubElement(0),2);
    }
  }
}
*/
/**Inicializa os coeficientes do par de nós do lado I do elemento de referencia*/
/*
void TPZGeoEl1d::SideMasterCo(int side,TPZVec<REAL> &IVec,TPZVec<REAL> &JVec) {
  // modified Philippe 28/7/97
  // if this method is called from a 3d element, IVec and JVec must be zeroed

  IVec.Fill(0.,0);
  JVec.Fill(0.,0);
  if (side==2) {
    IVec[0]=-1.;
    JVec[0]=1.;
  } else if(side==0) {
    PZError << "TPZGeoEl1d::SideMasterCo called for side 0\n";
    IVec[0]=-1.;
  } else if(side==1) {
    PZError << "TPZGeoEl1d::SideMasterCo called for side 1\n";
    IVec[0]=1.;
  }
}

void TPZGeoEl1d::SideMasterCo(int side,TPZFMatrix &coord) {
  int row = coord.Rows();
  if(side == 0 || side == 1) {
    coord.Redim(row,1);
    if(side == 0) coord(0,0) = -1.;
    else coord(0,0) = 1.;
  }
  else {
    coord.Redim(row,2);
    coord(0,0) = -1.;
    coord(0,1) = 1.;
  }
}
*/
/*
TPZGeoElSide TPZGeoEl1d::Father(int side) {
  if(!fFather) return TPZGeoElSide();
  int is;
  for(is=0; is<2; is++) if(fFather->SubElement(is) == this) break;
  if(is> 1) {
    cout << "TPZGeoEl1d::Father is fishy\n";
    return TPZGeoElSide();
  }
  TPZGeoElSide father;
  if(is == side || side == 2) {
    father = TPZGeoElSide(fFather,side);
  }
  return father;
}
*/

TPZGeoEl *TPZGeoEl1d::CreateBCGeoEl(int side, int bc) {
	TPZGeoEl *gel = TPZGeoLinear::CreateBCGeoEl(this,side,bc);
	return gel;
}

void TPZGeoEl1d::LowerDimensionSides(int side,TPZStack<int> &smallsides) {
	int nsidecon = TPZShapeLinear::NSideConnects(side);
	int is;
	for(is=0; is<nsidecon-1; is++) smallsides.Push(TPZShapeLinear::SideConnectLocId(side,is));


}

/*
TPZGeoElSide TPZGeoEl1d::HigherDimensionSides(int side,int targetdimension) {
  if(side > 1 && targetdimension < 1) {
    PZError << "TPZGeoEl1d::HigherDimensionSides called with side = " << side
	    << " targetdimension = " << targetdimension << endl;
    return TPZGeoElSide();
  }
  switch(targetdimension) {
	  case 0 :
	    return TPZGeoElSide(this,side);
	  case 1 :
	    return TPZGeoElSide(this,2);
	  default :
	    return TPZGeoElSide();
  }
}
*/
void TPZGeoEl1d::AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides){
	TPZStack<int> high;
	TPZShapeLinear::HigherDimensionSides(side,high);
	int cap = high.NElements(),s;
	for(s=0; s<cap; s++) {
		if(SideDimension(high[s]) == targetdimension) {
			elsides.Push(TPZGeoElSide(this,high[s]));
		}
	}
}


TPZGeoElSide TPZGeoEl1d::Father2(int side){
	if (!fFather) return TPZGeoElSide();
	int son = WhichSubel();
	return TPZGeoElSide(fFather,TPZRefLinear::FatherSide(side,son));
}

int TPZGeoEl1d::FatherSide(int side, int son){
	return TPZRefLinear::FatherSide(side,son);
}

void TPZGeoEl1d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){
	TPZRefLinear::GetSubElements(this,side,subel);
}

int TPZGeoEl1d::NSideSubElements2(int side) {
	return TPZRefLinear::NSideSubElements(side);
}


TPZTransform TPZGeoEl1d::BuildTransform2(int side, TPZGeoEl * father, TPZTransform &t){//Augusto:09/01/01
	if(side<0 || side>(TPZShapeLinear::NSides-1) || !fFather){
		PZError << "TPZGeoElement::BuildTransform2 side out of range or father null\n";
		return TPZTransform(0,0);
	}
	TPZGeoElSide fathloc = Father2(side);
  int son = WhichSubel();
  TPZTransform trans=fFather->GetTransform(side,son);
  trans = trans.Multiply(t);
  if(fathloc.Element() == father) return trans;
  trans = fFather->BuildTransform2(fathloc.Side(),father,trans);
  return trans;

}

static REAL MidSideNode[3][3] = {
/*00*/{-1.},/*01*/{1.},/*02*/{0.} };

int TPZGeoEl1d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3),pss(3),pf(3),pfs(3);
                    //point son, point side son, point father, point side father : elemento mestre
  pss[1] = 0.; pss[2] = 0.;//1d
  pfs[1] = 0.; pfs[2] = 0.;
  pf[1] = 0.; pf[2] = 0.;
  for(sn=0;sn<2;sn++){
    TPZGeoEl *son = subs[sn];
    for(sd=0;sd<3;sd++){
      ps[0] = MidSideNode[sd][0];//element master point
      ps[1] = MidSideNode[sd][1];// = 0
      ps[2] = MidSideNode[sd][2];// = 0
      if(son->WhichSide(ps) != sd) cout << "Lado nao bate\n";
      TPZTransform telsd = TPZShapeLinear::TransformElementToSide(sd);//2x2
      telsd.Apply(ps,pss);//son element -> side
      son->X(ps,x1);//ponto deformado filho
	  int dim = son->SideDimension(sd);
	  TPZTransform orig(dim);
      TPZTransform t = son->BuildTransform2(sd,gel,orig);
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      telsd = TPZShapeLinear::TransformSideToElement(sdfat);//2x2
      telsd.Apply(pfs,pf);//lado do pai -> pai
      son->Father2(2).Element()->X(pf,x2);
      if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 ){
      	PZError << "\nTransformacao furada\n";
        PZError << "son    = " << (son->Id()) << endl;
        PZError << "father = " << ((son->Father2(2).Element())->Id()) << endl;
        PZError << "side   = " << sd << endl << endl;
        int ok;
        cin >> ok;
      } else {
      	cout << "Transformacao OK!\n";
       	cout << "Filho/lado : " << son->Id() << "/" << sd << endl;
        cout << "Pai : " << son->Father2(2).Element()->Id() << endl << endl;
      }
    }
  }
  return 1;
}

/** Jorge 17/7/99 */
/** Compute the measure of the geometrical element */
/*REAL TPZGeoEl1d::Mesure(int dim) {
  if(dim!=1) return 0.;
  if(fMesure == 0.) {
    TPZGeoNode &nod1 = Mesh()->NodeVec()[fNodeIndexes[0]];
    REAL x0 = nod1.Coord(0), y0 = nod1.Coord(1), z0 = nod1.Coord(2);
    TPZGeoNode &nod2 = Mesh()->NodeVec()[fNodeIndexes[1]];
    x0 -= nod2.Coord(0);
    y0 -= nod2.Coord(1);
    z0 -= nod2.Coord(2);
    for(int i=0;i<3;i++) fMesure = sqrt(x0 * x0 + y0 * y0 + z0 * z0);
  }
  return fMesure;
}
*/
//void TPZGeoEl1d::Center(TPZVec<REAL> &center) {
//  center[0] = 0.;
//}

void TPZGeoEl1d::SetSubElement(int id, TPZGeoEl *el){
  if (id<0 || id >1){
    PZError << "TPZGeoEl1d::Trying do define subelement :" << id << endl;
    return;
  }
  fSubEl[id]=el;
  return;
}
    

TPZIntPoints * TPZGeoEl1d::CreateSideIntegrationRule(int side, int order)
{
	return TPZGeoLinear::CreateSideIntegrationRule(side,order);
}

TPZTransform TPZGeoEl1d::GetTransform(int side,int son) {
	return TPZRefLinear::GetTransform(side,son);
}

void TPZGeoEl1d::CenterPoint(int side, TPZVec<REAL> &masscent){

  if(side < 0 || side > NSides()-1){
    PZError << "TPZGeoEl1d::CenterPoint error side = " << side << endl;
    return;
  }
  TPZShapeLinear::CenterPoint(side,masscent);
}

REAL TPZGeoEl1d::RefElVolume(){
  return TPZShapeLinear::RefElVolume();
}
