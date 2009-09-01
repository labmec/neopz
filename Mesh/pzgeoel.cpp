//Id: $

// -*- c++ -*-
/**File : pzgeoel.c

Contains the methods definition for (abstract) base class TPZGeoEl.
*/
//#include "pztempmat.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include "pzgmesh.h"
#include "pzgnode.h"
#include "pzerror.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pztrnsform.h"
#include "pzvec.h"
#include "pzstack.h"
#include "pzquad.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "TPZRefPattern.h"
//#include "pzgeoelside.h"
//#include "pzshapetetra.h"
//#include "pzshapecube.h"
//#include "pzshapepiram.h"
//#include "pzshapeprism.h"
//#include "pzelg1d.h"
//#include "pzelgt2d.h"
//#include "pzelgq2d.h"
//#include "pzelgt3d.h"
//#include "pzelgpi3d.h"
//#include "pzelgpr3d.h"
//#include "pzelgc3d.h"
#include "pzlog.h"

#include <stdio.h>
#include <stdlib.h>

using namespace std;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzgeoel"));
#endif


class TPZRefPattern;

TPZFMatrix TPZGeoEl::gGlobalAxes;

/**Constructor and destructor*/
TPZGeoEl::~TPZGeoEl(){
  int index = Index();
  fMesh->ElementVec()[index] = 0;
  fMesh->ElementVec().SetFree(index);
};

TPZGeoEl::TPZGeoEl(int id,int materialid,TPZGeoMesh &mesh) {
  fMesh = &mesh;
  fId = id;
  mesh.SetElementIdUsed(id);
  fMatId = materialid;
  this->fReference = NULL;
  fFatherIndex = -1;
  int index = fMesh->ElementVec().AllocateNewElement();
  fIndex = index;
  fMesh->ElementVec()[index] = this;
  this->fNumInterfaces = 0;
//	fMesure = 0.;
}

TPZGeoEl::TPZGeoEl(const TPZGeoEl &el):TPZSaveable(el){
  fMesh = el.fMesh;
  fId = fMesh->CreateUniqueElementId();
  fMatId = el.fMatId;
  this->fReference = NULL;
  fFatherIndex = el.fFatherIndex;
  fIndex = fMesh->ElementVec().AllocateNewElement();
  fMesh->ElementVec()[fIndex] = this;
  this->fNumInterfaces = 0;
}



TPZGeoEl::TPZGeoEl(int materialid,TPZGeoMesh &mesh, int &index) {
  fId = mesh.CreateUniqueElementId();
  fMesh = &mesh;
  fMatId = materialid;
  fReference = NULL;
  fFatherIndex = -1;
  index = fMesh->ElementVec().AllocateNewElement();
  fIndex = index;
  fMesh->ElementVec()[index] = this;
  this->fNumInterfaces = 0;
//	fMesure = 0.;
}
TPZGeoEl::TPZGeoEl(int materialid,TPZGeoMesh &mesh) {
  this->fId = mesh.CreateUniqueElementId();
  this->fMesh = &mesh;
  this->fMatId = materialid;
  this->fReference = NULL;
  this->fFatherIndex = -1;
  const int index = this->Mesh()->ElementVec().AllocateNewElement();
  this->fIndex = index;
  this->Mesh()->ElementVec()[index] = this;
  this->fNumInterfaces = 0;
//	fMesure = 0.;
}

void TPZGeoEl::Initialize(int materialid, TPZGeoMesh &mesh, int &index){
  this->fId = mesh.CreateUniqueElementId();
  this->fMesh = &mesh;
  this->fMatId = materialid;
  this->fReference = NULL;
  this->fFatherIndex = -1;
  index = this->Mesh()->ElementVec().AllocateNewElement();
  this->Mesh()->ElementVec()[index] = this;
  this->fIndex = index;
}

void TPZGeoEl::Shape1d(double x,int num,TPZFMatrix &phi,TPZFMatrix &dphi){
  if(num != 2 && num != 3){
    PZError << "elcalc1d.shape, at this point only linear and quadratic elements\n";
    return;
  }

  if(num == 2) {
    phi(0,0) = (1-x)/2.;
    phi(1,0) = (1+x)/2.;
    dphi(0,0) = -0.5;
    dphi(0,1) = 0.5;
  } else {
    phi(0,0) = -x*(1.-x)*0.5;
    dphi(0,0) = x-0.5;
    phi(1,0) = (1.-x*x);
    dphi(0,1) = -2.*x;
    phi(2,0) = 0.5*x*(1.+x);
    dphi(0,2) = x+0.5;
  }
}

void TPZGeoEl::ShapePhi1d(double x,int num,TPZFMatrix &phi) {
  if(num != 2 && num != 3){
    PZError << "TPZGeoEl ShapePhi1d, at this point only linear and quadratic elements\n";
    return;
  }

  if(num == 2) {
    phi(0,0) = (1-x)/2.;
    phi(1,0) = (1+x)/2.;
  } else {
    phi(0,0) = -x*(1.-x)*0.5;
    phi(1,0) = (1.-x*x);
    phi(2,0) = 0.5*x*(1.+x);
  }
}

int TPZGeoEl::WhichSide(TPZVec<int> &SideNodeIds) {
  int cap = SideNodeIds.NElements();
  int nums = NSides();
  for(int side=0; side<nums; side++) {
    if(NSideNodes(side)==2 && cap == 2) {
      int isn1 = SideNodeIndex(side,0);
      int isn2 = SideNodeIndex(side,1);//sao = para side<3
      if(isn1 == SideNodeIds[0] && isn2 == SideNodeIds[1] ||
      isn2 == SideNodeIds[0] && isn1 == SideNodeIds[1])    return side;
    } else if(NSideNodes(side)== 1 && cap ==1) {
      if(SideNodeIndex(side,0) == SideNodeIds[0]) return side;
      //completar
    } else if(NSideNodes(side) == 3 && cap==3) {
         int sni[3],snx[3],k;
         for(k=0;k<3;k++) snx[k] = SideNodeIndex(side,k);//el atual
         for(k=0;k<3;k++) sni[k] = SideNodeIds[k];//el viz
         for(k=0;k<3;k++) {
           if(snx[0]==sni[k] && snx[1]==sni[(k+1)%3] && snx[2]==sni[(k+2)%3]) return side;
           if(snx[0]==sni[k] && snx[1]==sni[(k+2)%3] && snx[2]==sni[(k+1)%3]) return side;
         }//012 120 201 , 021 102 210
    } else if(NSideNodes(side) == 4 && cap == 4) {//face quadrilateral
         int sni[4],snx[4],k;
         for(k=0;k<4;k++) snx[k] = SideNodeIndex(side,k);//el atual
         for(k=0;k<4;k++) sni[k] = SideNodeIds[k];//vizinho
         if(snx[0]==sni[0]) {
           for(k=1;k<4;k++) {
            if(snx[1]==sni[k] && snx[2]==sni[k%3+1]     && snx[3]==sni[(k+1)%3+1]) return side;
            if(snx[1]==sni[k] && snx[2]==sni[(k+1)%3+1] && snx[3]==sni[k%3+1])     return side;
           }/* 0123 0231 0312 , 0132 0213 0321 */
         } else if(snx[1]==sni[0]) {
           for(k=1;k<4;k++) {
            if(snx[0]==sni[k] && snx[2]==sni[k%3+1]     && snx[3]==sni[(k+1)%3+1]) return side;
            if(snx[0]==sni[k] && snx[2]==sni[(k+1)%3+1] && snx[3]==sni[k%3+1])     return side;
           }/* 0123 0231 0312 , 0132 0213 0321 */                               /* 1023 1230 1302 , 1032 1203 1320 */
         } else if(snx[2]==sni[0]) {
           for(k=0;k<4;k++) {
            if(snx[0]==sni[k] && snx[1]==sni[k%3+1]     && snx[3]==sni[(k+1)%3+1]) return side;
            if(snx[0]==sni[k] && snx[1]==sni[(k+1)%3+1] && snx[3]==sni[k%3+1])     return side;
           }/* 0123 0231 0312 , 0132 0213 0321 */                               /* 2013 2130 2301 , 2031 2103 2310 */
         } else if(snx[3]==sni[0]) {
           for(k=0;k<4;k++) {
            if(snx[0]==sni[k] && snx[1]==sni[k%3+1]     && snx[2]==sni[(k+1)%3+1]) return side;
            if(snx[0]==sni[k] && snx[1]==sni[(k+1)%3+1] && snx[2]==sni[k%3+1])     return side;
           }/* 0123 0231 0312 , 0132 0213 0321 */                               /* 3012 3120 3201 , 3021 3102 3210 */
         }
    } else if(cap<1 || cap > 4) {
      PZError << "TPZGeoEl::WhichSide must be extended\n";
      exit(-1);
    }
  }
  return -1;
}

int TPZGeoEl::NeighbourExists(int side,const TPZGeoElSide &gel) {
  TPZGeoElSide thisside(this,side);
  if(gel == thisside) return 1;
  TPZGeoElSide neighbour = Neighbour(side);
  if(!neighbour.Exists()) return 0;
  while(neighbour != thisside) {
    if(gel == neighbour) return 1;
    neighbour = neighbour.Neighbour();
  }
  return 0;
}


void TPZGeoEl::Print(std::ostream & out) {

  out << "Element index      " << fIndex << endl;
  out << "Element id         " << fId << endl;
  out << "Is GeoBlend ? : ";
  if(this->IsGeoBlendEl())
  {
     out << "true" << endl;
  }
  else
  {
     out << "false" << endl;
  }
  out << "Is Linear Mapping ? : ";
  if(this->IsLinearMapping())
  {
     out << "true" << endl;
  }
  else
  {
     out << "false" << endl;
  }
  out << "Number of nodes    " << NNodes() << endl;
  out << "Corner nodes       " << NCornerNodes() << endl;
  out << "Nodes ids          ";
  int i;
  for (i = 0;i < NNodes();i++) out << NodeIndex(i) << " ";
  out << "\nNumber of sides    " << NSides() << endl;
  if (fMatId < 0) out << "boundary condition " << fMatId << endl;
  else out << "Material id        " << fMatId << endl;
  if (!Father()) out << "no father\n";
  else out << "Father id          " << Father()->Id() << endl;
  if (!SubElement(0)) out << "no subelements";
  else {
    out << "Subelements ids     ";
    for (i = 0;i < NSubElements();i++) {
      if (!SubElement(i) ) continue;
      out << SubElement(i)->Id() << ' ' ;
    }
  }
  out << endl;
  for (i = 0;i < NSides();i++) {
    out << "Neighbours for side   " << i << " : ";
    TPZGeoElSide neighbour = Neighbour(i);
    TPZGeoElSide thisside(this,i);
    if (!neighbour.Exists()) out << "No neighbour\n";
    else {
      int count = 0;
      while (neighbour != thisside && count++ < 10) {
        out << neighbour.Element()->Id() << "/" << neighbour.Side() << ' ';
        neighbour = neighbour.Neighbour();
      }
      out << endl;
    }
  }
  out << "Reference element: " << fReference << endl;
  if (this->Reference()) out << "Reference element index : " << this->Reference()->Index() << endl;
  out << "fNumInterfaces = " << fNumInterfaces << endl;
}

std::ostream &operator<<(std::ostream &out,TPZGeoEl & el) {
  el.Print(out);
  return out;
}

/**The root class TPZGeoEl looks now for a midsidenode*/
/*
TPZGeoNode *TPZGeoEl::MiddleSideNode(int side) {
  //Look for the midsidenode of a neighbour of the same level
  TPZGeoElSide neighbour = Neighbour(side);
  TPZGeoElSide thisside(this,side);
  if(!neighbour.Exists()) return 0;
  while(thisside != neighbour) {
    if(neighbour.HasSubElement())
      return neighbour.Element()->MiddleSideNode(neighbour.Side());
    neighbour = neighbour.Neighbour();
  }
  return NULL;
}
*/
//void TPZGeoEl::GetSubElement(int /*side*/,TPZVec<int> &/*referencenode*/,
//			     TPZVec<TPZGeoElSide> &/*subelements*/) {
//  PZError << "TPZGeoEl::GetSubEl should not be called\n";
//}



//TPZGeoElSide TPZGeoEl::Father(int /*side*/) {
//  PZError << "TPZGeoEl::Father should never be called\n";
//  return TPZGeoElSide();
//}

int TPZGeoEl::Level() {
  int level = 0;
  TPZGeoEl *father = Father();
  while(father) {
    level++;
    father = father->Father();
  }
  return level;
}


int TPZGeoEl::GetTransformId2dQ(TPZVec<int> &idfrom,TPZVec<int> &idto) {

  if(idfrom[0]==idto[0] && idfrom[1]==idto[1]) return 0;//sentido horario     : 0123
  if(idfrom[0]==idto[0] && idfrom[1]==idto[3]) return 1;//sentido antihorario : 0321
//0 esta na posicao 3 e 1 esta na posicao 0  i.e. 1,2,3,0
  if(idfrom[0]==idto[3] && idfrom[1]==idto[0]) return 2;//sentido horario     : 1230
  if(idfrom[0]==idto[1] && idfrom[1]==idto[0]) return 3;//sentido antihorario : 1032

  if(idfrom[0]==idto[2] && idfrom[1]==idto[3]) return 4;//sentido horario     : 2301
  if(idfrom[0]==idto[2] && idfrom[1]==idto[1]) return 5;//sentido antihorario : 2103

  if(idfrom[0]==idto[1] && idfrom[1]==idto[2]) return 6;//sentido horario     : 3012
  if(idfrom[0]==idto[3] && idfrom[1]==idto[2]) return 7;//sentido antihorario : 3210

  return 0;
}

int TPZGeoEl::GetTransformId2dT(TPZVec<int> &idfrom,TPZVec<int> &idto) {
//REVISAR
  if(idto[0]==idfrom[0] && idto[1]==idfrom[1]) return 0;//sentido horario
  if(idto[0]==idfrom[0] && idto[1]==idfrom[2]) return 1;//sentido antihorario

  if(idto[0]==idfrom[1] && idto[1]==idfrom[2]) return 2;//sentido horario
  if(idto[0]==idfrom[1] && idto[1]==idfrom[0]) return 3;//sentido antihorario

  if(idto[0]==idfrom[2] && idto[1]==idfrom[0]) return 4;//sentido horario
  if(idto[0]==idfrom[2] && idto[1]==idfrom[1]) return 5;//sentido antihorario

  return 0;
}

int TPZGeoEl::ElementExists(TPZGeoEl *elem,int id) {

	if(elem == 0 && id == 0) {
      PZError << "\nTPZGeoEl::ElementExists :  element id or element pointer will not be null\n";
      return -1;
   }
	int nelg = fMesh->ElementVec().NElements(),index;
   for(index=0;index<nelg;index++) {
      TPZGeoEl *gel = fMesh->ElementVec()[index];
      if(!gel) continue;
      if(gel == elem || gel->Id() == id) return index;
   }
   return -1;
}

TPZGeoElSide TPZGeoEl::Father2(int /*side*/){//Augusto:09/01/01
  PZError << "TPZGeoEl::Father2 should never be called\n";
  return TPZGeoElSide();
}

int TPZGeoEl::FatherSide(int side, int son){
  PZError << "TPZGeoEl::FatherSide should never be called\n";
  return -1;
}

TPZTransform TPZGeoEl::BuildTransform2(int /*side*/, TPZGeoEl * /*father*/, TPZTransform & /* tr */){//Augusto:09/01/01
  PZError << "TPZGeoEl::BuildTransform2 should never be called\n";
  return TPZTransform(0,0);
}


void TPZGeoEl::GetSubElements2(int /*side*/, TPZStack<TPZGeoElSide> &/*subel*/){//Augusto:09/01/01
  PZError << "TPZGeoEl::GetSubElements2 should never be called\n";
}

void TPZGeoEl::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel, int dimension){
	TPZStack<TPZGeoElSide> subel2;
	GetSubElements2(side,subel2);
	int cap = subel2.NElements();
	int s;
	for(s=0; s<cap; s++) {
		if(subel2[s].Dimension() == dimension) {
			subel.Push(subel2[s]);
		}
	}
}

int TPZGeoEl::WhichSubel(){

  if(!fFatherIndex == -1) {
    PZError << "TPZGeoEl::WhichSubel called with null element\n";
      return -1;
  }
  int son;
  TPZGeoEl *father = Father();
  int nsub = father->NSubElements();
  for(son=0;son<nsub;son++) if(father->SubElement(son) == this) break;
  if(son > (nsub-1)){
    PZError << "TPZGeoEl::WhichSubel son not exist\n";
    return -1;
  }
  return son;
}

int TPZGeoEl::WhichSide(TPZVec<REAL> &pt){
	int dim = Dimension();
	int nums = NSides();//n�mero de lados
	REAL tol = 1.e-06;//toler�ncia do zero -> x � zero se -o<=x<=+o
	int is;
	for(is=0; is<nums; is++) {
		int sdim = SideDimension(is);
		TPZTransform t1 = SideToSideTransform(nums-1,is);
		TPZTransform t2 = SideToSideTransform(is,nums-1);
		TPZVec<REAL> pts(sdim),pt2(dim);
		t1.Apply(pt,pts);
		t2.Apply(pts,pt2);
		REAL dif=0;
		int d;
		for(d=0; d<dim; d++) {
			dif += (pt[d]-pt2[d])*(pt[d]-pt2[d]);
		}
		//dif = sqrt(dif/dim);
		if(dif < tol) return is;
	}
  cout << "TPZGeoEl::WhichSide ERROR : side not found" << endl ;
/*
  if(nums==3){//LINHA
    if((pt[0]<-1.-O) || (pt[0]>1.+O) || pt[1]<-O || pt[1]>O || pt[2]<-O || pt[2]>O) return -1;//est� fora da linha
    if(pt[0]<=-1.+O) return 0;
    if(pt[0]>= 1.-O) return 1;
    return 3;//interior � linha
  } else
  if(nums == 7){//TRI�NGULO
  	int r0=0,r1=0,r2=0;
    if(pt[2]<-O || pt[2]>O) return -1;//est� fora do tri�ngulo
    if( (O<pt[0]) && (pt[0]<1.-O) && (O<pt[1]) && (pt[1]<1.-pt[0]-O) ) return 6;//interior ao tri�ngulo
    REAL ptz = pt[1]-(1.-pt[0]);
    if( (-O<=pt[0]) && (pt[0]<=1.+O) && (-O<=pt[1]) && (pt[1]<=O   ) ) r0 = 1;
    if( (-O<=pt[0]) && (pt[0]<=1.+O) && (-O<=ptz  ) && (ptz  <=O   ) ) r1 = 1;
    if( (-O<=pt[0]) && (pt[0]<=   O) && (-O<=pt[1]) && (pt[1]<=1.+O) ) r2 = 1;
    if(r0){
    	 if(r1) return 1;//canto c1=1
      if(r2) return 0;//canto c0=0
      return 3;//aresta r0=3
    }
    if(r1){
    	 if(r2) return 2;//canto c2=2
      return 4;//aresta r1=4
    }
    if(r2) return 5;//aresta r2=5
    return -1;//esta fora do tri�ngulo
  } else
  if(nums==9){//QUADRIL�TERO
  	 int r0=0,r1=0,r2=0,r3=0;
    if(pt[2]<-O || pt[2]>O) return -1;//est� fora do quadril�tero
    if( (-1.+O< pt[0]) && (pt[0]<  1.-O) && (-1.+O< pt[1]) && (pt[1]<  1.-O) ) return 9;
    if( (-1.-O<=pt[0]) && (pt[0]<= 1.+O) && (-1.-O<=pt[1]) && (pt[1]<=-1.+O) ) r0 = 1;
    if( ( 1.-O<=pt[0]) && (pt[0]<= 1.+O) && (-1.-O<=pt[1]) && (pt[1]<= 1.+O) ) r1 = 1;
    if( (-1.-O<=pt[0]) && (pt[0]<= 1.+O) && ( 1.-O<=pt[1]) && (pt[1]<= 1.+O) ) r2 = 1;
    if( (-1.-O<=pt[0]) && (pt[0]<=-1.+O) && (-1.-O<=pt[1]) && (pt[1]<= 1.+O) ) r3 = 1;
    if(r0){
      if(r3) return 0;//canto c0=0
    	 if(r1) return 1;//canto c1=1
      return 4;//aresta r0=4
    }
    if(r1){
    	 if(r2) return 2;//canto c2=2
      return 5;//aresta r1=5
    }
    if(r2){
      if(r3) return 3;//canto c3=3
      return 6;//aresta r2=6
    }
    if(r3) return 7;//aresta r3=7
    return -1;//esta fora do quadril�tero
  } else
  if(nums==15){//TETRAEDRO
    int f0=0,f1=0,f2=0,f3=0;//faces
    if(-O<=pt[2] && pt[2]<=O && -O<=pt[0] && pt[0]<=1.+O && -O<=pt[1] && pt[1]<=1.-pt[0]+O) f0 = 1;
    if(-O<=pt[1] && pt[1]<=O && -O<=pt[0] && pt[0]<=1.+O && -O<=pt[2] && pt[2]<=1.-pt[0]+O) f1 = 1;
    if(-O<=pt[0] && pt[0]<=O && -O<=pt[1] && pt[1]<=1.+O && -O<=pt[2] && pt[2]<=1.-pt[1]+O) f3 = 1;
    if(-O<=pt[0] && pt[0]<=1.+O && -O<=pt[1] && (pt[1]<=1.-pt[0]+O) && (1.-pt[0]-pt[1]-O<=pt[2]) && (pt[2]<=1.-pt[0]-pt[1]+O)) f2 = 1;
    if(f0){
      if(f1){
        if(f2) return 1;
        if(f3) return 0;
        return 4;
      }
      if(f2 && f3) return 2;
      if(f2) return 5;
      if(f3) return 6;
      return 10;
    }
    if(f1){
      if(f2 && f3) return 3;
      if(f2) return 8;
      if(f3) return 7;
      return 11;
    }
    if(f2){
      if(f3) return 9;
      return 12;
    }
    if(f3) return 13;
    if(O<pt[0] && pt[0]<1.-O && O<pt[1] && (pt[1]<1.-pt[0]-O) && O<pt[2] && (pt[2]<1.-pt[0]-pt[1]-O)) return 14;
    return -1;
  } else
  if(nums==27){//HEXAEDRO
    int f0=0,f1=0,f2=0,f3=0,f4=0,f5=0;
    if(-1.-O<=pt[2] && pt[2]<=-1.+O && -1.-O<=pt[0] && pt[0]<=1.+O && -1.-O<=pt[1] && pt[1]<=1.+O) f0 = 1;
    if(-1.-O<=pt[1] && pt[1]<=-1.+O && -1.-O<=pt[0] && pt[0]<=1.+O && -1.-O<=pt[2] && pt[2]<=1.+O) f1 = 1;
    if( 1.-O<=pt[0] && pt[0]<= 1.+O && -1.-O<=pt[1] && pt[1]<=1.+O && -1.-O<=pt[2] && pt[2]<=1.+O) f2 = 1;
    if( 1.-O<=pt[1] && pt[1]<= 1.+O && -1.-O<=pt[0] && pt[0]<=1.+O && -1.-O<=pt[2] && pt[2]<=1.+O) f3 = 1;
    if(-1.-O<=pt[0] && pt[0]<=-1.+O && -1.-O<=pt[1] && pt[1]<=1.+O && -1.-O<=pt[2] && pt[2]<=1.+O) f4 = 1;
    if( 1.-O<=pt[2] && pt[2]<= 1.+O && -1.-O<=pt[0] && pt[0]<=1.+O && -1.-O<=pt[1] && pt[1]<=1.+O) f5 = 1;
    if(f0){
      if(f1){
        if(f2) return 1;
        if(f4) return 0;
        return 8;
      }
      if(f3){
        if(f2) return 2;
        if(f4) return 3;
        return 10;
      }
      if(f2) return 9;
      if(f4) return 11;
      return 20;
    }
    if(f5){
      if(f1){
        if(f2) return 5;
        if(f4) return 4;
        return 16;
      }
      if(f3){
        if(f2) return 6;
        if(f4) return 7;
        return 18;
      }
      if(f2) return 17;
      if(f4) return 19;
      return 25;
    }
    if(f1){
      if(f2) return 13;
      if(f4) return 12;
      return 21;
    }
    if(f3){
      if(f2) return 14;
      if(f4) return 15;
      return 23;
    }
    if(f2) return 22;
    if(f4) return 24;
    if((-1.+O<pt[0]) && (pt[0]<1.-O) && (-1.+O<pt[1]) && (pt[1]<1.-O) && (-1.+O<pt[2]) && (pt[2]<1.-O)) return 26;
    return -1;
  } else
  if(nums==19){//PIR�MIDE
    int f0=0,f1=0,f2=0,f3=0,f4=0;
    if(-O<=pt[2] && pt[2]<=O && (-1.-O<=pt[0]) && (pt[0]<=1.+O) && (-1.-O<=pt[1]) && (pt[1]<=1.+O)) f0 = 1;
    if((-1.-O<=pt[0]) && (pt[0]<=1.+O) && (-1.-O<=pt[1]) && (pt[1]<=-fabs(pt[0])+O) && (1.+pt[1]-O<=pt[2]) && (pt[2]<=1.+pt[1]+O)) f1 = 1;
    if((fabs(pt[1])-O<=pt[0]) && (pt[0]<=1.+O) && (-1.-O<=pt[1]) && (pt[1]<=1.+O) && (1.-pt[0]-O<=pt[2]) && (pt[2]<=1.+pt[0]+O)) f2 = 1;
    if((-1.-O<=pt[0]) && (pt[0]<=1.+O) && (fabs(pt[0])-O<=pt[1]) && (pt[1]<=1.+O) && (1.-pt[1]-O<=pt[2]) && (pt[2]<=1.-pt[1]+O)) f3 = 1;
    if((-1.-O<=pt[0]) && (pt[0]<=-fabs(pt[1])+O) && (-1.-O<=pt[1]) && (pt[1]<=1.+O) && (1.+pt[0]-O<=pt[2]) && (pt[2]<=1.+pt[0]+O)) f4 = 1;
    if(f0){
      if(f1){
        if(f4) return 0;
        if(f2) return 1;
        return 5;
      }
      if(f3){
        if(f4) return 3;
        if(f2) return 2;
        return 7;
      }
      if(f2) return 6;
      if(f4) return 8;
      return 13;
    }
    if(f1){
      if(f2 && f4) return 4;
      if(f2) return 10;
      if(f4) return 9;
      return 14;
    }
    if(f3){
      if(f2) return 11;
      if(f4) return 12;
      return 16;
    }
    if(f2) return 15;
    if(f4) return 17;
    if((-1.+O<pt[0]) && (pt[0]<1.-O) && (-1.+O<pt[1]) && (pt[1]<=-fabs(pt[0])) && O<pt[2] && (pt[2]<1.+pt[1]-O)) return 18;
    if((fabs(pt[1])<=pt[0]) && (pt[0]<1.-O) && (-1.+O<pt[1]) && (pt[1]<=1.-O) && O<pt[2] && (pt[2]<1.-pt[0]-O)) return 18;
    if((-1.+O<pt[0]) && (pt[0]<1.-O) && (fabs(pt[0])<=pt[1]) && (pt[1]<1.-O) && O<pt[2] && (pt[2]<1.-pt[1]-O)) return 18;
    if((-1.+O<pt[0]) && (pt[0]<=-fabs(pt[1])) && (-1.+O<pt[1]) && (pt[1]<=1.-O) && O<pt[2] && (pt[2]<1.+pt[0]-O)) return 18;
    return -1;
  } else
  if(nums==21){//PRISMA
    int f0=0,f1=0,f2=0,f3=0,f4=0;
    if((-1.-O<=pt[2]) && (pt[2]<=-1.+O) && -O<=pt[0] && (pt[0]<=1.+O) && -O<=pt[1] && (pt[1]<=1.-pt[0]+O)) f0 = 1;
    if((-O<=pt[0]) && (pt[0]<=1.+O) && -O<=pt[1] && (pt[1]<=O) && (-1.-O<=pt[2]) && (pt[2]<=1.+O)) f1 = 1;
    if((-O<=pt[0]) && (pt[0]<=1.+O) && (1.-pt[0]-O<=pt[1]) && (pt[1]<=1.-pt[0]+O) && (-1.-O<=pt[2]) && (pt[2]<=1.+O)) f2 = 1;
    if((-O<=pt[0]) && pt[0]<=O && -O<=pt[1] && (pt[1]<=1.+O) && (-1.-O<=pt[2]) && (pt[2]<=1.+O)) f3 = 1;
    if((1.-O<=pt[2]) && (pt[2]<=1.+O) && -O<=pt[0] && (pt[0]<=1.+O) && -O<=pt[1] && (pt[1]<=1.-pt[0]+O)) f4 = 1;
    if(f0){
      if(f1){
        if(f2) return 1;
        if(f3) return 0;
        return 6;
      }
      if(f2){
        if(f3) return 2;
        return 7;
      }
      if(f3) return 8;
      return 15;
    }
    if(f1){
      if(f4){
        if(f2) return 4;
        if(f3) return 3;
        return 12;
      }
      if(f2) return 10;
      if(f3) return 9;
      return 16;
    }
    if(f4){
      if(f2){
        if(f3) return 5;
        return 13;
      }
      if(f3) return 14;
      return 19;
    }
    if(f2){
      if(f3) return 11;
      return 17;
    }
    if(f3) return 18;
    if((-1.+O<pt[2]) && (pt[2]<1.-O) && O<pt[0] && (pt[0]<1.-O) && O<pt[1] && (pt[1]<1.-pt[0]-O)) return 20;
    return -1;
  }
*/
	 return -1;
}

void TPZGeoEl::CheckSubelDataStructure(){

	TPZVec<TPZGeoEl *> sub;
  Divide(sub);
  int side,nsubside,isub,ok;            //no m�ximo 9 subelementos/lados
  TPZStack<TPZGeoElSide> subside;//4 faces+4 arestas+1 canto
  ofstream out("Checksubdata.dat");
  for(side=0;side<NSides();side++){//percorre-se os lados do elemento atual = pai
    GetSubElements2(side,subside);
    nsubside = subside.NElements();
    for(isub=0;isub<nsubside;isub++){
      TPZGeoElSide neigh = subside[isub];//1o o subelemento depois os seus vizinhos
      while(neigh.Element()){//existe vizinho!
        TPZGeoElSide fatside = neigh.Element()->Father2(neigh.Side());
        if(!fatside.Element()){
          cout << "Inconsistencia de dados : \n";
          cout << "Pai atual/lado : " << Id() << "/" << side << endl;
          cout << "Pai do sub elemento e' nulo!\n";
          cout << "Sub elemento/lado : " << neigh.Element()->Id() << "/" << neigh.Side() << endl;
          cin >> ok;
        } else //o pai existe � igual ao atual e a dimens�o dos lados � a mesma OK!
          if(fatside.Element() == this){
            cout << "Pai atual/lado        : " << Id() << "/" << side << endl;
            cout << "Sub elemento/lado     : " << neigh.Element()->Id() << "/" << neigh.Side() << endl;
            cout << "Pai sub elemento/lado : " << fatside.Element()->Id() << "/" << fatside.Side() << endl;
             out << "Pai atual/lado        : " << Id() << "/" << side << endl;
             out << "Sub elemento/lado     : " << neigh.Element()->Id() << "/" << neigh.Side() << endl;
             out << "Pai sub elemento/lado : " << fatside.Element()->Id() << "/" << fatside.Side() << endl;
          if(fatside.Side()!=side){
            cout << "Dados acima inconsistentes : lados distintos\n";
             out << "Dados acima inconsistentes : lados distintos\n";
            cin >> ok;
          }
        } else {
          cout << "OK!: Vizinho com pai diferente\n";
        }
        neigh = neigh.Neighbour();
        if(!neigh.Element()){
        	cout << "Vizinho nulo\n\n";
           out << "Vizinho nulo\n\n";
          break;
        }
        cout << "\n";
         out << "\n";
        if(neigh.Element()==subside[isub].Element()) break;
        //no pr�ximo while confia-se que o ciclo de conectividades � fechado e os viz sao filhos de alguem
        TPZGeoEl *gel = neigh.Element()->Father2(neigh.Element()->NSides()-1).Element();//pai pelo interior:sempre existe
        while(gel && !(gel== this)) neigh = neigh.Neighbour();
        if(neigh.Element()==subside[isub].Element()){
          cout << "Nao existe vizinho irmao!\n";
          cout << "Sub elemento/lado " << neigh.Element()->Id() << "/" << neigh.Side() << endl;
          cout << "OK!\n";
        }
      }
    }
  }
  out.flush();
  out.close();
}

/*
int TPZGeoEl::main(TPZGeoEl *gel,int type){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd,i;
  REAL x1store[3],x2store[3];
  TPZManVector<REAL> x1(3,x1store,2),x2(3,x2store,2);//x1 no filho deformado, x2 no pai deformado
  REAL ptstore0[3],ptstore1[3],ptstore2[3],ptstore3[3];
  TPZManVector<REAL> ps(3,ptstore0,3),pss(3,ptstore1,3),pf(3,ptstore2,3),pfs(3,ptstore3,3);
                    //point son, point side son, point father, point side father : elemento mestre
  pss[1] = 0.; pss[2] = 0.;//1d e 2d
  pfs[1] = 0.; pfs[2] = 0.;
  pf[1] = 0.; pf[2] = 0.;
  for(sn = 0; sn < NSubElements(); sn++){
    TPZGeoEl *son = subs[sn];
    for(sd=0; sd < NSides(); sd++){
      if(type==2) for(i=0;i<3;i++) ps[i] = TPZGeoEl1d::MidSideNode[sd][i];
      if(type==3) for(i=0;i<3;i++) ps[i] = TPZGeoElT2d::MidSideNode[sd][i];
      if(type==4) for(i=0;i<3;i++) ps[i] = TPZGeoElQ2d::MidSideNode[sd][i];
      if(type==7) for(i=0;i<3;i++) ps[i] = TPZGeoElT3d::MidSideNode[sd][i];
      if(type==5) for(i=0;i<3;i++) ps[i] = TPZGeoElPi3d::MidSideNode[sd][i];
      if(type==6) for(i=0;i<3;i++) ps[i] = TPZGeoElPr3d::MidSideNode[sd][i];
      if(type==8) for(i=0;i<3;i++) ps[i] = TPZGeoElC3d::MidSideNode[sd][i];
      TPZTransform telsd(0,0);
      if(type==2) telsd = TPZShapeLinear::TransformElementToSide(sd);//1x1
      if(type==3) telsd = TPZShapeTriang::TransformElementToSide(sd);//2x2
      if(type==4) telsd = TPZShapeQuad::TransformElementToSide(sd);//2x2
      if(type==7) telsd = TPZShapeTetra::TransformElementToSide(sd);//3x3
      if(type==5) telsd = TPZShapePiram::TransformElementToSide(sd);//3x3
      if(type==6) telsd = TPZShapePrism::TransformElementToSide(sd);//3x3
      if(type==8) telsd = TPZShapeCube::TransformElementToSide(sd);//3x3
      telsd.Apply(ps,pss);//son element -> side
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform(son->SideDimension(sd);
      TPZTransform t = son->BuildTransform2(sd,gel,trans);
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      if(type==2) telsd = TPZShapeLinear::TransformSideToElement(sdfat);//1x1
      if(type==3) telsd = TPZShapeTriang::TransformSideToElement(sdfat);//2x2
      if(type==4) telsd = TPZShapeQuad::TransformSideToElement(sdfat);//2x2
      if(type==7) telsd = TPZShapeTetra::TransformSideToElement(sdfat);//3x3
      if(type==5) telsd = TPZShapePiram::TransformSideToElement(sdfat);//3x3
      if(type==6) telsd = TPZShapePrism::TransformSideToElement(sdfat);//3x3
      if(type==8) telsd = TPZShapeCube::TransformSideToElement(sdfat);//3x3
      telsd.Apply(pfs,pf);//lado do pai -> pai
      son->Father2(son->NSides()-1).Element()->X(pf,x2);
      if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) ) > 1.e-10 ){
      	PZError << "\nTransformacao furada\n";
        PZError << "son    = " << (son->Id()) << endl;
        PZError << "father = " << ((son->Father2(son->NSides()-1).Element())->Id()) << endl;
        PZError << "side   = " << sd << endl << endl;
        int ok;
        cin >> ok;
      } else {
      	cout << "Transformacao OK!\n";
       	cout << "Filho/lado : " << son->Id() << "/" << sd << endl;
        cout << "Pai : " << son->Father2(son->NSides()-1).Element()->Id() << endl << endl;
      }
    }
  }
  return 1;
}
*/
/**Initializes the external connectivities of the subelements*/
void TPZGeoEl::SetSubElementConnectivities() {

//  this->Print(cout);
  int side;//,el;
  for(side=0; side<NCornerNodes(); side++) {
    TPZGeoElSide thisside(this,side);
    TPZStack<TPZGeoElSide> subel;
    this->GetSubElements2(side,subel);
    //Isso ocorre com padroes de refinamento direcional... Cesar 23-02-2004
    //if(subel.NElements() != 1) {
    //  cout << "TPZGeoEl::SetSubElementConnectivities wrong data\n";
    //} else {
    //for (el=0;el<subel.NElements();el++)
      if(!subel[0].NeighbourExists(thisside)) subel[0].SetConnectivity(thisside);
    //}
  }
  for(side=NCornerNodes(); side<NSides(); side++) {
	  TPZGeoElSide thisside(this,side);
	  TPZGeoElSide neighbour = this->Neighbour(side);
	  while(neighbour.Exists() && neighbour != thisside) {
		  if(neighbour.HasSubElement() && neighbour.NSubElements2() != 1) {
			  TPZStack<TPZGeoElSide> elvec,neighvec;
			  GetSubElements2(side,elvec);
			  neighbour.GetSubElements2(neighvec);
        // The currently divided element may have only one element as a son. In this case, the son is already
        // neighbour of the father
			  if(elvec.NElements() > 1) TPZGeoElSide::BuildConnectivities(elvec,neighvec);
			  break;
		  }
		  neighbour = neighbour.Neighbour();
	  }
  }
  InitializeNeighbours();
}
/*
void TPZGeoEl::ComputeXInverse(TPZVec<REAL> &XD, TPZVec<REAL> &ksi){

  TPZManVector<REAL,3> X0(3,0.);
   if(ksi.NElements()!=Dimension()) {
     PZError << "\nTPZGeoEl::ComputeXInverse vector dimension error\n";
     ksi.Resize(Dimension(),0.);//zero esta em todos os elementos mestres
    //return;
   }
  X(ksi,X0);//ksi deve ter dimensao do elemento atual
  TPZFNMatrix<3> DelX(3,1);
  int i;
  for(i=0; i<3; i++) DelX(i,0) = XD[i]-X0[i];
  int dim = Dimension();
  TPZFNMatrix<3> residual(dim,1),delksi(dim,1);
  REAL detJ;
  TPZFNMatrix<9> J(dim,dim,0.),axes(3,3,0.),Inv(dim,dim,0.);
  TPZFNMatrix<9> JXt(dim,3,0.),JX(3,dim,0.),JXtJX(dim,dim,0.);
  int nao = 0;
  if(NSides() == 19 && nao){
	  ksi.Resize(3,0.);
	  REAL epsilon = 0.002;
	  ofstream outp("JACOBIANO");
	  for(int l=0;l<10;l++){
		 ksi[0] = l*epsilon;
	 	 ksi[1] = l*epsilon/2.0;
		 outp << "ksi : " << ksi[0] << "  "	<< ksi[1] << endl;
	     Jacobian(ksi,J,axes,detJ,Inv);
		 outp << "\nJacobiano ";
		 J.Print("",outp);
		 outp << "\nInversa do Jacobiano ";
		 Inv.Print("",outp);
		 TPZFMatrix Unit = J*Inv;
		 outp << "Jac * InvJac = ";
		 Unit.Print("",outp);
		 outp << "\n det Jacobiano : " << detJ << "\n";
	  }
	  outp.flush();
	  //outp.close();
	  //exit(-1);
  }
  Jacobian(ksi,J,axes,detJ,Inv);
  TPZFNMatrix<9> axest;
  axes.Transpose(&axest);
  axest.Resize(3,dim);//casos 1D e 2D onde JX espacial � 1x3 e 2x3 respectivamente
  if(dim==1){
         JX(0,0) = axest(0,0)*J(0,0);
         JX(1,0) = axest(1,0)*J(0,0);
        JX(2,0) = axest(2,0)*J(0,0);
  } else {
       axest.Multiply(J,JX,0,1);
  }
  JX.Transpose(&JXt);
  JXt.Multiply(JX,JXtJX,0,1);//JXtJX = JXt*JX;
  JXt.Multiply(DelX,residual);//cout << "\nComputeXInverse: : \n";
  JXtJX.SolveDirect(residual,ELU);//cout << "Atual/dimensao : " << Id() << " / " << Dimension();
  for(i=0; i<dim; i++) ksi[i] += residual(i,0);
  X(ksi,X0);
  for(i=0; i<3; i++) DelX(i,0) = XD[i]-X0[i];
  REAL error = Norm(DelX);
  if(error > 1.0E-4) {
	  cout << "ComputeXInverse did not compute the inverse correctly\n";
          Jacobian(ksi,J,axes,detJ,Inv);
  }

}
*/
void TPZGeoEl::ComputeXInverse(TPZVec<REAL> &XD, TPZVec<REAL> &ksi, double Tol){

  REAL error = 10.;
  int iter = 0;
  const int nMaxIter = 1000;
  while(error > Tol && iter < nMaxIter)
  {
      iter++;
      TPZManVector<REAL,3> X0(3,0.);
      if(ksi.NElements()!=Dimension()) {
        PZError << "\nTPZGeoEl::ComputeXInverse vector dimension error\n";
        ksi.Resize(Dimension(),0.);//zero esta em todos os elementos mestres
        //return;
      }
      X(ksi,X0);//ksi deve ter dimensao do elemento atual
      TPZFNMatrix<9> DelX(3,1);
      int i;
      for(i=0; i<3; i++) DelX(i,0) = XD[i]-X0[i];
      int dim = Dimension();
      TPZFNMatrix<9> residual(dim,1),delksi(dim,1);
      REAL detJ;
      TPZFNMatrix<9> J(dim,dim,0.),axes(3,3,0.),Inv(dim,dim,0.);
      TPZFNMatrix<9> JXt(dim,3,0.),JX(3,dim,0.),JXtJX(dim,dim,0.);
      int nao = 0;
      if(NSides() == 19 && nao){
              ksi.Resize(3,0.);
              REAL epsilon = 0.002;
              ofstream outp("JACOBIANO");
              for(int l=0;l<10;l++){
                    ksi[0] = l*epsilon;
                    ksi[1] = l*epsilon/2.0;
                    outp << "ksi : " << ksi[0] << "  "	<< ksi[1] << endl;
                Jacobian(ksi,J,axes,detJ,Inv);
                    outp << "\nJacobiano ";
                    J.Print("",outp);
                    outp << "\nInversa do Jacobiano ";
                    Inv.Print("",outp);
                    TPZFMatrix Unit = J*Inv;
                    outp << "Jac * InvJac = ";
                    Unit.Print("",outp);
                    outp << "\n det Jacobiano : " << detJ << "\n";
              }
              outp.flush();
              //outp.close();
              //exit(-1);
      }
      Jacobian(ksi,J,axes,detJ,Inv);
      TPZFNMatrix<9> axest;
      axes.Transpose(&axest);
      axest.Resize(3,dim);//casos 1D e 2D onde JX espacial � 1x3 e 2x3 respectivamente
      if(dim==1){
            JX(0,0) = axest(0,0)*J(0,0);
            JX(1,0) = axest(1,0)*J(0,0);
            JX(2,0) = axest(2,0)*J(0,0);
      } else {
          axest.Multiply(J,JX,0,1);
      }
      JX.Transpose(&JXt);
      JXt.Multiply(JX,JXtJX,0,1);//JXtJX = JXt*JX;
      JXt.Multiply(DelX,residual);//cout << "\nComputeXInverse: : \n";
      JXtJX.SolveDirect(residual,ELU);//cout << "Atual/dimensao : " << Id() << " / " << Dimension();
      for(i=0; i<dim; i++) ksi[i] += residual(i,0);
      X(ksi,X0);
      for(i=0; i<3; i++) DelX(i,0) = XD[i]-X0[i];
      error = Norm(DelX);
  }

#ifdef DEBUG
  if(iter == nMaxIter){
    std::stringstream sout;
    sout << "Error at " << __PRETTY_FUNCTION__ << " - nMaxIter was reached before tolerance is achieved";
    PZError << "\n" << sout.str() << "\n";
#ifdef LOG4CXX
    LOGPZ_ERROR(logger,sout.str().c_str());
#endif
    DebugStop();
  }///if

  if(this->IsInParametricDomain(ksi,Tol) == false){
    std::stringstream sout;
    sout << "Error at " << __PRETTY_FUNCTION__ << " - Ksi parameter found is outside parametric element domain. Element type = " << this->TypeName();
    sout << " - ksi = ";
    for(int i = 0; i < ksi.NElements(); i++) sout << ksi[i] << "\t";
    sout << "iter = " << iter << ", nMaxIter = " << nMaxIter << "\t";
    sout << "X = ";
    for(int i = 0; i < XD.NElements(); i++) sout << XD[i] << "\t";
    sout << "\n";this->Print(sout);
    PZError << "\n" << sout.str() << "\n";
#ifdef LOG4CXX
    LOGPZ_ERROR(logger,sout.str().c_str());
#endif
    DebugStop();
  }///if

#endif

}

TPZTransform TPZGeoEl::ComputeParamTrans(TPZGeoEl *fat,int fatside, int sideson){

  //transformacao do lado de elemento pequeno para elemento grande que o contem
  int dimf = fat->Dimension();
  int dimsf = fat->SideDimension(fatside);
  int dim = Dimension();
  int dimss = SideDimension(sideson);
  if(dimsf < dimss){
    PZError << "\nTPZGeoEl::ComputeParamTrans called with sides error\n";
    exit(-1);
  }

  /**para o canto do pai n�o existe transformac�o definida*/
  if(!fat->SideDimension(fatside)) return TPZTransform(0,0);

  REAL weight;
  TPZFNMatrix<9> jac(dim,dim),axes(3,3,0.);
  TPZFNMatrix<9> jacinv(dim,dim);
  TPZManVector<REAL,3> x(3,0.);
  TPZManVector<REAL,3> intpoint(dimss,0.);
  int tam = (dimss+1);
  TPZFNMatrix<16> hess(tam,tam,0.),grad0(tam,1,0.);
  TPZIntPoints *intrule = CreateSideIntegrationRule(sideson,2);
  TPZManVector<int,2> order(dimss,2);
  intrule->SetOrder(order);
  //integra��o sobre o lado-filho contido no lado-pai
  int ij,ik,indp;
  REAL D2Edaikdaij,D2Edcidaij,D2Edci2;
  D2Edcidaij = 0.;
  D2Edci2 = 0.;
  for(ij=0;ij<dimss;ij++){
    for(ik=ij;ik<dimss;ik++){
      D2Edaikdaij = 0.;
      D2Edcidaij  = 0.;
      for(indp = 0; indp < intrule->NPoints(); ++indp){
        intrule->Point(indp,intpoint,weight);
        //Jacobian(intpoint,jac,axes,detjac,jacinv);/**estes passos n�o s�o precisos pois:*/
        //weight *= fabs(detjac);/**ambas as matrizes, A e b, multiplicam a mesma constante detjac*/
        D2Edaikdaij += intpoint[ik]*intpoint[ij]*weight;
        if(ik==ij) D2Edcidaij += intpoint[ij]*weight;
        if(ij==0 && ik==0) D2Edci2 += weight;
      }
      hess(ij,ik) = 2.*D2Edaikdaij;
      hess(ik,ij) = hess(ij,ik);/**basta repassar sendo ik>ij*/
      if(ik==ij) {
      	hess(ij,dimss) = 2.*D2Edcidaij;
        hess(dimss,ij) = hess(ij,dimss);
      }
      if(ij==0 && ik==0) hess(dimss,dimss) = 2.*D2Edci2;
    }
  }// final do integral hess
  //do lado sideson para o elemento atual (filho)
  TPZTransform tsidetoson(Dimension());//identidade
  if(dimss<Dimension()) tsidetoson = SideToSideTransform(sideson,NSides()-1);
  TPZTransform fatelside = fat->SideToSideTransform(fat->NSides()-1,fatside);
  TPZManVector<REAL,3> sidepoint(Dimension());//dimensao do dominio da transformacao X do filho
  int j;//transf. para o lado do pai
  TPZFNMatrix<9> A(dimsf,dimss,0.),sol(dimsf,1,0.);
  for(int ifat=0;ifat<dimsf;ifat++){//numero de variaveis do pai
    REAL DEdci = 0.;
    for(j=0;j<(dimss+1);j++){
      REAL DEdaij = 0.;
      for(indp = 0; indp < intrule->NPoints(); ++indp){
        intrule->Point(indp,intpoint,weight);
	tsidetoson.Apply(intpoint,sidepoint);//lado do filho para o seu interior: mestre
        X(sidepoint,x);//ponto do mestre do filho para o filho deformado, 3 coordenadas
        TPZVec<REAL> csi(dimf,0.);/**o seguinte passo n�o � preciso dado que fat = elemento mestre*/
        fat->ComputeXInverse(x,csi);//csi ponto no pai mestre
        TPZVec<REAL> outcsi(dimsf,0.);
		fatelside.Apply(csi,outcsi);
//        fat->ProjectPointElToSide(fatside,csi,outcsi);
        if(j<dimss) DEdaij += outcsi[ifat]*intpoint[j]*weight;
        if(j==0) DEdci += outcsi[ifat]*weight;
      }
      if(j<dimss) grad0(j,0) = 2.*DEdaij;
      if(j==0) grad0(dimss,0) = 2.*DEdci;
      if(!dimss) grad0(0,0) = DEdci;
    }//final integral gradiente
    //resolu��o do sistema para cada variavel ifat do pai
    if(dimss) hess.SolveDirect(grad0,ELU);
    for(int k=0;k<dimss;k++) A(ifat,k) = grad0(k,0);
    sol(ifat,0) = grad0(dimss,0);
  }//fim sistema ifat
  delete intrule;
  TPZTransform t(dimsf,dimss);
  t.SetMatrix(A,sol);
  return t;
}

REAL TPZGeoEl::Distance(TPZVec<REAL> &centel,TPZVec<REAL> &centface){

  if(centel.NElements() != 3 || centface.NElements() != 3)
    PZError << "TPZGeoEl::Distance point dimension error\n";

  int i;
  REAL distance = 0.;
  for(i=0;i<3;i++)
    distance += (centel[i]-centface[i])*(centel[i]-centface[i]);

  return sqrt(distance);
}

REAL TPZGeoEl::ElementRadius(){

   switch( this->Dimension() )
   {
      case 0:
	 return 0.;
 //	 break;

      case 1:
      case 2:{
	 TPZManVector<REAL, 3> centel  (this->Dimension(), 0.);
	 TPZManVector<REAL, 3> centface(this->Dimension(), 0.);
	 TPZManVector<REAL, 3> masscent(3,0.), massface(3, 0.);
	 REAL mindist = 1000.;
	 REAL dist;

	 CenterPoint(NSides()-1,centel);
	 X(centel,masscent);

	 int nsides = NSides();
	 for(int iside = 0; iside < nsides - 1; iside++){
	    CenterPoint(iside, centface);
	    X(centface,massface);
	    dist = Distance(masscent,massface);
	    if(mindist > dist) mindist = dist;
	 }
	 return mindist;
//	 break;
      }

      case 3:{
	 TPZManVector<REAL, 3> centel(3,0.),centface(3,0.);
	 TPZManVector<REAL, 3> masscent(3,0.),massface(3,0.);
	 REAL mindist = 1000.,dist;
	 CenterPoint(NSides()-1,centel);
	 X(centel,masscent);
	 int nfaces,face,face0;
	 if(NSides()==15) {nfaces = 14; face0 = 10;}//tetrahedron
	 else if(NSides()==19) {nfaces = 18; face0 = 13;}//pyramid
	 else if(NSides()==21) {nfaces = 20; face0 = 15;}//prism
	 else if(NSides()==27) {nfaces = 26; face0 = 20;}//hexahedro
	 else return 0.;//line, point, triangle, quadrilateral
	 for(face=face0;face<nfaces;face++){
	    CenterPoint(face,centface);
	    X(centface,massface);
	    dist = Distance(masscent,massface);
	    if(mindist > dist) mindist = dist;
	 }
	 return mindist;
//	 break;
      }

      default:
	 PZError <<  "TPZGeoEl::ElementRadius - Dimension not implemented." << endl;
	 return 0.;

   }//end of switch
//   return 0.;
} //end of method

REAL TPZGeoEl::TriangleArea(TPZVec<TPZGeoNode *> &nodes){

  if(nodes.NElements() != 3 && nodes.NElements() != 4){
    PZError << "TPZGeoEl::AreaFromTheFaceT argument error size: nodes\n";
    nodes.Resize(0);
    return 0.;
  }

  REAL cb0 = nodes[2]->Coord(0) - nodes[1]->Coord(0);
  REAL cb1 = nodes[2]->Coord(1) - nodes[1]->Coord(1);
  REAL cb2 = nodes[2]->Coord(2) - nodes[1]->Coord(2);
  REAL ab0 = nodes[0]->Coord(0) - nodes[1]->Coord(0);
  REAL ab1 = nodes[0]->Coord(1) - nodes[1]->Coord(1);
  REAL ab2 = nodes[0]->Coord(2) - nodes[1]->Coord(2);

  //produto vetorial
  REAL coord0 = cb1*ab2-ab1*cb2;
  REAL coord1 = ab0*cb2-cb0*ab2;
  REAL coord2 = cb0*ab1-ab0*cb1;
  //norma da metade do vetor
  return ( 0.5*sqrt(coord0*coord0+coord1*coord1+coord2*coord2) );
}

REAL TPZGeoEl::QuadArea(TPZVec<TPZGeoNode *> &nodes){

  if(nodes.NElements() != 4){
    PZError << "TPZGeoEl::AreaFromTheFaceQ argument error size: nodes\n";
    nodes.Resize(0);
    return 0.;
  }

  REAL areat1 = TriangleArea(nodes);

  nodes[1] = nodes[0];

  nodes[0] = nodes[2];
  nodes[2] = nodes[1];//antigo node0
  nodes[1] = nodes[3];

  REAL areat2 = TriangleArea(nodes);

  return (areat1+areat2);
}

REAL TPZGeoEl::Volume(){

  TPZVec<REAL> param(3,0.);
  REAL detjac;
  TPZFMatrix jacinv(3,3),jacobian(3,3),axes(3,3);
  //supondo jacobiano constante: X linear
  CenterPoint(NSides()-1,param);
  Jacobian(param,jacobian,axes,detjac,jacinv);
  return (RefElVolume()*detjac);//RefElVolume(): volume do elemento de refer�ncia
}

REAL TPZGeoEl::SideArea(int side){

  if(side < 0 || side > NSides()-1)
    PZError << "TPZGeoEl::AreaFromTheFace side error, side = " << side << endl;

  if(SideDimension(side) != 2)
    cout << "TPZGeoEl::SideArea not implemented for side = " << side << endl;

  if(SideDimension(side) == 2){

    int nsn = NSideNodes(side);

    TPZVec<TPZGeoNode *> nodes(nsn);
    int i;

    for(i=0;i</*3*/nsn;i++)
      nodes[i] = &Mesh()->NodeVec()[  SideNodeIndex(side,i) ];
    if (nsn==4)
      return ( QuadArea(nodes) );
    else
      return ( TriangleArea(nodes) );
  }
  return 0.;
}

TPZCompEl *TPZGeoEl::CreateBCCompEl(int side,int bc,TPZCompMesh &cmesh) {
	TPZGeoEl *gel = CreateBCGeoEl(side,bc);
	int index;
	return gel->CreateCompEl(cmesh,index);
}

void TPZGeoEl::RemoveConnectivities(){

  int nsides = NSides(),side;
  for(side=0;side<nsides;side++){
    TPZGeoElSide thisside(this,side);
    thisside.RemoveConnectivity();
  }
}

void TPZGeoEl::InitializeNeighbours(){
  int i,j;
  for (i=0;i<NSides();i++){

    TPZStack <TPZGeoElSide> subel;
    if (HasSubElement()){
      GetSubElements2(i,subel);
      for (j=0;j<subel.NElements();j++){
	TPZGeoEl *el = subel[j].Element();
	el->InitializeNeighbours();
      }
    }
    TPZGeoElSide neighside = Neighbour(i);
    if (!neighside.Element() || neighside.Side() == -1){
      TPZGeoElSide thisside (this,i);
      SetNeighbour(i,thisside);
    }
  }
}

void TPZGeoEl::MidSideNodeIndices(int side,TPZVec<int> &indices) {
  indices.Resize(1);
  MidSideNodeIndex(side,indices[0]);
  if(indices[0] == -1) indices.Resize(0);
}

/** Defines the refinement pattern. It's used only in TPZGeoElRefPattern objects. */
void TPZGeoEl::SetRefPattern(TPZAutoPointer<TPZRefPattern> ){
  PZError << "TPZGeoEl::SetRefPattern ERROR : Should not be called in TPZGeoEl" << endl;
}

void TPZGeoEl::Read(TPZStream &buf, void *context) {
  TPZSaveable::Read(buf, context);
  this->fMesh = static_cast<TPZGeoMesh *>(context);
  buf.Read(&fId,1);
  buf.Read(&fIndex,1);
  buf.Read(&fFatherIndex,1);
  buf.Read(&fMatId,1);
}

void TPZGeoEl::Write(TPZStream &buf, int withclassid) {
  TPZSaveable::Write(buf,withclassid);
  buf.Write(&fId,1);
  buf.Write(&fIndex,1);
  buf.Write(&fFatherIndex,1);
  buf.Write(&fMatId,1);
}

TPZGeoEl::TPZGeoEl(TPZGeoMesh & DestMesh, const TPZGeoEl &cp):TPZSaveable(cp){
  this->fMesh = &DestMesh;
  this->fId = cp.fId;
  this->fMatId = cp.fMatId;
  this->fReference = cp.fReference;
  this->fFatherIndex = cp.fFatherIndex;
  this->fIndex = cp.fIndex;
  this->fMesh->ElementVec()[this->fIndex] = this;
  this->fNumInterfaces = 0;
}

TPZGeoEl::TPZGeoEl(TPZGeoMesh & DestMesh, const TPZGeoEl &cp, std::map<int,int> &org2clnMap):TPZSaveable(cp){
  this->fMesh = &DestMesh;
  this->fId = cp.fId;
  this->fMatId = cp.fMatId;
  this->fReference = 0;
  if ( cp.fFatherIndex == -1) this->fFatherIndex = -1;
  else if (org2clnMap.find(cp.fFatherIndex) == org2clnMap.end())
  {
    std::stringstream sout;
    sout << "ERROR in - " << __PRETTY_FUNCTION__
        << " original father element index: " << cp.fIndex << " is not mapped!";
    LOGPZ_ERROR (logger,sout.str().c_str());
    exit(-1);
  }
  else this->fFatherIndex = org2clnMap[cp.fFatherIndex];

  if (org2clnMap.find(cp.fIndex) == org2clnMap.end())
  {
    std::stringstream sout;
    sout << "ERROR in - " << __PRETTY_FUNCTION__
         << " original element index: " << cp.fIndex << " is not mapped!";
    LOGPZ_ERROR (logger,sout.str().c_str());
    exit(-1);
  }

  this->fIndex = org2clnMap[cp.fIndex];
  this->fMesh->ElementVec()[this->fIndex] = this;
  this->fNumInterfaces = 0;
}

  /// return the refinement pattern associated with the element
TPZAutoPointer<TPZRefPattern> TPZGeoEl::GetRefPattern()
{
  TPZAutoPointer<TPZRefPattern> result;
  return result;
}

void TPZGeoEl::VerifyNodeCoordinates(REAL tol){
  const int nnodes = this->NCornerNodes();
  TPZManVector<REAL,3> qsi(this->Dimension());
  TPZManVector<REAL,3> MappedX(3), NodeX(3);
  for(int inode = 0; inode < nnodes; inode++){
    for(int dim = 0; dim < 3; dim++){
      NodeX[dim] = this->NodePtr(inode)->Coord(dim);
    }///dim
    this->CenterPoint(inode,qsi);
    this->X(qsi,MappedX);
    double error = 0.;
    for(int dim = 0; dim < 3; dim++){
      error += (NodeX[dim]-MappedX[dim])*(NodeX[dim]-MappedX[dim]);
    }///dim
    error = sqrt(error);
    if(error > tol){
      std::stringstream mess;
      mess << "FATAL ERROR AT " << __PRETTY_FUNCTION__ << " - Node coordinate differs from mapped node.\n";
      this->Print(mess);
      PZError << mess.str() << "\n";
#ifdef LOG4CXX
      LOGPZ_ERROR(logger,mess.str().c_str());
#endif
      DebugStop();
    }
  }///for i
}///method


#include "tpzgeoelrefpattern.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzgeopoint.h"
#include "pzrefpoint.h"
#include "pzshapepoint.h"

using namespace pzgeom;
using namespace pzrefine;
using namespace pzshape;

TPZGeoEl *TPZGeoEl::CreateGeoElement(MElementType type,
                                       TPZVec<int>& nodeindexes,
                                       int matid,
                                       int& index)
{
  TPZGeoMesh &mesh = *Mesh();
  if(!&mesh) return 0;
  
  switch( type ){
    case 0://point
    {
      TPZGeoElRefPattern<TPZGeoPoint> * gel =
          new TPZGeoElRefPattern<TPZGeoPoint> (nodeindexes, matid, mesh, index);
      return gel;
    }
    case 1://line
    {
      TPZGeoElRefPattern < TPZGeoLinear > *gel =
          new TPZGeoElRefPattern < TPZGeoLinear >
          (nodeindexes, matid, mesh, index);
      return gel;
    }
    case 2://triangle
    {
      TPZGeoElRefPattern < TPZGeoTriangle > *gel =
          new TPZGeoElRefPattern < TPZGeoTriangle >
          (nodeindexes, matid, mesh, index);
      return gel;
    }
    case 3://quadrilatera
    {
      TPZGeoElRefPattern < TPZGeoQuad > * gel =
          new TPZGeoElRefPattern < TPZGeoQuad >
          (nodeindexes, matid, mesh, index);
      return gel;
    }
    case 4://tetraedra
    {
      TPZGeoElRefPattern < TPZGeoTetrahedra > *gel =
          new TPZGeoElRefPattern < TPZGeoTetrahedra >
          (nodeindexes, matid, mesh, index);
      return gel;
    }
    case 5://pyramid
    {
      TPZGeoElRefPattern < TPZGeoPyramid > *gel =
          new TPZGeoElRefPattern < TPZGeoPyramid >
          (nodeindexes, matid, mesh, index);
      return gel;
    }
    case 6://prism
    {
      TPZGeoElRefPattern < TPZGeoPrism > *gel =
          new TPZGeoElRefPattern < TPZGeoPrism >
          (nodeindexes, matid, mesh, index);
      return gel;
    }
    case 7://cube
    {
      TPZGeoElRefPattern < TPZGeoCube > *gel =
          new TPZGeoElRefPattern < TPZGeoCube >
          (nodeindexes, matid, mesh, index);
      return gel;
    }
    default:
    {
      PZError << "TPZGeoMesh::CreateGeoElement type element not exists:"
          << " type = " << type << std::endl;
      return NULL;
    }
  }
}

int ConjugateSide(TPZGeoEl *gel, int side, TPZStack<int> &allsides, int dimension); 

void NormalVector(TPZGeoElSide &LC, TPZGeoElSide &LS, TPZVec<REAL> &normal)
{
}

void Normalize(TPZVec<REAL> &normlow, TPZVec<REAL> &normal)
{
}

void TPZGeoEl::ComputeNormals(TPZMatrix &normals)
{
	int numbernormals = 0;
	int dimension = Dimension();
	int is;
	int nsides = NSides();
	for(is=0; is<nsides; is++)
	{
		if(SideDimension(is) == dimension-1)
		{
			TPZStack<int> lowdim;
			LowerDimensionSides(is,lowdim);
			numbernormals += lowdim.NElements();
		}
	}
	normals.Redim(3, numbernormals);
	int counter = 0;
	for(is=0; is<nsides; is++)
	{
		if(SideDimension(is) == dimension-1)
		{
			TPZStack<int> lowdim;
			LowerDimensionSides(is,lowdim);
			lowdim.Push(is);
			int nlowdim = lowdim.NElements();
			int lowis;
			for(lowis=0; lowis < nlowdim; lowis++)
			{
				int lowsidedimension = SideDimension(lowdim[lowis]);
				int conj_side = ConjugateSide(this,lowdim[lowis],lowdim,lowsidedimension+1);
				TPZGeoElSide LC(this,conj_side);
				TPZGeoElSide LS(this,lowdim[lowis]);
				TPZManVector<REAL> normal(3,0.);
				NormalVector(LC,LS,normal);
				int d;
				for(d=0; d<3; d++) normals(d,counter) = normal[d];
				counter++;
			}
			TPZManVector<REAL> normal(3,0.);
			int d;
			for(d=0; d<3; d++) normal[d] = normals(d,counter-1);
			counter -= nlowdim;
			for(lowis = counter; lowis < counter+nlowdim; lowis++)
			{
				TPZManVector<REAL> normlow(3,0.);
				for(d=0; d<3; d++) normlow[d] = normals(d,lowis);
				Normalize(normlow,normal);
				for(d=0; d<3; d++) normals(d,lowis) = normlow[d];
			}
		}		
	}
}

int ConjugateSide(TPZGeoEl *gel, int side, TPZStack<int> &allsides, int dimension)
{
	std::set<int> allside;
	allside.insert(&allsides[0],&allsides[allsides.NElements()-1]);
	TPZStack<TPZGeoElSide> highsides;
	int targetdimension = gel->SideDimension(side)+1;
	gel->AllHigherDimensionSides(side,targetdimension,highsides);
	int nhigh = highsides.NElements();
	int is;
	for(is=0; is<nhigh; is++)
	{
		int highside = highsides[is].Side();
		if(allside.find(highside) == allside.end())
		{
			return highside;
		}
	}
	return -1;
}
void TPZGeoEl::SetNeighbourForBlending(int side){
  if( !this->IsGeoBlendEl() ) return;

  TPZGeoElSide ElemSide(this,side);
  TPZGeoElSide NextSide(this,side);

  ///trying to get a self non-linear neighbour
  while(NextSide.Neighbour().Element() != ElemSide.Element())
  {
    if(NextSide.Neighbour().Exists() && !NextSide.Neighbour().Element()->IsLinearMapping() && !NextSide.Neighbour().Element()->IsGeoBlendEl())
    {
      if(NextSide.Neighbour().IsRelative(ElemSide) == false){
        if(NextSide.Neighbour().Element()->IsGeoElMapped() == false){
          TPZGeoElSide NeighSide = NextSide.Neighbour();
          TPZTransform NeighTransf(NeighSide.Dimension(),NeighSide.Dimension());
          ElemSide.SideTransform3(NeighSide,NeighTransf);
          this->SetNeighbourInfo(side,NeighSide,NeighTransf);
          return;;
        }/// !TPZGeoElMapped
      }///if IsRelative == false
    }
    NextSide = NextSide.Neighbour();
  }

  ///now TPZGeoElMapped are accepted
  while(NextSide.Neighbour().Element() != ElemSide.Element())
  {
    if(NextSide.Neighbour().Exists() && !NextSide.Neighbour().Element()->IsLinearMapping() && !NextSide.Neighbour().Element()->IsGeoBlendEl())
    {
      if(NextSide.Neighbour().IsRelative(ElemSide) == false){
        TPZGeoElSide NeighSide = NextSide.Neighbour();
        TPZTransform NeighTransf(NeighSide.Dimension(),NeighSide.Dimension());
        ElemSide.SideTransform3(NeighSide,NeighTransf);
        this->SetNeighbourInfo(side,NeighSide,NeighTransf);
        return;;
      }///if IsRelative == false
    }
    NextSide = NextSide.Neighbour();
  }

}///method

void TPZGeoEl::BuildBlendConnectivity(){
  if( !this->IsGeoBlendEl() ) return;
  const int nsides = this->NSides();
  for(int byside = this->NNodes(); byside < nsides; byside++)
  {
    this->SetNeighbourForBlending(byside);
  }///for byside
}///method

