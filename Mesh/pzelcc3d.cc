// -*- c++ -*-
#include "pzelcc3d.h"
#include "pzelgc3d.h"
#include "pzmatrix.h"
#include "pzelmat.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzerror.h"
#include "pzconnect.h"
#include "pzgraphel.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapecube.h"
#include "pzgraphelq3dd.h"
#include <math.h>
#include <stdio.h>

//TShortMatrix TPZCompElC3d::gEqNumbers;
//linha i é filho i, {a,b,c} = {lado do filho atual,irmão vizinho,lado do vizinho}
int TPZCompElC3d::InNeigh[8][19][3] = { {{1,1,0},{2,1,3},{3,3,0},{4,4,0},{5,1,4},{6,1,7},{7,3,4},{9,1,11},{10,3,8},{13,1,12},{14,1,15},{15,3,12},{16,4,8},{17,1,19},{18,3,16},{19,4,11},{22,1,24},{23,3,21},{25,4,20}},
					{{0,0,1},{2,2,1},{3,2,0},{4,5,0},{5,5,1},{6,2,5},{7,2,4},{10,2,8},{11,0,9},{12,0,13},{14,2,13},{15,2,12},{16,5,8},{17,5,9},{18,2,16},{19,5,11},{23,2,21},{24,0,22},{25,5,20}},
                                        {{0,3,1},{1,1,2},{3,3,2},{4,3,5},{5,6,1},{6,6,2},{7,6,3},{8,1,10},{11,3,9},{12,3,13},{13,1,14},{15,3,14},{16,6,8},{17,6,9},{18,6,10},{19,6,11},{21,1,23},{24,3,22},{25,6,20}},
                                        {{0,0,3},{1,0,2},{2,2,3},{4,7,0},{5,4,2},{6,2,7},{7,7,3},{8,0,10},{9,2,11},{12,0,15},{13,0,14},{14,2,15},{16,7,8},{17,2,19},{18,7,10},{19,7,11},{21,0,23},{22,2,24},{25,7,20}},
					{{0,0,4},{1,0,5},{2,5,3},{3,0,7},{5,5,4},{6,5,7},{7,7,4},{8,0,16},{9,0,17},{10,0,18},{11,0,19},{13,5,12},{14,5,15},{15,7,12},{17,5,19},{18,7,16},{20,0,25},{22,5,24},{23,7,21}},
                                        {{0,4,1},{1,1,5},{2,1,6},{3,6,0},{4,4,5},{6,6,5},{7,6,4},{8,1,16},{9,1,17},{10,1,18},{11,4,9},{12,4,13},{14,6,13},{15,6,12},{18,6,16},{19,4,17},{20,1,25},{23,6,21},{24,4,22}},
                                        {{0,7,1},{1,5,2},{2,2,6},{3,7,2},{4,7,5},{5,5,6},{7,7,6},{8,5,10},{9,2,17},{10,2,18},{11,7,9},{12,7,13},{13,5,14},{15,7,14},{16,5,18},{19,7,17},{20,2,25},{21,5,23},{24,7,22}},
					{{0,4,3},{1,0,6},{2,3,6},{3,3,7},{4,4,7},{5,4,6},{6,6,7},{8,4,10},{9,3,17},{10,3,18},{11,3,19},{12,4,15},{13,4,14},{14,6,15},{16,4,18},{17,6,19},{20,3,25},{21,4,23},{22,6,24}} };

int TPZCompElC3d::CornerSons[8][8] = { {0,8,20,11,12,21,26,24},{8,1,9,20,21,13,22,26},{20,9,2,10,26,22,14,23},
				       {11,20,10,3,24,26,23,15},{12,21,26,24,4,16,25,19},{21,13,22,26,16,5,17,25},
				       {26,22,14,23,25,17,6,18},{24,26,23,15,19,25,18,7} };
//{2,3,7,6}
int TPZCompElC3d::FaceSons[6][4] = { {0,1,2,3},{0,1,5,4},{1,2,6,5},{3,2,6,7},{0,3,7,4},{4,5,6,7} };
//{2,3,7,6}
int TPZCompElC3d::FaceNodes[6][4]  = { {0,1,2,3},{0,1,5,4},{1,2,6,5},{3,2,6,7},{0,3,7,4},{4,5,6,7} };

int TPZCompElC3d::SideNodes[12][2]  = { {0,1},{1,2},{2,3},{3,0},{0,4},{1,5},
					{2,6},{3,7},{4,5},{5,6},{6,7},{7,4} };
//nós medios dos lados como pertencentes a algum filho : (filho,nó)
int TPZCompElC3d::MidSideNodes[19][2]  = { {0,1},{1,2},{2,3},{3,0},{4,0},{5,1},
					   {6,2},{7,3},{4,5},{5,6},{6,7},{7,4},  //lados
					   {0,2},{0,5},{1,6},{2,7},{0,7},{4,6},  //faces
					   {0,6} };//centro

int TPZCompElC3d::FaceSides[6][4] = { {8,9,10,11},{8,13,16,12},{9,14,17,13},
				      {10,14,18,15},{11,15,19,12},{16,17,18,19} };
//{10,15,18,14}
REAL TPZCompElC3d::MidCoord[19][3] = { {0.,-1.,-1.},{1.,0.,-1.},{0.,1.,-1.},{-1.,0.,-1.},{-1.,-1.,0.},{1.,-1.,0.},
				       {1.,1.,0.},{-1.,1.,0.},{0.,-1.,1.},{1.,0.,1.},{0.,1.,1.},{-1.,0.,1.},
				       {0.,0.,-1.},{0.,-1.,0.},{1.,0.,0.},{0.,1.,0.},{-1.,0.,0.},{0.,0.,1.},{0.,0.,0.} };

REAL TPZCompElC3d::MasterCoord[8][3] = { {-1.,-1.,-1.},{1.,-1.,-1.},{1.,1.,-1.},{-1.,1.,-1.},
					 {-1.,-1.,1.} ,{1.,-1.,1.} ,{1.,1.,1.} ,{-1.,1.,1.}  };

int TPZCompElC3d::ShapeFaceId[6][2] = { {0,2},{0,5},{1,6},{3,6},{0,7},{4,6} };

int TPZCompElC3d::FaceConnectLocId[6][9] = { {0,1,2,3,8,9,10,11,20},{0,1,5,4,8,13,16,12,21},
					     {1,2,6,5,9,14,17,13,22},{3,2,6,7,10,14,18,15,23},//{2,3,7,6,10,15,18,14,23}
					     {0,3,7,4,11,15,19,12,24},{4,5,6,7,16,17,18,19,25} };

TPZCompElC3d::TPZCompElC3d(TPZCompMesh &mesh,TPZGeoEl *ref,int &index) :
  TPZInterpolatedElement(mesh,ref,index), fIntRule(2,2,2) {
  int i;                                     //dois pontos por eixo
  for(i=0; i<27; i++) fConnectIndexes[i]=-1;
  for(i=0;i<19;i++) {
    fSideOrder[i] = gOrder;
    fPreferredSideOrder[i] = gOrder;
  }
  //RemoveSideRestraintsII(EInsert);
  ref->SetReference(this);
  for(i=0;i<8;i++) {//8 cantos
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
  }
  for(;i<27;i++) {//12 lados e 6 faces 1 centro
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
    IdentifySideOrder(i);
  }
  TPZVec<int> order(3,3*fSideOrder[18]);//integra variavel ao cubo em cada direção
  fIntRule.SetOrder(order);
}

TPZCompElC3d::TPZCompElC3d(TPZCompMesh &mesh,TPZGeoEl *ref,int &index,int /*noconnects*/) :
  TPZInterpolatedElement(mesh,ref,index), fIntRule(2,2,2) {
  int i;
  for(i=0; i<27; i++) fConnectIndexes[i]=-1;
  for(i=0;i<19;i++) {
    fSideOrder[i] = gOrder;
    fPreferredSideOrder[i] = gOrder;
  }
  RemoveSideRestraintsII(EInsert);
}

void TPZCompElC3d::Shape(TPZVec<REAL> &x, TPZFMatrix &phi, TPZFMatrix &dphi) {
  TPZVec<int> id(8);
  for(int i=0;i<8;i++) id[i] = fReference->NodePtr(i)->Id();// fReference->NodeIndex(i);
  TPZManVector<int> ord(19);
  int i;
  for(i=0; i<19; i++ ) ord[i] = fSideOrder[i];
  TPZShapeCube::ShapeCube(x,id,ord,phi,dphi);
}

int TPZCompElC3d::NConnectShapeF(int side) {
  //fSideOrder[0 a 18] : lados 8 a 26
  if(side<8) return 1;//0 a 7
  int s = side-8;
  if(side<20) return fSideOrder[s]-1;//8 a 19
  if(side<26) return (fSideOrder[s]-1)*(fSideOrder[s]-1);//20 a 25
  if(side==26) return (fSideOrder[18]-1)*(fSideOrder[18]-1)*(fSideOrder[18]-1);
  PZError << "TPZCompElQ2d::NConnectShapeF, bad parameter side " << side << endl;
  return 0;
}

int TPZCompElC3d::NSideConnects(int side) {
	return TPZShapeCube::NSideConnects(side);
}

// Pronto 23/04/98
int TPZCompElC3d::SideConnectLocId(int node, int side) {
	return TPZShapeCube::SideConnectLocId(side,node);
}

int TPZCompElC3d::SideOrder(int side) {
  //0 <= side <= 26
  if(side<8 || side>26) return 0;//cantos ou side ruim
  return fSideOrder[side-8];//deve tirar os cantos
}

void TPZCompElC3d::FaceOrder(int face,int &ord1,int &ord2) {

  if (face<0 || face>5) cout << "\nTPZCompElC3d::FaceOrder called whit face = " << face << endl;
  //ordem minima para cada direção do plano da face
  ord1 = fSideOrder[12+face];
  ord2 = ord1;
  /*
    if (face == 0) {
    ord1 = (fSideOrder[0] <= fSideOrder[ 2]) ? fSideOrder[0] : fSideOrder[2];
    ord2 = (fSideOrder[1] <= fSideOrder[ 3]) ? fSideOrder[1] : fSideOrder[3];
    } else
    if (face == 1) {
    ord1 = (fSideOrder[0] <= fSideOrder[ 8]) ? fSideOrder[0] : fSideOrder[8];
    ord2 = (fSideOrder[4] <= fSideOrder[ 5]) ? fSideOrder[4] : fSideOrder[5];
    } else
    if (face == 2) {
    ord1 = (fSideOrder[3] <= fSideOrder[11]) ? fSideOrder[3] : fSideOrder[11];
    ord2 = (fSideOrder[4] <= fSideOrder[ 7]) ? fSideOrder[4] : fSideOrder[ 7];
    } else
    if (face == 3) {
    ord1 = (fSideOrder[2] <= fSideOrder[10]) ? fSideOrder[2] : fSideOrder[10];
    ord2 = (fSideOrder[6] <= fSideOrder[ 7]) ? fSideOrder[6] : fSideOrder[ 7];
    } else
    if (face == 4) {
    ord1 = (fSideOrder[1] <= fSideOrder[ 9]) ? fSideOrder[1] : fSideOrder[9];
    ord2 = (fSideOrder[5] <= fSideOrder[ 6]) ? fSideOrder[5] : fSideOrder[6];
    } else
    if (face == 5) {
    ord1 = (fSideOrder[8] <= fSideOrder[10]) ? fSideOrder[8] : fSideOrder[10];
    ord2 = (fSideOrder[9] <= fSideOrder[11]) ? fSideOrder[9] : fSideOrder[11];
    }
  */
}

void TPZCompElC3d::Print(ostream &out) {
  TPZInterpolatedElement::Print(out);
  /*return;
    out << "TPZCompElC3d element:\n";
    if(fMaterial) {
    out << "material number = " << fMaterial->Id() << "\n";
    } else {
    out << "no material defined\n";
    }
    out << "number of dof nodes = " << 27 << "\n"
    << "nodes : ";
    for (int nod=0; nod<27; nod++) {
    out << "\nNode number = " << nod << " number of shape functions = "
    << NConnectShapeF(nod) << "\n";
    Mesh()->ConnectVec()[fConnectIndexes[nod]].Print(*Mesh(),out);
    }
    out << "\nlocal equation matrix = \n";*/

  //EqNumber(gEqNumbers);
  //gEqNumbers.Print("",out);
}

// compare the level of two elements
//short TPZCompElC3d::CompareLevel(TPZCompElC3d &element,TPZCompElC3d &neighbour)	{
  // if the level of elements is equal return 1, else return 0
//  if (   ( (element.Reference()) -> Level() )	==
//	 ( (neighbour.Reference()) -> Level() )   )	return 1;
//  else return 0;
//}

void TPZCompElC3d::NormalVector(int side,TPZVec<REAL> &int_point,
				TPZVec<REAL> &normal, TPZFMatrix &axes, TPZFMatrix &norm_r3){
//  Reference()->NormalVector(side,int_point,normal,axes,norm_r3);
}

void TPZCompElC3d::EqNumber(TPZShortMatrix &/*mat*/) {

  /*	mat.Redim(SideOrder(0),SideOrder(1));//ordem x e ordem y da face 0
	mat *= 0;
	mat -= 1;
	mat(0,0) = 0;
	mat(1,0) = 1;
	mat(1,1) = 2;
	mat(0,1) = 3;
	int ieq = 4;
	for(int is =0; is < 4; is++) {
	int OrSide = SideOrder(is);
	for(int i=2 ; i < OrSide; i++) {
	switch (is) {
	case 0 : mat(i,0) = ieq;
	break;
	case 1 : mat(1,i) = ieq;
	break;
	case 2 : mat(i,1) = ieq;
	break;
	case 3 : mat(0,i) = ieq;
	break;
	default :
	break;
	}
	ieq++;
	}
	}
	for(int i=2;i<fSideOrder[0];i++) {
	for(int j=2;j<fSideOrder[1];j++) {
	mat(i,j) = ieq++;
	}
	} */
}

void TPZCompElC3d::VarRange(int var,REAL &min,REAL &max) {
  PZError << "TPZCompElC3d::VarRange is not defined.\n";
  if(var>-1) max = min = 0.;
}

void TPZCompElC3d::Load() {
  PZError << "TPZCompElC3d::Load is called.\n";
}

void TPZCompElC3d::SetConnectIndex(int i,int connectindex) {
  if(i>-1 && i<27) fConnectIndexes[i] = connectindex;
  else {
    PZError << "TPZCompElC3d::SetConnectIndex. Bad parameter i.\n";
    PZError.flush();
  }
}

int TPZCompElC3d::ConnectIndex(int i) {
  if(i<0 || i>26) {
    PZError << "TCompElT2d::ConnectIndex. Bad parameter i.\n";
    return -1;
  }
  return fConnectIndexes[i];
}

void TPZCompElC3d::SetInterpolationOrder(TPZVec<int> &ord) {
  if(ord.NElements()!=19)
    PZError << "TPZCompElC3d::SetInterpolationOrder. ord has bad number of elements.\n";
  for(int i=0;i<19;i++) fPreferredSideOrder[i] = ord[i];
}

void TPZCompElC3d::GetInterpolationOrder(TPZVec<int> &ord) {
  ord.Resize(19);
  for(int i=0;i<19;i++) ord[i] = fSideOrder[i];
}

TPZIntPoints *TPZCompElC3d::CreateSideIntegrationRule(int side) {
  if(side<0 || side>26) {
    PZError << "TPZCompElC3d::CreateSideIntegrationRule. bad side number.\n";
    return new TPZInt1d(0);
  }
  //SideOrder corrige sides de 8 a 26 para 0 a 18
  if(side<8)   return new TPZInt1Point();//cantos 0 a 7
  if(side<20)  return new TPZInt1d(2*SideOrder(side));//lados 8 a 19
  if(side<26)  {//faces : 20 a 25
    int intord = 2*SideOrder(side);
    return new TPZIntQuad(intord,intord);
  }
  if(side==26) {//integração do elemento
    int intord = 3*SideOrder(26);
    return new TPZIntCube3D(intord,intord,intord);
  }
  return new TPZInt1d(0);
}

int TPZCompElC3d::PreferredSideOrder(int side) {
  if(side>-1 && side<8) return 0;//cantos

  if(side<27) {
    int order = fPreferredSideOrder[side-8];//lados,faces e centro (ou interior)
    return AdjustPreferredSideOrder(side,order);
  }
  PZError << "TPZCompElC3d::PreferredSideOrder called for side = " << side << "\n";
  return 0;
}

void TPZCompElC3d::SetPreferredSideOrder(int side, int order) {
  if(side>-1 && side<8) return;//cantos
  if(side<27) fPreferredSideOrder[side-8] = order;//sides 8 a 26
  else PZError << "TPZCompElC3d::SetPreferredSideOrder called for side = " << side << "\n";
  return;
}

void TPZCompElC3d::SetSideOrder(int side, int order) {
  //  if(fConnectIndexes[side] ==-1) {
  //  cout << "TPZCompElC3d::SetSideOrder called for uninitialized connect\n";
  //  return;
  // }
	if(side<0 || side>26 || order<1) {
      PZError << "TPZCompElC3d::SetSideOrder. Bad paramenter side: side " << side << " order: " << order << ".\n";
	  return;
	}

  if(side>7) {
    fSideOrder[side-8] = order;
    if(fConnectIndexes[side] != -1){
      TPZConnect &c = Connect(side);
      int seqnum = c.SequenceNumber();
      int nvar = 1;
      TPZMaterial *mat = Material();
      if(mat) nvar = mat->NStateVariables();
      Mesh()->Block().Set(seqnum,NConnectShapeF(side)*nvar);
      if (side == NConnects()-1){
	TPZVec<int> ord(3,2*order+2);
	fIntRule.SetOrder(ord);
      }
    }
  }
  cout.flush();
}

void TPZCompElC3d::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi) {

  if(side<0 || side>26) PZError << "TPZCompElC3d::SideShapeFunction. Bad paramenter side.\n";
  else if(side==26) Shape(point,phi,dphi);
  else if(side<8) phi(0,0)=1.;
  else if(side<20) {//8 a 19
    TPZVec<int> id(2);
    int s = side-8;                 //nó local
    //       id[0] = fReference->NodeIndex(SideNodes[s][0]);
    //       id[1] = fReference->NodeIndex(SideNodes[s][1]);//nó global
    id[0] = fReference->NodePtr(SideNodes[s][0])->Id();
    id[1] = fReference->NodePtr(SideNodes[s][1])->Id();
    TPZShapeLinear::Shape1d(point[0],SideOrder(side),phi,dphi,id);
  }
  else if(side<26) {//faces
    TPZVec<int> id(4);
    int face = side-20;
    //       id[0] = fReference->NodeIndex(FaceNodes[face][0]);
    //       id[1] = fReference->NodeIndex(FaceNodes[face][1]);//nó global
    //       id[2] = fReference->NodeIndex(FaceNodes[face][2]);
    //       id[3] = fReference->NodeIndex(FaceNodes[face][3]);
    id[0] = fReference->NodePtr(FaceNodes[face][0])->Id();
    id[1] = fReference->NodePtr(FaceNodes[face][1])->Id();
    id[2] = fReference->NodePtr(FaceNodes[face][2])->Id();
    id[3] = fReference->NodePtr(FaceNodes[face][3])->Id();
    
    TPZManVector<int> ord(5);
    ord[0] = fSideOrder[FaceSides[face][0]-8];//lado
    ord[1] = fSideOrder[FaceSides[face][1]-8];
    ord[2] = fSideOrder[FaceSides[face][2]-8];
    ord[3] = fSideOrder[FaceSides[face][3]-8];
    ord[4] = fSideOrder[side-8];//face
    int ns = (ord[0]+ord[1]+ord[2]+ord[3])+(ord[4]-1)*(ord[4]-1);
    if(ns != phi.Rows()) {
      cout << "TPZCompElC3d::SideShapeFunction achei um bug!!\n";
      cout.flush();
    }
    TPZShapeQuad::ShapeQuad(point,id,ord,phi,dphi);
  }
}

//Cedric 19/03/99
void TPZCompElC3d::SetIntegrationRule(int order) {
  TPZIntCube3D intcube(order,order,order);
  SetIntegrationRule(intcube);
}

void TPZCompElC3d::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
  if(dimension == 3) new TPZGraphElQ3dd(this,&grmesh);
}

/*
void TPZCompElC3d::FaceIdsCube(int face,TPZVec<int> &ids,TPZVec<int> &id,int &id0,int &id1) {

  switch(face) {

  case 0:
    ids[0] = id[0];
    ids[1] = id[1];
    ids[2] = id[2];
    ids[3] = id[3];
    id0 = 0;
    id1 = 2;
    break;

  case 1:
    ids[0] = id[0];
    ids[1] = id[1];
    ids[2] = id[5];
    ids[3] = id[4];
    id0 = 0;
    id1 = 5;
    break;

  case 2:
    ids[0] = id[0];
    ids[1] = id[3];
    ids[2] = id[7];
    ids[3] = id[4];
    id0 = 0;
    id1 = 7;
    break;

  case 3:
    ids[0] = id[2];
    ids[1] = id[3];
    ids[2] = id[7];
    ids[3] = id[6];
    id0 = 3;
    id1 = 6;
    break;

  case 4:
    ids[0] = id[1];
    ids[1] = id[2];
    ids[2] = id[6];
    ids[3] = id[5];
    id0 = 1;
    id1 = 6;
    break;

  case 5:
    ids[0] = id[4];
    ids[1] = id[5];
    ids[2] = id[6];
    ids[3] = id[7];
    id0 = 4;
    id1 = 6;
  }
}
*/
