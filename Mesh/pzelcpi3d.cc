// -*- c++ -*-
//$Id: pzelcpi3d.cc,v 1.9 2003-11-05 16:02:21 tiago Exp $
#include "pzelcpi3d.h"
#include "pzelct3d.h"
//#include "pzelgt3d.h"
//#include "pzelgpi3d.h"
#include "pzgeoel.h"
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
#include "pzshapetriang.h"
#include "pzshapepiram.h"
#include <math.h>
#include <stdio.h>



//TShortMatrix TPZCompElPi3d::gEqNumbers;
//linha i é filho i, {a,b,c} = {lado do filho atual,irmão vizinho,lado do vizinho}
int TPZCompElPi3d::InNeigh[10][18][3] = { {{1,6,1},{2,6,2},{3,9,3},{4,6,0},{6,6,5},{7,9,7},{10,6,4},{11,6,6},{12,9,8},{15,6,10},{16,9,11},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
					  {{0,0,1},{2,7,0},{3,7,3},{4,7,1},{7,7,7},{8,0,6},{9,6,8},{11,7,4},{12,7,8},{16,7,11},{17,6,12},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
                                          {{0,8,1},{1,1,2},{3,8,2},{4,8,3},{5,1,7},{8,8,5},{9,8,8},{10,7,6},{12,8,9},{14,7,13},{17,8,12},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
                                          {{0,0,3},{1,9,0},{2,2,3},{4,9,2},{5,0,7},{6,2,8},{9,9,9},{10,9,6},{11,8,6},{14,9,13},{15,8,10},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
					  {{0,0,4},{1,6,3},{2,5,3},{3,3,4},{5,6,7},{6,7,5},{7,8,7},{8,9,5},{13,5,13},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
                                          {{0,4,1},{1,9,1},{2,8,0},{3,7,2},{4,0,2},{5,4,5},{6,4,8},{7,4,7},{8,4,6},{9,6,9},{10,9,4},{11,8,4},{12,7,9},{13,4,13},{14,6,13},{15,9,10},{16,8,11},{17,7,12}},
                                          {{0,5,1},{1,1,0},{2,1,3},{3,1,4},{4,0,10},{5,1,8},{6,5,10},{7,5,5},{8,1,9},{9,1,12},{10,0,15},{12,1,17},{13,5,14},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
                                          {{0,2,1},{1,5,0},{2,2,4},{3,2,0},{4,1,11},{5,5,8},{6,2,10},{7,2,5},{8,5,9},{9,2,9},{11,1,16},{12,5,17},{13,2,14},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
                                          {{0,4,3},{1,3,1},{2,3,2},{3,4,2},{4,3,10},{5,3,6},{6,3,11},{7,5,7},{8,5,12},{9,2,12},{10,3,15},{11,5,16},{12,2,17},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
                                          {{0,5,4},{1,4,0},{2,5,2},{3,3,0},{4,0,11},{5,5,6},{6,5,11},{7,3,5},{8,0,12},{9,3,9},{10,5,15},{11,0,16},{13,3,14},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}} };

int TPZCompElPi3d::CornerSons[10][5] = { {0,5,13,8,9},{5,1,6,13,10},{13,6,2,7,11},{8,13,7,3,12},{9,10,11,12,4},{10,9,12,11,13},{9,5,13,10,-1},{6,10,11,13,-1},{12,13,7,11,-1},{13,9,12,8,-1} };

int TPZCompElPi3d::FaceSons[5][4] = { {0,1,2,3},{0,1,4,6},{1,2,4,7},{3,2,4,8},{0,3,4,9} };

int TPZCompElPi3d::FaceInRib[5][4] = { {6,7,8,5},{10,9,5,-1},{11,10,6,-1},{11,12,7,-1},{12,9,8,-1} };

int TPZCompElPi3d::MiddleFace[4] = {11,10,13,12 };//face do elemento associada a face do subelemento pirâmide contida nela

int TPZCompElPi3d::FaceNodes[5][4]  = { {0,1,2,3},{0,1,4,-1},{1,2,4,-1},{3,2,4,-1},{0,3,4,-1} };
//{3,2},{0,3}
int TPZCompElPi3d::SideNodes[8][2]  = { {0,1},{1,2},{2,3},{3,0},{0,4},{1,4},{2,4},{3,4} };
//nós medios dos lados como pertencentes a algum filho : (filho,nó)
int TPZCompElPi3d::MidSideNodes[9][2]  = { {0,1},{1,2},{2,3},{3,0},{4,0},{4,1},{4,2},{4,3},{0,2} };

int TPZCompElPi3d::FaceSides[5][4] = { {5,6,7,8},{5,10,9,-1},{6,11,10,-1},{7,11,12,-1},{8,12,9,-1} };

REAL TPZCompElPi3d::MidCoord[9][3] = { {0.,-1.,0.},{1.,0.,0.},{0.,1.,0.},{-1.,0.,0.},{-.5,-.5,.5},{.5,-.5,.5},{.5,.5,.5},{-.5,.5,.5},{0.,0.,0.} };

REAL TPZCompElPi3d::MasterCoord[5][3] = { {-1.,-1.,0.},{1.,-1.,0.},{1.,1.,0.},{-1.,1.,0.},{0.,0.,1.} };

int TPZCompElPi3d::ShapeFaceId[5][4] = { {0,1,2,3},{0,1,4,-1},{1,2,4,-1},{3,2,4,-1},{0,3,4,-1} };

//int TPZCompElPi3d::FaceConnectLocId[5][9] = { {0,1,2,3,5,6,7,8,13},{0,1,4,5,10,9,14,-1,-1},
//					      {1,2,4,6,11,10,15,-1,-1},{3,2,4,7,11,12,16,-1,-1},{0,3,4,8,12,9,17,-1,-1} };

TPZCompElPi3d::TPZCompElPi3d(TPZCompMesh &mesh,TPZGeoEl *ref,int &index) :
  TPZInterpolatedElement(mesh,ref,index), fIntRule(2) {
  int i;                                     //dois pontos por eixo
  for(i=0; i<19; i++) fConnectIndexes[i]=-1;
  for(i=0;i<14;i++) {
    fSideOrder[i] = gOrder;
    fPreferredSideOrder[i] = gOrder;
  }
  ref->SetReference(this);
  for(i=0;i<5;i++) {//5 cantos
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
  }
  for(;i<19;i++) {//8 lados e 5 faces 1 centro
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
    IdentifySideOrder(i);
  }
  TPZVec<int> order(3,2*fSideOrder[18]);//integra variavel ao quadrado em cada direção
  fIntRule.SetOrder(order);
}

TPZCompElPi3d::TPZCompElPi3d(TPZCompMesh &mesh,TPZGeoEl *ref,int &index,int /*noconnects*/) :
  TPZInterpolatedElement(mesh,ref,index), fIntRule(2) {
  int i;
  for(i=0; i<19; i++) fConnectIndexes[i]=-1;
  for(i=0;i<14;i++) {
    fSideOrder[i] = gOrder;
    fPreferredSideOrder[i] = gOrder;
  }
}

TPZCompElPi3d::TPZCompElPi3d(TPZCompMesh &mesh, const TPZCompElPi3d &copy) :
  TPZInterpolatedElement(mesh,copy), fIntRule(copy.fIntRule) {
  int i;
  for(i=0; i<19; i++) fConnectIndexes[i]=copy.fConnectIndexes[i];
  for(i=0;i<14;i++) {
    fSideOrder[i] = copy.fSideOrder[i];
    fPreferredSideOrder[i] = copy.fPreferredSideOrder[i];
  }
}

void TPZCompElPi3d::Shape(TPZVec<REAL> &x, TPZFMatrix &phi, TPZFMatrix &dphi) {
  TPZVec<int> id(5);
  for(int i=0;i<5;i++) id[i] = fReference->NodePtr(i)->Id();//fReference->NodeIndex(i);
  int i;
  TPZManVector<int> ord(14);
  for(i=0; i<14; i++) ord[i] = fSideOrder[i];
  TPZShapePiram::Shape(x,id,ord,phi,dphi);
}

int TPZCompElPi3d::NConnectShapeF(int side) {
  if(side<5) return 1;//0 a 4
  int s = side-5;//s = 0 a 14 ou side = 5 a 18
  if(side<13) return fSideOrder[s]-1;//5 a 12
  if(side==13) return (fSideOrder[s]-1)*(fSideOrder[s]-1);//13
  if(side<18) {//14 a 17
    int sum = 0;
    for(int i=0;i<fSideOrder[s]-1;i++) sum += i;
    return sum;
  }
  if(side==18) {
    int totsum = 0,sum;
    for(int i=1;i<fSideOrder[s]-1;i++) {
      sum = i*(i+1) / 2;
      totsum += sum;
    }
    return totsum;
  }
  PZError << "TPZCompElPi3d::NConnectShapeF, bad parameter side " << side << endl;
  return 0;
}

int TPZCompElPi3d::NSideConnects(int side) {
	return TPZShapePiram::NSideConnects(side);
}

int TPZCompElPi3d::SideConnectLocId(int node, int side) {
	return TPZShapePiram::SideConnectLocId(side,node);
}

void TPZCompElPi3d::FaceOrder(int face,int &ord1,int &ord2) {

  if (face<0 || face>4) cout << "\nTPZCompElPi3d::FaceOrder called whit face = " << face << endl;
  //ordem minima para cada direção do plano da face
  ord1 = fSideOrder[8+face];
  ord2 = ord1;
  /*
    ord2 = 0;
    if (face == 0 || face == 13) {//13
    ord1 = (fSideOrder[0] <= fSideOrder[2]) ? fSideOrder[0] : fSideOrder[2];
    ord2 = (fSideOrder[1] <= fSideOrder[3]) ? fSideOrder[1] : fSideOrder[3];
    ord1 = (ord1 <= ord2) ? ord1 : ord2;
    } else
    if (face == 1 || face == 14) {//14
    ord1 = (fSideOrder[0] <= fSideOrder[5]) ? fSideOrder[0] : fSideOrder[5];
    ord1 = (   ord1       <= fSideOrder[4]) ?     ord1      : fSideOrder[4];
    } else
    if (face == 2 || face == 15) {//15
    ord1 = (fSideOrder[1] <= fSideOrder[6]) ? fSideOrder[1] : fSideOrder[6];
    ord1 = (    ord1       <= fSideOrder[5]) ?     ord1      : fSideOrder[5];
    } else
    if (face == 3 || face == 16) {//16
    ord1 = (fSideOrder[2] <= fSideOrder[7]) ? fSideOrder[2] : fSideOrder[7];
    ord1 = (   ord1       <= fSideOrder[6]) ?     ord1      : fSideOrder[6];
    } else
    if (face == 4 || face == 17) {//17
    ord1 = (fSideOrder[3] <= fSideOrder[7]) ? fSideOrder[3] : fSideOrder[7];
    ord1 = (   ord1       <= fSideOrder[4]) ?     ord1      : fSideOrder[4];
    }
  */
}

/*void TPZCompElPi3d::SetFaceOrder(int face, int order) {
  fFaceOrder[face] = order;//bota ordem de interpolação de face
  }  */

void TPZCompElPi3d::Print(ostream &out) {
  TPZInterpolatedElement::Print(out);
  /*   return;
       out << "TPZCompElPi3d element:\n";
       if(fMaterial) {
       out << "material number = " << fMaterial->Id() << "\n";
       } else {
       out << "no material defined\n";
       }
       out << "number of dof nodes = " << 15 << "\n"
       << "nodes : ";
       for (int nod=0; nod<15; nod++) {
       out << "\nNode number = " << nod << " number of shape functions = "
       << NConnectShapeF(nod) << "\n";
       Mesh()->ConnectVec()[fConnectIndexes[nod]].Print(*Mesh(),out);
       }
       out << "\nlocal equation matrix = \n";  */

  //EqNumber(gEqNumbers);
  //gEqNumbers.Print("",out);
}

// compare the level of two elements
/*
short TPZCompElPi3d::CompareLevel(TPZCompElPi3d &element,TPZCompElPi3d &neighbour)	{
  // if the level of elements is equal return 1, else return 0
  if (   ( (element.Reference()) -> Level() )	==
	 ( (neighbour.Reference()) -> Level() )   )	return 1;
  else return 0;
}
*/
void TPZCompElPi3d::NormalVector(int side,TPZVec<REAL> &int_point,
				 TPZVec<REAL> &normal, TPZFMatrix &axes, TPZFMatrix &norm_r3){
//  Reference()->NormalVector(side,int_point,normal,axes,norm_r3);
}

void TPZCompElPi3d::EqNumber(TPZShortMatrix &/*mat*/) {

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

void TPZCompElPi3d::VarRange(int var,REAL &min,REAL &max) {
  PZError << "TPZCompElPi3d::VarRange is not defined.\n";
  if(var>-1) max = min = 0.;
}

void TPZCompElPi3d::Load() {
  PZError << "TPZCompElPi3d::Load is called.\n";
}

void TPZCompElPi3d::SetConnectIndex(int i,int connectindex) {
  if(i>-1 && i<19) fConnectIndexes[i] = connectindex;
  else {
    PZError << "TPZCompElPi3d::SetConnectIndex. Bad parameter i.\n";
    PZError.flush();
  }
}

int TPZCompElPi3d::ConnectIndex(int i) {
  if(i<0 || i>19) {
    PZError << "TCompElT2d::ConnectIndex. Bad parameter i.\n";
    return -1;
  }
  return fConnectIndexes[i];
}

void TPZCompElPi3d::SetInterpolationOrder(TPZVec<int> &ord) {
  if(ord.NElements()!=14)
    PZError << "TPZCompElPi3d::SetInterpolationOrder. ord has bad number of elements.\n";
  for(int i=0;i<14;i++) fPreferredSideOrder[i] = ord[i];
}

void TPZCompElPi3d::GetInterpolationOrder(TPZVec<int> &ord) {
  ord.Resize(14);
  for(int i=0;i<14;i++) ord[i] = fSideOrder[i];
}

TPZIntPoints *TPZCompElPi3d::CreateSideIntegrationRule(int side) {
  if(side<0 || side>18) {
    PZError << "TPZCompElPi3d::CreateSideIntegrationRule. bad side number.\n";
    return new TPZInt1d(0);
  }
  //SideOrder corrige sides de 5 a 18 para 0 a 13
  if(side<5)   return new TPZInt1Point();//cantos 0 a 3
  if(side<13)  return new TPZInt1d(2*SideOrder(side));//lados 5 a 12
  if(side==13) return new TPZIntQuad(2*SideOrder(13),2*SideOrder(13));
  if(side<18)  {//faces : 14 a 17
    return new TPZIntTriang(2*SideOrder(side));
  }
  if(side==18) {//integração do elemento
    return new TPZIntPyram3D(3*SideOrder(18));
  }
  return new TPZInt1Point();
}

int TPZCompElPi3d::PreferredSideOrder(int side) {
  if(side>-1 && side<5) return 0;//cantos
  if(side<19) {
    int order = fPreferredSideOrder[side-5];//lados,faces e centro (ou interior)
    return AdjustPreferredSideOrder(side,order);
  }
  PZError << "TPZCompElPi3d::PreferredSideOrder called for side = " << side << "\n";
  return 0;
}

void TPZCompElPi3d::SetPreferredSideOrder(int order) {
  for(int side = 5; side<19; side++) fPreferredSideOrder[side-5] = order;//sides 5 a 18
}

int TPZCompElPi3d::SideOrder(int side) {
	if(side < 5) return 0;
	if(side <19) return fSideOrder[side-5];
    PZError << "TPZCompElPi3d::SideOrder. Bad paramenter side.\n";
	return -1;
}

void TPZCompElPi3d::SetSideOrder(int side, int order) {
/*   if(fConnectIndexes[side] ==-1) { */
/*     cout << "TPZCompElPi3d::SetSideOrder called for uninitialized connect\n"; */
/*     return; */
/*   } */
	if(side<0 || side>18 || order<1) {
	    PZError << "TPZCompElPi3d::SetSideOrder. Bad paramenter side.\n";
		return;
	}
  if(side>4) {
    fSideOrder[side-5] = order;
    if(fConnectIndexes[side] !=-1) { 
      TPZConnect &c = Connect(side);
      c.SetOrder(order);
      int seqnum = c.SequenceNumber();
      int nvar = 1;
      TPZMaterial *mat = Material();
      if(mat) nvar = mat->NStateVariables();
      Mesh()->Block().Set(seqnum,NConnectShapeF(side)*nvar);
/*       if (side == NConnects()-1){ */
/* 	TPZVec<int> ord(3,2*order+2); */
/* 	fIntRule.SetOrder(ord); */
/*       } */
    }
  }
  int s;
  int maxorder = 0;
  for (s = 5; s<NConnects()-1; s++){
    if (fSideOrder[s]>maxorder)maxorder=fSideOrder[s];	
  }
}

//void TPZCompElPi3d::CornerShape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {
//  ShapeCornerCube(pt,phi,dphi);
//}

void TPZCompElPi3d::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi) {

  if(side<0 || side>18) PZError << "TPZCompElPi3d::SideShapeFunction. Bad paramenter side.\n";
  else if(side==18) Shape(point,phi,dphi);
  else if(side<5) phi(0,0)=1.;
  else if(side<13) {//5 a 12
    TPZVec<int> id(2);
    TPZManVector<int,1> sideord(1,SideOrder(side));
    int s = side-5;                 //nó local
    //       id[0] = fReference->NodeIndex(SideNodes[s][0]);
    //       id[1] = fReference->NodeIndex(SideNodes[s][1]);//nó global
    id[0] = fReference->NodePtr(SideNodes[s][0])->Id();
    id[1] = fReference->NodePtr(SideNodes[s][1])->Id();
    TPZShapeLinear::Shape(point,id,sideord,phi,dphi);
  }
  else if(side<18) {//faces 13  a 17
    TPZVec<int> id(4);
    int face = side-13;
    id[0] = fReference->NodePtr(FaceNodes[face][0])->Id();
    id[1] = fReference->NodePtr(FaceNodes[face][1])->Id();//nó global
    id[2] = fReference->NodePtr(FaceNodes[face][2])->Id();
    if(face == 0) {//face quadrilateral
      TPZManVector<int> ord(5);
      id[3] = fReference->NodePtr(FaceNodes[face][3])->Id();//fReference->NodeIndex(FaceNodes[face][3]); <-> ERRO DETECTADO PELO TESTE
      ord[0] = fSideOrder[FaceSides[face][0]-5];//arestas
      ord[1] = fSideOrder[FaceSides[face][1]-5];
      ord[2] = fSideOrder[FaceSides[face][2]-5];
      ord[3] = fSideOrder[FaceSides[face][3]-5];
      ord[4] = fSideOrder[side-5];//face
      TPZShapeQuad::Shape(point,id,ord,phi,dphi);
    } else {
      id.Resize(3);
      TPZManVector<int> ord(4);
      ord[0] = fSideOrder[FaceSides[face][0]-5];//arestas
      ord[1] = fSideOrder[FaceSides[face][1]-5];
      ord[2] = fSideOrder[FaceSides[face][2]-5];
      ord[3] = fSideOrder[side-5];//face
      TPZShapeTriang::Shape(point,id,ord,phi,dphi);
    }
  }
}

void TPZCompElPi3d::SetIntegrationRule(int order) {
  TPZIntPyram3D intpiram(order);
  SetIntegrationRule(intpiram);
}

void TPZCompElPi3d::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
  //if(dimension == 3) new TPZGraphEl(this,&grmesh);
}

