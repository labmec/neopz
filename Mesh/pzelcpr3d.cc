// -*- c++ -*-
//$Id: pzelcpr3d.cc,v 1.9 2003-11-06 19:15:19 cesar Exp $
#include "pzelcpr3d.h"
#include "pzelct3d.h"
#include "pzgeoel.h"
//#include "pzelgpr3d.h"
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
#include "pzshapeprism.h"
#include <math.h>
#include <stdio.h>

//TShortMatrix TPZCompElPr3d::gEqNumbers;
//linha i é filho i, {a,b,c} = {lado do filho atual,irmão vizinho,lado do vizinho}
int TPZCompElPr3d::InNeigh[8][19][3] =  { {{1,1,0},{2,2,0},{3,4,0},{4,4,1},{5,4,2},{7,3,14},{10,1,9},{11,2,9},{12,4,6},{13,4,7},{14,4,8},{17,3,18},{19,4,15},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
					  {{0,3,5},{2,2,1},{3,7,5},{4,5,1},{5,2,4},{8,3,13},{9,3,11},{11,2,10},{12,5,6},{13,5,7},{14,3,7},{18,3,17},{19,5,15},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
                                          {{0,3,3},{1,3,4},{3,7,3},{4,3,1},{5,6,2},{6,3,12},{9,3,9},{10,3,10},{12,3,6},{13,6,7},{14,6,8},{16,3,16},{19,6,15},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
                                          {{0,0,5},{1,5,2},{2,0,4},{3,0,2},{4,1,2},{5,0,1},{6,6,6},{7,5,8},{8,0,13},{9,0,11},{10,1,11},{11,0,10},{12,2,6},{13,1,8},{14,0,7},{15,7,19},{16,2,16},{17,1,18},{18,0,17}},
                                          {{0,0,3},{1,5,0},{2,6,0},{4,5,3},{5,6,3},{6,0,12},{7,7,14},{8,0,14},{10,5,9},{11,6,9},{13,7,8},{15,0,19},{17,7,18},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
                                          {{0,1,3},{1,1,4},{2,6,1},{3,7,2},{5,6,4},{6,1,12},{7,1,13},{8,7,13},{9,7,11},{11,6,10},{14,7,7},{15,1,19},{18,7,17},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
                                          {{0,2,3},{1,7,4},{2,2,5},{3,7,0},{4,7,1},{6,7,12},{7,2,13},{8,2,14},{9,7,9},{10,7,10},{12,7,6},{15,2,19},{16,7,16},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
                                          {{0,4,5},{1,5,5},{2,4,4},{3,3,0},{4,1,5},{5,3,2},{6,6,12},{7,5,14},{8,4,13},{9,4,11},{10,5,11},{11,4,10},{12,2,12},{13,1,14},{14,3,8},{16,6,16},{17,5,18},{18,4,17},{19,3,15}} };

//int TPZCompElPr3d::CornerSons[8][6] = { {0,9,11,6,12,14},{9,1,10,12,7,13},{11,10,2,14,13,8},{14,13,12,11,10,9},{6,12,14,3,15,17},{12,7,13,15,4,16},{14,13,8,17,16,5},{17,16,15,14,13,12}};
//filhos : 0                1                2                  3                 4                  5                 6                   7
int TPZCompElPr3d::CornerSons[8][6] = { {0,6,8,9,15,17},{6,1,7,15,10,16},{8,7,2,17,16,11},{17,16,15,8,7,6},{9,15,17,3,12,14},{15,10,16,12,4,13},{17,16,11,14,13,5},{14,13,12,17,16,15}};
int TPZCompElPr3d::FaceSons[5][4] = { {0,1,2,3},{0,1,5,4},{1,2,6,5},{0,2,6,4},{4,5,6,7} };//de cada sub da face
//F15       F16       F17       F18       F19
int TPZCompElPr3d::FaceInRib[5][4] = { {7,8,6,-1},{10,12,9,6},{11,13,10,7},{11,14,9,8},{13,14,12,-1} };
//F15       F16           F17         F18           F19
//int TPZCompElPr3d::MiddleFace[4] = {11,10,13,12 };//face do elemento associada a face do subelemento pirâmide contida nela
int TPZCompElPr3d::RibSons[9][2] = { {0,1},{1,2},{2,0},{0,4},{1,5},{2,6},{4,5},{5,6},{6,4} };

int TPZCompElPr3d::FaceNodes[5][4]  = { {0,1,2,-1},{0,1,4,3},{1,2,5,4},{0,2,5,3},{3,4,5,-1} };
//F15        F16       F17       F18        F19
int TPZCompElPr3d::SideNodes[9][2]  = { {0,1},{1,2},{2,0},{0,3},{1,4},{2,5},{3,4},{4,5},{5,3} };
//arestas   6     7      8    9     10    11    12    13    14
//nós medios dos lados como pertencentes a algum filho : (filho,nó)
int TPZCompElPr3d::MidSideNodes[12][2]  = { {0,1},{1,2},{0,2},{0,3},{1,4},{2,5},{4,4},{5,5},{4,5},{0,4},{1,5},{2,3} };
//6    7     8     9    10     11   12     13    14
int TPZCompElPr3d::FaceSides[5][4] = { {6,7,8,-1},{6,10,12,9},{7,11,13,10},{8,11,14,9},{12,13,14,-1} };
//F15         F16         F17          F18          F19
REAL TPZCompElPr3d::MidCoord[12][3] = { {.5,0.,-1.},{.5,.5,-1.},{0.,.5,-1.},{0.,0.,0.},{1.,0.,0.},{0.,1.,0.},{.5,0.,1.},{.5,.5,1.},{0.,.5,1.},{.5,0.,0.},{.5,.5,0.},{0.,.5,0.} };
//6            7          8           9          10         11         12        13          14          15         16       17
REAL TPZCompElPr3d::MasterCoord[6][3] = { {0.,0.,-1.},{1.,0.,-1.},{0.,1.,-1.},{0.,0.,1.},{1.,0.,1.},{0.,1.,1.} };
//c0            c1          c2          c3          c4         c5
int TPZCompElPr3d::ShapeFaceId[5][4] = { {0,1,2,-1},{0,1,4,3},{1,2,5,4},{0,2,5,3},{3,4,5,-1} };
//F15        F16       F17       F18       F19
//int TPZCompElPr3d::FaceConnectLocId[5][9] = { {0,1,2,6,7,8,15,-1,-1},{0,1,4,3,6,10,12,9,16},
//					      {1,2,5,4,7,11,13,10,17},{0,2,5,3,8,11,14,9,18},{3,4,5,12,13,14,19,-1,-1} };

TPZCompElPr3d::TPZCompElPr3d(TPZCompMesh &mesh,TPZGeoEl *ref,int &index) :
  TPZInterpolatedElement(mesh,ref,index), fIntRule(2) {
  int i;                                     //dois pontos por eixo
  for(i=0; i<21; i++) fConnectIndexes[i]=-1;
  for(i=0;i<15;i++) {
    fSideOrder[i] = gOrder;
    fPreferredSideOrder[i] = gOrder;
  }
  ref->SetReference(this);
  for(i=0;i<6;i++) {//6 cantos
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
  }
  for(;i<21;i++) {//9 arestas e 5 faces 1 centro
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
    IdentifySideOrder(i);
  }
  TPZVec<int> order(3,2*fSideOrder[20]+2);//integra variavel ao quadrado em cada direção
  fIntRule.SetOrder(order);
}

TPZCompElPr3d::TPZCompElPr3d(TPZCompMesh &mesh,TPZGeoEl *ref,int &index,int /*noconnects*/) :
  TPZInterpolatedElement(mesh,ref,index), fIntRule(2) {
  int i;
  for(i=0; i<21; i++) fConnectIndexes[i]=-1;
  for(i=0;i<15;i++) {
    fSideOrder[i] = gOrder;
    fPreferredSideOrder[i] = gOrder;
  }
}

TPZCompElPr3d::TPZCompElPr3d(TPZCompMesh &mesh, const TPZCompElPr3d &copy) :
  TPZInterpolatedElement(mesh,copy), fIntRule(copy.fIntRule) {
  int i;
  for(i=0; i<21; i++) fConnectIndexes[i]= copy.fConnectIndexes[i];
  for(i=0;i<15;i++) {
    fSideOrder[i] = copy.fSideOrder[i];
    fPreferredSideOrder[i] = copy.fPreferredSideOrder[i];
  }
}

void TPZCompElPr3d::Shape(TPZVec<REAL> &x, TPZFMatrix &phi, TPZFMatrix &dphi) {
  TPZVec<int> id(6);
  int i;
  for(i=0;i<6;i++) id[i] = fReference->NodePtr(i)->Id();// fReference->NodeIndex(i);

  TPZManVector<int> ord(15);
  for(i=0; i<15; i++) ord[i] = fSideOrder[i];
  TPZShapePrism::Shape(x,id,ord,phi,dphi);
}

int TPZCompElPr3d::NConnectShapeF(int side) {
  if(side<6) return 1;//0 a 4
  int s = side-6;//s = 0 a 14 ou side = 6 a 20
  if(side<15) return (fSideOrder[s]-1);//6 a 14
  if(side==15 || side==19) {
    return ((fSideOrder[s]-2)*(fSideOrder[s]-1)/2);
  }
  if(side>15 && side<19) {//16,17,18
    return ((fSideOrder[s]-1)*(fSideOrder[s]-1));
  }
  if(side==20) {
    return ((fSideOrder[s]-2)*(fSideOrder[s]-1)*(fSideOrder[s]-1)/2);
  }
  PZError << "TPZCompElPr3d::NConnectShapeF, bad parameter side " << side << endl;
  return 0;
}

int TPZCompElPr3d::NSideConnects(int side) {
	return TPZShapePrism::NSideConnects(side);
}
//até aqui até aqui até aqui até aqui até aqui até aqui até aqui até aqui até aqui
int TPZCompElPr3d::SideConnectLocId(int node, int side) {
	return TPZShapePrism::SideConnectLocId(side,node);
}

int TPZCompElPr3d::SideOrder(int side) {
  //0 <= side <= 20
  if(side<6 || side>20) return 0;//cantos ou side ruim
  return fSideOrder[side-6];//deve tirar os cantos

}

void TPZCompElPr3d::FaceOrder(int face,int &ord1,int &ord2) {

  if (face<0 || face>4) cout << "\nTPZCompElPr3d::FaceOrder called whit face = " << face << endl;
  //ordem minima para cada direção do plano da face
  ord1 = fSideOrder[face+9];
  ord2 = ord1;
  /*
    ord2 = 0;
    if (face == 0 || face == 15) {//15
    ord1 = (fSideOrder[0] <= fSideOrder[1]) ? fSideOrder[0] : fSideOrder[1];
    ord1 = (   ord1       <= fSideOrder[2]) ?     ord1      : fSideOrder[2];
    } else
    if (face == 1 || face == 16) {//16
    ord1 = (fSideOrder[0] <= fSideOrder[6]) ? fSideOrder[0] : fSideOrder[6];
    ord2 = (fSideOrder[4] <= fSideOrder[3]) ? fSideOrder[4] : fSideOrder[3];
    ord1 = (ord1 <= ord2) ? ord1 : ord2;
    } else
    if (face == 2 || face == 17) {//17
    ord1 = (fSideOrder[1] <= fSideOrder[7]) ? fSideOrder[1] : fSideOrder[7];
    ord2 = (fSideOrder[4] <= fSideOrder[5]) ? fSideOrder[4] : fSideOrder[5];
    ord1 = (ord1 <= ord2) ? ord1 : ord2;
    } else
    if (face == 3 || face == 18) {//18
    ord1 = (fSideOrder[8] <= fSideOrder[2]) ? fSideOrder[8] : fSideOrder[2];
    ord2 = (fSideOrder[5] <= fSideOrder[3]) ? fSideOrder[5] : fSideOrder[3];
    ord1 = (ord1 <= ord2) ? ord1 : ord2;
    } else
    if (face == 4 || face == 19) {//19
    ord1 = (fSideOrder[6] <= fSideOrder[7]) ? fSideOrder[6] : fSideOrder[7];
    ord1 = (   ord1       <= fSideOrder[8]) ?     ord1      : fSideOrder[8];
    }
  */
}

void TPZCompElPr3d::Print(ostream &out) {
  TPZInterpolatedElement::Print(out);
  /*   return;
       out << "TPZCompElPr3d element:\n";
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

/*
// compare the level of two elements
short TPZCompElPr3d::CompareLevel(TPZCompElPr3d &element,TPZCompElPr3d &neighbour)	{
  // if the level of elements is equal return 1, else return 0
  if (   ( (element.Reference()) -> Level() )	==
	 ( (neighbour.Reference()) -> Level() )   )	return 1;
  else return 0;
}
*/

void TPZCompElPr3d::NormalVector(int side,TPZVec<REAL> &int_point,
				 TPZVec<REAL> &normal, TPZFMatrix &axes, TPZFMatrix &norm_r3){
//  Reference()->NormalVector(side,int_point,normal,axes,norm_r3);
}

void TPZCompElPr3d::EqNumber(TPZShortMatrix &/*mat*/) {

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

void TPZCompElPr3d::VarRange(int var,REAL &min,REAL &max) {
  PZError << "TPZCompElPr3d::VarRange is not defined.\n";
  if(var>-1) max = min = 0.;
}

void TPZCompElPr3d::Load() {
  PZError << "TPZCompElPr3d::Load is called.\n";
}

void TPZCompElPr3d::SetConnectIndex(int i,int connectindex) {
  if(i>-1 && i<21) fConnectIndexes[i] = connectindex;
  else {
    PZError << "TPZCompElPr3d::SetConnectIndex. Bad parameter i.\n";
    PZError.flush();
  }
}

int TPZCompElPr3d::ConnectIndex(int i) {
  if(i<0 || i>21) {
    PZError << "TCompElT2d::ConnectIndex. Bad parameter i.\n";
    return -1;
  }
  return fConnectIndexes[i];
}

void TPZCompElPr3d::SetInterpolationOrder(TPZVec<int> &ord) {
  if(ord.NElements()!=15)
    PZError << "TPZCompElPr3d::SetInterpolationOrder. ord has bad number of elements.\n";
  for(int i=0;i<14;i++) fPreferredSideOrder[i] = ord[i];
}

void TPZCompElPr3d::GetInterpolationOrder(TPZVec<int> &ord) {
  ord.Resize(15);
  for(int i=0;i<15;i++) ord[i] = fSideOrder[i];
}

TPZIntPoints *TPZCompElPr3d::CreateSideIntegrationRule(int side) {
  if(side<0 || side>20) {
    PZError << "TPZCompElPr3d::CreateSideIntegrationRule. bad side number.\n";
    return 0;
  }
  //SideOrder corrige sides de 5 a 18 para 0 a 13
  if(side<6)   return new TPZInt1Point();//cantos 0 a 5
  if(side<15)  return new TPZInt1d(2*SideOrder(side));//lados 6 a 14
  if(side==15 || side==19) return new TPZIntTriang(2*SideOrder(side));
  if(side<20)  {//faces : 16 a 18
    return new TPZIntQuad(2*SideOrder(15),2*SideOrder(side));
  }
  if(side==20) {//integração do elemento

    int ord = 2*SideOrder(20);
    return new TPZIntPrism3D(ord,ord);
  }
  return new TPZInt1Point();
}

int TPZCompElPr3d::PreferredSideOrder(int side) {
  if(side>-1 && side<6) return 0;//cantos
  if(side<21) {
    int order = fPreferredSideOrder[side-6];//lados,faces e centro (ou interior)
    return AdjustPreferredSideOrder(side,order);
  }
  PZError << "TPZCompElPr3d::PreferredSideOrder called for side = " << side << "\n";
  return 0;
}

void TPZCompElPr3d::SetPreferredSideOrder(int order) {
  for(int side=6; side<21; side++){
    fPreferredSideOrder[side-6] = order;
  }
}

void TPZCompElPr3d::SetSideOrder(int side, int order) {
/*   if(fConnectIndexes[side] ==-1) { */
/*     cout << "TPZCompElPr3d::SetSideOrder called for uninitialized connect\n"; */
/*     return; */
/*   } */
  if(side<0 || side>20 || order<1) {
    PZError << "TPZCompElPr3d::SetSideOrder. Bad paramenter side.\n";
    return;
  }
  if(side>5) fSideOrder[side-6] = order;
  if(fConnectIndexes[side] !=-1) {
    TPZConnect &c = Connect(side);
    c.SetOrder(order);
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

void TPZCompElPr3d::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi) {

  if(side<0 || side>20) PZError << "TPZCompElPr3d::SideShapeFunction. Bad paramenter side.\n";
  else if(side==20) {
    Shape(point,phi,dphi);
    return;
  } else if(side<6) {
    phi(0,0)=1.;
    return;
  }
  int nconnect = NSideConnects(side);
  int nnodes = TPZShapePrism::NSideNodes(side);
  int norder = nconnect - TPZShapePrism::NSideNodes(side);
  TPZVec<int> ids(nnodes);
  TPZVec<int> ord(norder);
  int in;
  for(in=0; in<nnodes; in++) {
    ids[in] = fReference->NodePtr(TPZShapePrism::SideNodeLocId(side,in))->Id();
  }
  int ors;
  for(ors=0; ors<norder; ors++) {
    ord[ors] = fSideOrder[TPZShapePrism::SideConnectLocId(side,nnodes+ors)-TPZShapePrism::NNodes];
  }
  if(side<15) {
/*     TPZVec<int> id(2); */
/*     int s = side-6;                 //nó local */
/*     //       id[0] = fReference->NodeIndex(SideNodes[s][0]); */
/*     //       id[1] = fReference->NodeIndex(SideNodes[s][1]);//nó global */
/*     id[0] = fReference->NodePtr(SideNodes[s][0])->Id(); */
/*     id[1] = fReference->NodePtr(SideNodes[s][1])->Id(); */
/*     TPZShapeLinear::Shape1d(point[0],SideOrder(side),phi,dphi,id); */
    TPZShapeLinear::Shape(point,ids,ord,phi,dphi);
  }
  else if(side == 15 || side == 19) {
    TPZShapeTriang::Shape(point,ids,ord,phi,dphi);
  } else {
    TPZShapeQuad::Shape(point,ids,ord,phi,dphi);
  }
/*   else if(side<20) { */
/*     TPZVec<int> id(4); */
/*     int face = side-15; */
/*     id[0] = fReference->NodePtr(FaceNodes[face][0])->Id();//??? */
/*     id[1] = fReference->NodePtr(FaceNodes[face][1])->Id();//nó global??? */
/*     id[2] = fReference->NodePtr(FaceNodes[face][2])->Id();//??? */
/*     if(face && face<4) {//face quadrilateral */
/*       int faceorder[5];//4 arestas e uma face no quadrilatero */
/*       id[3] = fReference->NodeIndex(FaceNodes[face][3]); */
/*       faceorder[0] = fSideOrder[FaceSides[face][0]];//arestas */
/*       faceorder[1] = fSideOrder[FaceSides[face][1]]; */
/*       faceorder[2] = fSideOrder[FaceSides[face][2]]; */
/*       faceorder[3] = fSideOrder[FaceSides[face][3]]; */
/*       faceorder[4] = fSideOrder[side];//face */
/*       TPZManVector<int> ord(5,faceorder,5); */
/*       TPZShapeQuad::ShapeQuad(point,id,ord,phi,dphi); */
/*     }  */
/*     else { */
/*       id.Resize(3); */
/*       int faceorder[4]; */
/*       faceorder[0] = fSideOrder[FaceSides[face][0]];//arestas */
/*       faceorder[1] = fSideOrder[FaceSides[face][1]]; */
/*       faceorder[2] = fSideOrder[FaceSides[face][2]]; */
/*       faceorder[3] = fSideOrder[side];//face */
/*       TPZManVector<int> ord(4,faceorder,4); */
/*       TPZShapeTriang::ShapeTriang(point,id,ord,phi,dphi); */
/*     } */
/*   } */
}

void TPZCompElPr3d::SetIntegrationRule(int order) {
  TPZIntPrism3D intprisma(order);
  SetIntegrationRule(intprisma);
}

void TPZCompElPr3d::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
  //if(dimension == 3) new TPZGraphEl(this,&grmesh);
}

