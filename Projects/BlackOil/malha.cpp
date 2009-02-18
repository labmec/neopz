//$Id: malha.cpp,v 1.4 2009-02-18 11:54:03 fortiago Exp $

#include "malha.h"
#include "pzblackoil2p3d.h"

#include "TPZGeoElement.h"
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
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"

#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "pzbndcond.h"


TPZCompMesh *Unidimensional(int h, double deltaT){

  const REAL Length = 1000.;
  const REAL Base   = 30.;
  const REAL Height = 30.; 
  const REAL PressaoInjecao = 40.e6;
  const REAL VazaoInj = 2.96083333333333e-06;
  const REAL PressaoProducao = 30.e6;
  const REAL SoInj = 0.;
  const REAL SoProd = 1e12; ///a saturacao na producao nao eh prescrita por ser um outflow

/*
  const int nnode = 12;
  REAL co[nnode][3] = {{0.,0.,0.},    {0.,0., Base},    {0.,Height, Base},    {0., Height, 0.}, 
                   {Length,0.,0.},{Length,0., Base},{Length,Height, Base},{Length, Height, 0.},
                   {0.5*Length,0.,0.},{0.5*Length,0., Base},{0.5*Length,Height, Base},{0.5*Length, Height, 0.}};
  const int nelem = 2;
  int indices[nelem][8] = {{0,1,2,3,8,9,10,11},{8,9,10,11,4,5,6,7}};*/

  const int nnode = 8;
  REAL co[nnode][3] = {{0.,0.,0.},    {0.,0., Base},    {0.,Height, Base},    {0., Height, 0.}, 
                   {Length,0.,0.},{Length,0., Base},{Length,Height, Base},{Length, Height, 0.}
                   };
  const int nelem = 1;
  int indices[nelem][8] = {{0,1,2,3,4,5,6,7}};
  
  TPZGeoEl *elvec[1];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nod;
  for(nod=0; nod < nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,3> coord(3);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    coord[2] = co[nod][2];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el=0; el<nelem; el++) {
    TPZManVector<int,8> nodind(8);
    for(nod=0; nod < 8; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
  }

  int index;
  TPZManVector<int,4> bcincid(4);
  bcincid[0] = 0;
  bcincid[1] = 1;
  bcincid[2] = 2;
  bcincid[3] = 3;
  gmesh->CreateGeoElement(EQuadrilateral,bcincid,-1,index);

  bcincid[0] = 4;
  bcincid[1] = 5;
  bcincid[2] = 6;
  bcincid[3] = 7;  
  gmesh->CreateGeoElement(EQuadrilateral,bcincid,-2,index);

  gmesh->BuildConnectivity();

  TPZVec<TPZGeoEl*> filhos;
  for(int i = 0; i < h; i++){
    int n = gmesh->NElements();
    for(int j = 0; j < n; j++){
      gmesh->ElementVec()[j]->Divide(filhos);
    }
  }

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(3);

  TPZAutoPointer<TPZMaterial> mat = new TPZBlackOil2P3D(1, deltaT);

//   TPZFMatrix val1(2,2,0.),val2(2,1,0.);
//   val2(0,0) = PressaoInjecao;
//   val2(1,0) = SoInj;
//   TPZAutoPointer<TPZMaterial> bcI = mat->CreateBC(mat,-1, 0,val1,val2);

  TPZFMatrix val1(2,2,0.),val2(2,1,0.);
  val2(0,0) = VazaoInj;
  val2(1,0) = SoInj;
  TPZAutoPointer<TPZMaterial> bcI = mat->CreateBC(mat,-1, 3,val1,val2);

  val2(0,0) = PressaoProducao;
  val2(1,0) = SoProd;
  TPZAutoPointer<TPZMaterial> bcP = mat->CreateBC(mat,-2, 1,val1,val2);

  cmesh->InsertMaterialObject(mat); 
  cmesh->InsertMaterialObject(bcI);
  cmesh->InsertMaterialObject(bcP);

  cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetDefaultOrder(0);
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
//   cmesh->CleanUpUnconnectedNodes();
//   cmesh->ExpandSolution();

  ofstream file("malha.txt");
  cmesh->Reference()->Print(file);
  cmesh->Print(file);
  return cmesh;
}

TPZCompMesh *UnidimensionalGravidade(int h, double deltaT){
  const REAL H = 1000.;
  const REAL B   = 20.;
  const REAL L = 30.;

  const int nnode = 12;
  const int nelem = 2;
  REAL co[nnode][3] = { {0,0,0},{B,0,0},{B,L,0},{0,L,0},
                        {0,0,H},{B,0,H},{B,L,H},{0,L,H},
                        {0,0,H/2.},{B,0,H/2.},{B,L,H/2.},{0,L,H/2.}
                      };
  int indices[nelem][8] = {{0,1,2,3,8,9,10,11},{8,9,10,11,4,5,6,7}};
  TPZGeoEl *elvec[1];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nod;
  for(nod=0; nod < nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,3> coord(3);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    coord[2] = co[nod][2];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el=0; el<nelem; el++) {
    TPZManVector<int,8> nodind(8);
    for(nod=0; nod < 8; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
  }

  int index;
  TPZManVector<int,4> bcincid(4);
  bcincid[0] = 0;
  bcincid[1] = 1;
  bcincid[2] = 2;
  bcincid[3] = 3;
  gmesh->CreateGeoElement(EQuadrilateral,bcincid,-1,index);

  bcincid[0] = 4;
  bcincid[1] = 5;
  bcincid[2] = 6;
  bcincid[3] = 7;  
  gmesh->CreateGeoElement(EQuadrilateral,bcincid,-2,index);

  gmesh->BuildConnectivity();

  TPZVec<TPZGeoEl*> filhos;
  for(int i = 0; i < h; i++){
    int n = gmesh->NElements();
    for(int j = 0; j < n; j++){
      gmesh->ElementVec()[j]->Divide(filhos);
    }
  }

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(3);

  TPZAutoPointer<TPZMaterial> mat = new TPZBlackOil2P3D(1, deltaT);

  cmesh->InsertMaterialObject(mat);

  cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetDefaultOrder(0);
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
//   cmesh->CleanUpUnconnectedNodes();
//   cmesh->ExpandSolution();

  ofstream file("malha.txt");
  cmesh->Reference()->Print(file);
  cmesh->Print(file);
  return cmesh;
}

double AreaPoco(TPZGeoMesh * gmesh, int matid){
  const int n = gmesh->NElements();
  double result = 0.;
  for(int i = 0; i < n; i++){
    TPZGeoEl * gel = gmesh->ElementVec()[i];
    if(!gel) continue;
    if(gel->MaterialId() != matid) continue;
    if(gel->HasSubElement()) continue;
    result += gel->SideArea( gel->NSides()-1 );
  }///for i
  return result;
}///method

TPZCompMesh * QuarterFiveSpot(TPZGeoMesh * gmesh, double deltaT){

  const double AreaPocoInjetor = AreaPoco(gmesh,-1);
  const REAL PressaoInjecao = 34813607.;
  const REAL SoInj = 0.;
  const REAL VazaoInjecao = 0.0073605229/AreaPocoInjetor;///para levar de Qsc para Qreal
  bool ImporPressaoInjecao = false;///true pressao de injecao, false vazao de injecao

  const REAL PressaoProducao = 27866541.;
  const REAL SoProd = 1.; ///a saturacao na producao nao eh prescrita por ser um outflow

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(3);

  TPZAutoPointer<TPZMaterial> mat = new TPZBlackOil2P3D(1, deltaT);

  TPZFMatrix val1(2,2,0.),val2(2,1,0.);
  TPZAutoPointer<TPZMaterial> bcI;
  if(ImporPressaoInjecao){
    val2(0,0) = PressaoInjecao;
    val2(1,0) = SoInj;
    bcI = mat->CreateBC(mat,-1, 0,val1,val2);
  }
  else{
    val2(0,0) = VazaoInjecao;
    val2(1,0) = SoInj;
    bcI = mat->CreateBC(mat,-1, 3,val1,val2);
  }

  val2(0,0) = PressaoProducao;
  val2(1,0) = SoProd;
  TPZAutoPointer<TPZMaterial> bcP = mat->CreateBC(mat,-2, 1,val1,val2);

  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bcI);
  cmesh->InsertMaterialObject(bcP);

  cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetDefaultOrder(0);
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
//   cmesh->CleanUpUnconnectedNodes();
//   cmesh->ExpandSolution();

  ofstream file("malha.txt");
  cmesh->Reference()->Print(file);
  cmesh->Print(file);
  return cmesh;
}///method

void DivideMalha(TPZGeoMesh * gmesh){
  gmesh->InitializeRefPatterns();
  TPZVec<TPZGeoEl*> filhos;
  const int n = gmesh->NElements();
  for(int iel = 0; iel < n; iel++){
    TPZGeoEl * gel = gmesh->ElementVec()[iel];
    if(!gel) continue;
    gel->Divide(filhos);
  }///for

}///void

void DivideTornoPocos(TPZGeoMesh * gmesh){
  gmesh->InitializeRefPatterns();
  TPZVec<TPZGeoEl*> filhos;
  const int n = gmesh->NElements();
  set<TPZGeoEl *> dividir;
  for(int iel = 0; iel < n; iel++){
    TPZGeoEl * gel = gmesh->ElementVec()[iel];
    if(!gel) continue;
    if((gel->MaterialId() != -1) && (gel->MaterialId() != -2)) continue;
    TPZGeoElSide gelside(gel, gel->NSides()-1);
    TPZGeoElSide neigh = gelside.Neighbour();
    dividir.insert(gel);
    dividir.insert(neigh.Element());
  }///for

  for(set<TPZGeoEl*>::iterator w = dividir.begin(); w != dividir.end(); w++){
    (*w)->Divide(filhos);
  }
}///void

/// CAJU ///

#include <sstream>

#include "pzvec.h"
#include "pzcmesh.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"

#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "TPZCopySolve.h"
#include "TPZStackEqnStorage.h"

#include "pzbstrmatrix.h"
#include "pzstepsolver.h"

#include "pzbndcond.h"
#include "pzpoisson3d.h"

#include "pzvisualmatrix.h"

#include "TPZGeoElement.h"
#include "pzgeoel.h"

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
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzskylstrmatrix.h"

#include <time.h>
#include <stdio.h>
#include "pzl2projection.h"
#include "tpzgeoelmapped.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <cstdlib>

#include "pzgeotriangle.h"
#include "tpzarc3d.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoelrefpattern.h.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"
#include "TPZGeoCube.h"
#include <pzcompel.h>

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("main"));
#endif


#include <sstream>
using namespace std;
using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;






void ScaleVec(TPZVec<REAL> &NodeIni, TPZVec<REAL> &NodeFin, double Norm, TPZVec<REAL> &OriginVec, TPZVec<REAL> &OutputVec)
{
  TPZVec<REAL> vec(3);
  for(int i = 0; i < 3; i++) vec[i] = NodeFin[i] - NodeIni[i];

  double N = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

  for(int j = 0; j < 3; j++)
  {
    OutputVec[j] = OriginVec[j] + (vec[j]/N)*Norm;
  }
}

TPZGeoMesh * QuarterFiveSpot(int ndiv_xy, int ndiv_z, double lxy, double lz, double wellDiam, double factor)
{
  if(ndiv_z < 1)
  {
  cout << "Quantity of layer(s) must be higher or equal to 1!\n";
  exit(-1);
  }
  
  TPZGeoMesh * Mesh = new TPZGeoMesh;

  int Qnodes = (23 + 18*ndiv_xy) * (ndiv_z + 1);
  TPZVec < TPZVec <REAL> > NodesCoords(Qnodes);
  for(int i = 0; i < Qnodes; i++) NodesCoords[i].Resize(3);

  for(int lay = 0; lay < (ndiv_z + 1); lay++) // Begin NodesGeneration by Layer
  {
      int QTDnodesByLayer = (23 + 18*ndiv_xy);
      int firstNode = QTDnodesByLayer * lay; //First node ID on this layer

      /// Nohs da diagonal
      NodesCoords[0 + firstNode][0] = -lxy/2.;
      NodesCoords[0 + firstNode][1] =  lxy/2.;
      NodesCoords[0 + firstNode][2] = -lz/2. + lay*(lz/ndiv_z);

      NodesCoords[1 + firstNode][0] = -lxy/4.;
      NodesCoords[1 + firstNode][1] =  lxy/4.;
      NodesCoords[1 + firstNode][2] = -lz/2. + lay*(lz/ndiv_z);

      NodesCoords[2 + firstNode][0] =  0.;
      NodesCoords[2 + firstNode][1] =  0.;
      NodesCoords[2 + firstNode][2] = -lz/2. + lay*(lz/ndiv_z);

      NodesCoords[3 + firstNode][0] =  lxy/4.;
      NodesCoords[3 + firstNode][1] = -lxy/4.;
      NodesCoords[3 + firstNode][2] = -lz/2. + lay*(lz/ndiv_z);

      NodesCoords[4 + firstNode][0] =  lxy/2.;
      NodesCoords[4 + firstNode][1] = -lxy/2.;
      NodesCoords[4 + firstNode][2] = -lz/2. + lay*(lz/ndiv_z);
      //

      /// Nohs dos pocos
      TPZVec<REAL> orig(3);
      TPZVec<REAL> temp(3);

      orig[0] = -lxy/2.;
      orig[1] = -lxy/2.;
      orig[2] = -lz/2. + lay*(lz/ndiv_z);

      double R = wellDiam/2.;

      ScaleVec(orig, NodesCoords[0 + firstNode], R, orig, NodesCoords[5 + firstNode]);
      NodesCoords[6 + firstNode][0] = -NodesCoords[5 + firstNode][1];
      NodesCoords[6 + firstNode][1] = -NodesCoords[5 + firstNode][0];
      NodesCoords[6 + firstNode][2] =  NodesCoords[5 + firstNode][2];

      temp[0] = (NodesCoords[0 + firstNode][0] + NodesCoords[1 + firstNode][0])/2.;
      temp[1] = (NodesCoords[0 + firstNode][1] + NodesCoords[1 + firstNode][1])/2.;
      temp[2] = (NodesCoords[0 + firstNode][2]);
      ScaleVec(orig, temp, R, orig, NodesCoords[7 + firstNode]);
      NodesCoords[8 + firstNode][0] = -NodesCoords[7 + firstNode][1];
      NodesCoords[8 + firstNode][1] = -NodesCoords[7 + firstNode][0];
      NodesCoords[8 + firstNode][2] =  NodesCoords[7 + firstNode][2];

      ScaleVec(orig, NodesCoords[1 + firstNode], R, orig, NodesCoords[9 + firstNode]);
      NodesCoords[10 + firstNode][0] = -NodesCoords[9 + firstNode][1];
      NodesCoords[10 + firstNode][1] = -NodesCoords[9 + firstNode][0];
      NodesCoords[10 + firstNode][2] =  NodesCoords[9 + firstNode][2];

      temp[0] = (NodesCoords[1 + firstNode][0] + NodesCoords[2 + firstNode][0])/2.;
      temp[1] = (NodesCoords[1 + firstNode][1] + NodesCoords[2 + firstNode][1])/2.;
      temp[2] = (NodesCoords[1 + firstNode][2]);
      ScaleVec(orig, temp, R, orig, NodesCoords[11 + firstNode]);
      NodesCoords[12 + firstNode][0] = -NodesCoords[11 + firstNode][1];
      NodesCoords[12 + firstNode][1] = -NodesCoords[11 + firstNode][0];
      NodesCoords[12 + firstNode][2] =  NodesCoords[11 + firstNode][2];

      ScaleVec(orig, NodesCoords[2 + firstNode], R, orig, NodesCoords[13 + firstNode]);
      NodesCoords[14 + firstNode][0] = -NodesCoords[13 + firstNode][1];
      NodesCoords[14 + firstNode][1] = -NodesCoords[13 + firstNode][0];
      NodesCoords[14 + firstNode][2] =  NodesCoords[13 + firstNode][2];

      temp[0] = (NodesCoords[2 + firstNode][0] + NodesCoords[3 + firstNode][0])/2.;
      temp[1] = (NodesCoords[2 + firstNode][1] + NodesCoords[3 + firstNode][1])/2.;
      temp[2] = (NodesCoords[2 + firstNode][2]);
      ScaleVec(orig, temp, R, orig, NodesCoords[15 + firstNode]);
      NodesCoords[16 + firstNode][0] = -NodesCoords[15 + firstNode][1];
      NodesCoords[16 + firstNode][1] = -NodesCoords[15 + firstNode][0];
      NodesCoords[16 + firstNode][2] =  NodesCoords[15 + firstNode][2];

      ScaleVec(orig, NodesCoords[3 + firstNode], R, orig, NodesCoords[17 + firstNode]);
      NodesCoords[18 + firstNode][0] = -NodesCoords[17 + firstNode][1];
      NodesCoords[18 + firstNode][1] = -NodesCoords[17 + firstNode][0];
      NodesCoords[18 + firstNode][2] =  NodesCoords[17 + firstNode][2];

      temp[0] = (NodesCoords[3 + firstNode][0] + NodesCoords[4 + firstNode][0])/2.;
      temp[1] = (NodesCoords[3 + firstNode][1] + NodesCoords[4 + firstNode][1])/2.;
      temp[2] = (NodesCoords[3 + firstNode][2]);
      ScaleVec(orig, temp, R, orig, NodesCoords[19 + firstNode]);
      NodesCoords[20 + firstNode][0] = -NodesCoords[19 + firstNode][1];
      NodesCoords[20 + firstNode][1] = -NodesCoords[19 + firstNode][0];
      NodesCoords[20 + firstNode][2] =  NodesCoords[19 + firstNode][2];

      ScaleVec(orig, NodesCoords[4 + firstNode], R, orig, NodesCoords[21 + firstNode]);
      NodesCoords[22 + firstNode][0] = -NodesCoords[21 + firstNode][1];
      NodesCoords[22 + firstNode][1] = -NodesCoords[21 + firstNode][0];
      NodesCoords[22 + firstNode][2] =  NodesCoords[21 + firstNode][2];
      //

      // Nohs dos Arcos
      for(int divxy = 1; divxy <= ndiv_xy; divxy++)
      {
          R = (lxy - wellDiam/2.)/pow(factor,divxy) + wellDiam/2.;

          for(int th = 0; th < 9; th++)
          {
              int NodeId = firstNode + (5 + 18*divxy + 2*th);

              ScaleVec(orig, NodesCoords[NodeId-18], R, orig, NodesCoords[NodeId]);
              NodesCoords[NodeId+1][0] = -NodesCoords[NodeId][1];
              NodesCoords[NodeId+1][1] = -NodesCoords[NodeId][0];
              NodesCoords[NodeId+1][2] =  NodesCoords[NodeId][2];
          }
      }
      //
  } // End NodesGeneration by Layer

  Mesh->NodeVec().Resize(Qnodes);
  TPZVec <TPZGeoNode> Node(Qnodes);
  for(int n = 0; n < Qnodes; n++)
  {
    Node[n].SetNodeId(n);
    Node[n].SetCoord(&NodesCoords[n][0]);
    Mesh->NodeVec()[n] = Node[n];
  }

  TPZVec <int> Topol;

  // Materials IDs
  int reservMat = 1;
  int InjectionBC = -1;
  int ProductionBC = -2;
  int arcMat = -3;

  int id = 0;

  for(int lay = 0; lay < ndiv_z; lay++) // Begin Elements Generation by Layer
  {
      Topol.Resize(8);

      int QTDnodesByLayer = (23 + 18*ndiv_xy);
      int firstNode = QTDnodesByLayer * lay; //First node ID on this layer

      ///Reservoir Cubes
      int p = 0;
      if(ndiv_xy > 0) p = 1;

//el0
      Topol[0] = 0 + firstNode;
      Topol[1] = 5 + firstNode + 18*p;
      Topol[2] = 9 + firstNode + 18*p;
      Topol[3] = 1 + firstNode;
      Topol[4] = 0 + firstNode + QTDnodesByLayer;
      Topol[5] = 5 + firstNode + 18*p + QTDnodesByLayer ;
      Topol[6] = 9 + firstNode + 18*p + QTDnodesByLayer ;
      Topol[7] = 1 + firstNode + QTDnodesByLayer;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
    id++;

//el1
      Topol[0] = 0 + firstNode;
      Topol[1] = 1 + firstNode;
      Topol[2] = 10 + firstNode + 18*p;
      Topol[3] = 6 + firstNode + 18*p;
      Topol[4] = 0 + firstNode + QTDnodesByLayer;
      Topol[5] = 1 + firstNode + QTDnodesByLayer;
      Topol[6] = 10 + firstNode + 18*p + QTDnodesByLayer ;
      Topol[7] = 6 + firstNode + 18*p + QTDnodesByLayer ;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
    id++;

//el2
      Topol[0] = 1 + firstNode;
      Topol[1] = 9 + firstNode + 18*p;
      Topol[2] = 13 + firstNode + 18*p;
      Topol[3] = 2 + firstNode;
      Topol[4] = 1 + firstNode + QTDnodesByLayer;
      Topol[5] = 9 + firstNode + 18*p + QTDnodesByLayer ;
      Topol[6] = 13 + firstNode + 18*p + QTDnodesByLayer ;
      Topol[7] = 2 + firstNode + QTDnodesByLayer;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
    id++;

//el3
      Topol[0] = 1 + firstNode;
      Topol[1] = 2 + firstNode;
      Topol[2] = 14 + firstNode + 18*p;
      Topol[3] = 10 + firstNode + 18*p;
      Topol[4] = 1 + firstNode + QTDnodesByLayer;
      Topol[5] = 2 + firstNode + QTDnodesByLayer;
      Topol[6] = 14 + firstNode + 18*p + QTDnodesByLayer ;
      Topol[7] = 10 + firstNode + 18*p + QTDnodesByLayer ;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
    id++;

//el4
      Topol[0] = 2 + firstNode;
      Topol[1] = 13 + firstNode + 18*p;
      Topol[2] = 17 + firstNode + 18*p;
      Topol[3] = 3 + firstNode;
      Topol[4] = 2 + firstNode + QTDnodesByLayer;
      Topol[5] = 13 + firstNode + 18*p + QTDnodesByLayer ;
      Topol[6] = 17 + firstNode + 18*p + QTDnodesByLayer ;
      Topol[7] = 3 + firstNode + QTDnodesByLayer;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
    id++;

//el5
      Topol[0] = 2 + firstNode;
      Topol[1] = 3 + firstNode;
      Topol[2] = 18 + firstNode + 18*p;
      Topol[3] = 14 + firstNode + 18*p;
      Topol[4] = 2 + firstNode + QTDnodesByLayer;
      Topol[5] = 3 + firstNode + QTDnodesByLayer;
      Topol[6] = 18 + firstNode + 18*p + QTDnodesByLayer ;
      Topol[7] = 14 + firstNode + 18*p + QTDnodesByLayer ;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
    id++;

//el6
      Topol[0] = 3 + firstNode;
      Topol[1] = 17 + firstNode + 18*p;
      Topol[2] = 21 + firstNode + 18*p;
      Topol[3] = 4 + firstNode;
      Topol[4] = 3 + firstNode + QTDnodesByLayer;
      Topol[5] = 17 + firstNode + 18*p + QTDnodesByLayer ;
      Topol[6] = 21 + firstNode + 18*p + QTDnodesByLayer ;
      Topol[7] = 4 + firstNode + QTDnodesByLayer;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
    id++;

//el7
      Topol[0] = 3 + firstNode;
      Topol[1] = 4 + firstNode;
      Topol[2] = 22 + firstNode + 18*p;
      Topol[3] = 18 + firstNode + 18*p;
      Topol[4] = 3 + firstNode + QTDnodesByLayer;
      Topol[5] = 4 + firstNode + QTDnodesByLayer;
      Topol[6] = 22 + firstNode + 18*p + QTDnodesByLayer ;
      Topol[7] = 18 + firstNode + 18*p + QTDnodesByLayer ;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
    id++;

      if(p == 1) //if ndiv_xy == 1, so there is more cubes in this layer to be created!
      {
      int nd = 23 + 18*(ndiv_xy-1);

      //el8
      Topol[0] = 5 + firstNode;
      Topol[1] = 9 + firstNode;
      Topol[2] = (nd+4) + firstNode;
      Topol[3] = nd + firstNode;
      Topol[4] = 5 + firstNode + QTDnodesByLayer;
      Topol[5] = 9 + firstNode + QTDnodesByLayer ;
      Topol[6] = (nd+4) + firstNode + QTDnodesByLayer ;
      Topol[7] = nd + firstNode + QTDnodesByLayer;
      new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
      id++;

      //el9
      Topol[0] = 6 + firstNode;
      Topol[1] = (nd+1) + firstNode;
      Topol[2] = ((nd+1)+4) + firstNode;
      Topol[3] = 10 + firstNode;
      Topol[4] = 6 + firstNode + QTDnodesByLayer;
      Topol[5] = (nd+1) + firstNode + QTDnodesByLayer;
      Topol[6] = ((nd+1)+4) + firstNode + QTDnodesByLayer ;
      Topol[7] = 10 + firstNode + QTDnodesByLayer ;
      new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
      id++;

      //el10
      Topol[0] = 9 + firstNode;
      Topol[1] = 13 + firstNode;
      Topol[2] = (nd+8) + firstNode;
      Topol[3] = (nd+4) + firstNode;
      Topol[4] = 9 + firstNode + QTDnodesByLayer;
      Topol[5] = 13 + firstNode + QTDnodesByLayer ;
      Topol[6] = (nd+8) + firstNode + QTDnodesByLayer ;
      Topol[7] = (nd+4) + firstNode + QTDnodesByLayer;
      new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
      id++;

      //el11
      Topol[0] = 10 + firstNode;
      Topol[1] = ((nd+1)+4) + firstNode;
      Topol[2] = ((nd+1)+8) + firstNode;
      Topol[3] = 14 + firstNode;
      Topol[4] = 10 + firstNode + QTDnodesByLayer;
      Topol[5] = ((nd+1)+4) + firstNode + QTDnodesByLayer;
      Topol[6] = ((nd+1)+8) + firstNode + QTDnodesByLayer ;
      Topol[7] = 14 + firstNode + QTDnodesByLayer ;
      new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
      id++;

      //el12
      Topol[0] = 13 + firstNode;
      Topol[1] = 17 + firstNode;
      Topol[2] = (nd+12) + firstNode;
      Topol[3] = (nd+8) + firstNode;
      Topol[4] = 13 + firstNode + QTDnodesByLayer;
      Topol[5] = 17 + firstNode + QTDnodesByLayer ;
      Topol[6] = (nd+12) + firstNode + QTDnodesByLayer ;
      Topol[7] = (nd+8) + firstNode + QTDnodesByLayer;
      new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
      id++;

      //el13
      Topol[0] = 14 + firstNode;
      Topol[1] = ((nd+1)+8) + firstNode;
      Topol[2] = ((nd+1)+12) + firstNode;
      Topol[3] = 18 + firstNode;
      Topol[4] = 14 + firstNode + QTDnodesByLayer;
      Topol[5] = ((nd+1)+8) + firstNode + QTDnodesByLayer;
      Topol[6] = ((nd+1)+12) + firstNode + QTDnodesByLayer ;
      Topol[7] = 18 + firstNode + QTDnodesByLayer ;
      new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
      id++;

      //el14
      Topol[0] = 17 + firstNode;
      Topol[1] = 21 + firstNode;
      Topol[2] = (nd+16) + firstNode;
      Topol[3] = (nd+12) + firstNode;
      Topol[4] = 17 + firstNode + QTDnodesByLayer;
      Topol[5] = 21 + firstNode + QTDnodesByLayer ;
      Topol[6] = (nd+16) + firstNode + QTDnodesByLayer ;
      Topol[7] = (nd+12) + firstNode + QTDnodesByLayer;
      new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
      id++;

      //el15
      Topol[0] = 18 + firstNode;
      Topol[1] = ((nd+1)+12) + firstNode;
      Topol[2] = ((nd+1)+16) + firstNode;
      Topol[3] = 22 + firstNode;
      Topol[4] = 18 + firstNode + QTDnodesByLayer;
      Topol[5] = ((nd+1)+12) + firstNode + QTDnodesByLayer;
      Topol[6] = ((nd+1)+16) + firstNode + QTDnodesByLayer ;
      Topol[7] = 22 + firstNode + QTDnodesByLayer ;
      new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
      id++;
      }
      for(int nd = ndiv_xy-1; nd > 0; nd--) // if ndiv_xy > 1, so there is more cubes in this layer to be created! (id >=16)
      {
      int ini = 23 + 18*nd;
      for(int t = 0; t < 4; t++)
      {
        Topol[0] = ini + 4*t + firstNode;
        Topol[1] = ini + 4*t + 4 + firstNode;
        Topol[2] = ini + 4*t + 4 + firstNode - 18;
        Topol[3] = ini + 4*t + 4 + firstNode - 22;
        Topol[4] = ini + 4*t + firstNode + QTDnodesByLayer;
        Topol[5] = ini + 4*t + 4 + firstNode + QTDnodesByLayer;
        Topol[6] = ini + 4*t + 4 + firstNode - 18 + QTDnodesByLayer;
        Topol[7] = ini + 4*t + 4 + firstNode - 22 + QTDnodesByLayer;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
        id++;

        Topol[0] = (ini+1) + 4*t + firstNode;
        Topol[1] = (ini+1) + 4*t + 4 + firstNode - 22;
        Topol[2] = (ini+1) + 4*t + 4 + firstNode - 18;
        Topol[3] = (ini+1) + 4*t + 4 + firstNode;
        Topol[4] = (ini+1) + 4*t + firstNode + QTDnodesByLayer;
        Topol[5] = (ini+1) + 4*t + 4 + firstNode - 22 + QTDnodesByLayer;
        Topol[6] = (ini+1) + 4*t + 4 + firstNode - 18 + QTDnodesByLayer;
        Topol[7] = (ini+1) + 4*t + 4 + firstNode + QTDnodesByLayer;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
        id++;
      }
      }

      Topol.Resize(4);
      //Quadrilaterals InjectionBC
      for(int t = 0; t < 4; t++)
      {
              Topol[0] = 5 + firstNode + 4*t;
              Topol[1] = 9 + firstNode + 4*t;
              Topol[2] = 9 + firstNode + 4*t + QTDnodesByLayer;
              Topol[3] = 5 + firstNode + 4*t + QTDnodesByLayer;
              new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,InjectionBC,*Mesh);
              id++;
      }

      //Quadrilaterals ProductionBC
      for(int t = 0; t < 4; t++)
      {
              Topol[0] = 10 + firstNode + 4*t;
              Topol[1] =   6 + firstNode + 4*t;
              Topol[2] =   6 + firstNode + 4*t + QTDnodesByLayer;
              Topol[3] = 10 + firstNode + 4*t + QTDnodesByLayer;
              new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ProductionBC,*Mesh);
              id++;
      }
  } // End Elements Generation by Layer

  Mesh->BuildConnectivity();

// #ifdef LOG4CXX
//   {
//     std::stringstream sout;
//     Mesh->Print(sout);
//     LOGPZ_DEBUG(logger,sout.str());
//   }
// #endif

  return Mesh;
}

TPZGeoMesh * QuarterFiveSpotReg(double lxy, double lz)
{ 
  TPZGeoMesh * Mesh = new TPZGeoMesh;

  int Qnodes = 484;
  TPZVec < TPZVec <REAL> > NodesCoords(Qnodes);
  for(int i = 0; i < Qnodes; i++) NodesCoords[i].Resize(3);

  int QTDnodesByLayer = 121;

  for(int lay = 0; lay < 4; lay++) // Begin NodesGeneration by Layer
  {
    int firstNode = QTDnodesByLayer * lay; //First node ID on this layer

    for(int nx = 0; nx < 11; nx++)
    {
        for(int ny = 0; ny < 11; ny++)
        {
            NodesCoords[11*nx+ny + firstNode][0] = -lxy/2. + lxy/10.*nx;
            NodesCoords[11*nx+ny + firstNode][1] = -lxy/2. + lxy/10.*ny;
            NodesCoords[11*nx+ny + firstNode][2] = -lz/2. + lz/3.*lay;
        }
    }
  }

  Mesh->NodeVec().Resize(Qnodes);
  TPZVec <TPZGeoNode> Node(Qnodes);
  for(int n = 0; n < Qnodes; n++)
  {
      Node[n].SetNodeId(n);
      Node[n].SetCoord(&NodesCoords[n][0]);
      Mesh->NodeVec()[n] = Node[n];
  }

  TPZVec <int> Topol;

  // Materials IDs
  int reservMat = 1;
  int InjectionBC = -1;
  int ProductionBC = -2;
  int arcMat = -3;

  int id = 0;

  for(int lay = 0; lay < 3; lay++) // Begin NodesGeneration by Layer
  {
    int firstNode = QTDnodesByLayer * lay; //First node ID on this layer
    Topol.Resize(8);
    for(int ny = 0; ny < 10; ny++)
    {
        for(int nx = 0; nx < 10; nx++)
        {
            if((nx != 0 || ny != 0) && (nx != 9 || ny != 9))
            {
                ///Reservoir Cubes
                Topol[0] = 11*ny+nx + firstNode;
                Topol[1] = 11*ny+nx+1 + firstNode;
                Topol[2] = 11*(ny+1)+nx+1 + firstNode;
                Topol[3] = 11*(ny+1)+nx + firstNode;
                Topol[4] = 11*ny+nx + firstNode + QTDnodesByLayer;
                Topol[5] = 11*ny+nx+1 + firstNode + QTDnodesByLayer ;
                Topol[6] = 11*(ny+1)+nx+1 + firstNode + QTDnodesByLayer ;
                Topol[7] = 11*(ny+1)+nx + firstNode + QTDnodesByLayer;
                new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
                id++;
            }
        }
    }

    ///Quadrilaterals InjectionBC
    Topol.Resize(4);
    Topol[0] = firstNode + 12;
    Topol[1] = firstNode + 1;
    Topol[2] = firstNode + 1  + QTDnodesByLayer;
    Topol[3] = firstNode + 12 + QTDnodesByLayer;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,InjectionBC,*Mesh);
    id++;

    Topol[0] = firstNode + 11;
    Topol[1] = firstNode + 12;
    Topol[2] = firstNode + 12  + QTDnodesByLayer;
    Topol[3] = firstNode + 11 + QTDnodesByLayer;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,InjectionBC,*Mesh);
    id++;

    ///Quadrilaterals ProductionBC
    Topol[0] = firstNode + 109;
    Topol[1] = firstNode + 108;
    Topol[2] = firstNode + 108 + QTDnodesByLayer;
    Topol[3] = firstNode + 109 + QTDnodesByLayer;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ProductionBC,*Mesh);
    id++;

    Topol[0] = firstNode + 108;
    Topol[1] = firstNode + 119;
    Topol[2] = firstNode + 119 + QTDnodesByLayer;
    Topol[3] = firstNode + 108 + QTDnodesByLayer;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ProductionBC,*Mesh);
    id++;
  }


  Mesh->BuildConnectivity();

#ifdef LOG4CXX
  {
    std::stringstream sout;
    Mesh->Print(sout);
    LOGPZ_DEBUG(logger,sout.str());
  }
#endif

  return Mesh;
}




#include "pzgeoelrefless.h.h"

///CreateGeoElement -> TPZArc3D
template< >
TPZGeoEl *TPZGeoElRefLess<TPZArc3D >::CreateGeoElement(MElementType type, TPZVec<int>& nodeindexes, int matid, int& index)
{
    TPZGeoMesh &mesh = *(this->Mesh());
    if(!&mesh) return 0;
    return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

    #define TPZGEOELEMENTARC3DID 300
    template<>
    int TPZGeoElRefPattern<TPZArc3D>::ClassId() const {
        return TPZGEOELEMENTARC3DID;
    }
    template class 
    TPZRestoreClass< TPZGeoElRefPattern<TPZArc3D>, TPZGEOELEMENTARC3DID>;

    template<>
    TPZCompEl *(*TPZGeoElRefLess<TPZArc3D>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = TPZCompElDisc::CreateDisc;

    template class TPZGeoElRefLess<TPZArc3D>;



///CreateGeoElement -> TPZGeoBlend
#define IMPLEMENTBLEND(TGEO,CLASSID,CREATEFUNCTION) \
template< > \
TPZGeoEl *TPZGeoElRefLess<TPZGeoBlend<TGEO> >::CreateGeoElement(MElementType type, TPZVec<int>& nodeindexes, int matid, int& index) \
{ \
    TPZGeoMesh &mesh = *(this->Mesh()); \
    if(!&mesh) return 0; \
    return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index); \
} \
\
    template<> \
    int TPZGeoElRefPattern<TPZGeoBlend<TGEO>  >::ClassId() const { \
        return CLASSID; \
    } \
    template class \
    TPZRestoreClass< TPZGeoElRefPattern<TPZGeoBlend<TGEO> >, CLASSID>; \
\
    template<> \
    TPZCompEl *(*TPZGeoElRefLess<TPZGeoBlend<TGEO> >::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CREATEFUNCTION; \
\
    template class TPZGeoElRefLess<TPZGeoBlend<TGEO> >;

#define TPZGEOBLENDPOINTID 303
#define TPZGEOBLENDLINEARID 304
#define TPZGEOBLENDQUADID 305
#define TPZGEOBLENDTRIANGLEID 306
#define TPZGEOBLENDCUBEID 307
#define TPZGEOBLENDPRISMID 308
#define TPZGEOBLENDPYRAMIDID 309
#define TPZGEOBLENDTETRAHEDRAID 310

IMPLEMENTBLEND(pzgeom::TPZGeoPoint,TPZGEOBLENDPOINTID,TPZCompElDisc::CreateDisc)
IMPLEMENTBLEND(pzgeom::TPZGeoLinear,TPZGEOBLENDLINEARID,TPZCompElDisc::CreateDisc)
IMPLEMENTBLEND(pzgeom::TPZGeoQuad,TPZGEOBLENDQUADID,TPZCompElDisc::CreateDisc)
IMPLEMENTBLEND(pzgeom::TPZGeoTriangle,TPZGEOBLENDTRIANGLEID,TPZCompElDisc::CreateDisc)
IMPLEMENTBLEND(pzgeom::TPZGeoCube,TPZGEOBLENDCUBEID,TPZCompElDisc::CreateDisc)
IMPLEMENTBLEND(pzgeom::TPZGeoPrism,TPZGEOBLENDPRISMID,TPZCompElDisc::CreateDisc)
IMPLEMENTBLEND(pzgeom::TPZGeoPyramid,TPZGEOBLENDPYRAMIDID,TPZCompElDisc::CreateDisc)
IMPLEMENTBLEND(pzgeom::TPZGeoTetrahedra,TPZGEOBLENDTETRAHEDRAID,TPZCompElDisc::CreateDisc)





#include "pznoderep.h.h"
template class pzgeom::TPZNodeRep<3,TPZArc3D>;
template class pzgeom::TPZNodeRep<8,TPZGeoBlend<TPZGeoCube> >;



