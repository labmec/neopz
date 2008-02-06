//Classes utilit�ias
#include "pzvec.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzfmatrix.h"

//Bibliotecas std, math etc
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "TPZCompElDisc.h"
#include "TPZShapeDisc.h"
#include "pzl2projection.h"
#include "tpzautopointer.h"
#include "pzanalysis.h"
#include "TExtFunction.h"
using namespace std;
using namespace pzshape;

int main(){
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  REAL co[4][2] = {{-0.5,-0.2},{2.,-0.2},{2.,3.},{-0.5,3.}};
  for(int nod = 0; nod < 4; nod++){
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }
  cout << "\nTriangulo (0) ou Quadrado (1)? ";
  int opcao;
  cin >> opcao;
  int nnodes;
  if (opcao == 0) nnodes = 3;
  if (opcao == 1) nnodes = 4;
  TPZVec<int> nodind(nnodes);
  for(int i = 0; i < nnodes; i++) nodind[i] = i;
  int index;
  TPZGeoEl * geo;
  if (nnodes == 3) geo = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
  if (nnodes == 4) geo = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
  gmesh->BuildConnectivity();
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  TPZVec<REAL> sol(1,0.);
  TPZAutoPointer<TPZMaterial> material = new TPZL2Projection(1,2,1,sol);
  cmesh->InsertMaterialObject(material);
  TPZCompMesh::SetAllCreateFunctionsDiscontinuous();
  TPZCompEl::SetgOrder(1);
  cmesh->SetDefaultOrder(1);
  cmesh->AutoBuild();
  cout << "\ncmesh->NElements() = " << cmesh->NElements() << "\n";
  TPZCompElDisc * cel = dynamic_cast<TPZCompElDisc*>(cmesh->ElementVec()[0]);
  TPZAutoPointer<TPZFunction> extShapes = new TExtFunction();
  cel->SetExternalShapeFunction(extShapes);  
  TPZVec<int> ord(3,10);
  cel->GetIntegrationRule().SetOrder(ord);
  const int npoints = cel->GetIntegrationRule().NPoints();
  TPZVec<REAL> qsi(2);
  REAL w;
  TPZFMatrix phi, dphi;
  ofstream saidaX("saidaX.txt");
  ofstream saidaY("saidaY.txt");
  saidaX << "{";
  saidaY << "{";
  for(int i = 0; i < npoints; i++){
    cel->GetIntegrationRule().Point(i, qsi, w);
    cel->Shape(qsi, phi, dphi);
    TPZVec<REAL> X(3);
    cel->Reference()->X(qsi,X);
    int pos = phi.Rows()-1;
    saidaX << "{" << X[0] << "," << phi(pos,0) << "},";
    saidaY << "{" << X[1] << "," << phi(pos,0) << "},";
  }
  bool dx = true;
  if (dx){
    cmesh->ExpandSolution();
    cmesh->Solution().Zero();
    int pos = phi.Rows()-1;
    TPZConnect & np = cel->Connect(0);
    int blocknumber = np.SequenceNumber();
    int firsteq = cmesh->Block().Position(blocknumber);
    int ndf = cmesh->Block().Size(blocknumber);
    cmesh->Solution()(firsteq+ndf-1,0) = 1.0;
    
    TPZAnalysis an(cmesh);
    an.Solution() = cmesh->Solution();
    TPZVec<char *> scalnames(1);
    scalnames[0] = "Solution";
    TPZVec<char *> vecnames(0);
    char plotfile[] = "singular.dx";
    an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
    an.PostProcess(8);  
  }
  delete cmesh;
  delete gmesh;
  return 0;
}

