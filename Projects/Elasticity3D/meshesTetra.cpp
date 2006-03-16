//$Id: meshesTetra.cpp,v 1.1 2006-03-16 13:25:38 tiago Exp $

#include "meshesTetra.h"

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
#include "pzelast3d.h"

void SetInterpOrder(int p){
  TPZCompEl::gOrder = p;
}

TPZCompMesh * VigaEngastadaTetra(int h, int p){
  SetInterpOrder(p);

  const REAL Length = 1.;
  const REAL Base   = 0.2;
  const REAL Height = 0.5; 
//  const REAL STRESS = 100.;
  const REAL EYoung = 205000.;
  const REAL Poisson= 0.0;
  
  REAL co[9][3] = {{0.,0.,0.},    {0.,0., Base},    {0.,Height, Base},    {0., Height, 0.}, 
                   {Length,0.,0.},{Length,0., Base},{Length,Height, Base},{Length, Height, 0.}, 
                   {Length/2., Height/2., Base/2.} };
  int indices[12][4] = {{0,1,3,8}, {1,2,3,8}, {3,6,7,8}, {3,2,6,8}, {0,4,3,8}, {3,4,7,8}, {0,1,4,8}, {1,5,4,8}, {1,2,5,8}, {2,5,6,8},
                        {4,5,7,8}, {5,6,7,8}};
  TPZGeoEl *elvec[1];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 9;
  int nelem = 12;
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
    TPZManVector<int,4> nodind(4);
    for(nod=0; nod < 4; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(ETetraedro,nodind,1,index);
  }

  int index;
  TPZManVector<int,3> bcincid(3);
  bcincid[0] = 0;
  bcincid[1] = 1;
  bcincid[2] = 3;
  gmesh->CreateGeoElement(ETriangle,bcincid,-1,index);
  
  bcincid[0] = 1;
  bcincid[1] = 2;
  bcincid[2] = 3;
  gmesh->CreateGeoElement(ETriangle,bcincid,-1,index);  
  
  bcincid[0] = 4;
  bcincid[1] = 5;
  bcincid[2] = 7;
  gmesh->CreateGeoElement(ETriangle,bcincid,-2,index);
  
  bcincid[0] = 5;
  bcincid[1] = 6;
  bcincid[2] = 7;
  gmesh->CreateGeoElement(ETriangle,bcincid,-2,index);

  gmesh->BuildConnectivity();
  
  TPZVec<TPZGeoEl*> filhos;
  for(int i = 0; i < h; i++){
    int n = gmesh->NElements();
    for(int j = 0; j < n; j++){
      if (gmesh->ElementVec()[j]->Dimension() == 3) gmesh->ElementVec()[j]->Divide(filhos);
    }
  }    
  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(3);
  
  TPZManVector<REAL,3> NullForce(3); 
  NullForce.Fill(0.);
  TPZElasticity3D * mat = new TPZElasticity3D( 1, EYoung, Poisson, NullForce );

  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  val2.Zero();
  TPZBndCond * bcD = mat->CreateBC(-1, 0,val1,val2);
  
  val2(1,0) = 0.; //*STRESS;
  TPZBndCond * bcN = mat->CreateBC(-2, 1,val1,val2);
  bcN->SetForcingFunction(MomentoExtremidadeTetra);
  
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bcD);
  cmesh->InsertMaterialObject(bcN);
  
  //cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetAllCreateFunctionsContinuous();
  
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
//   cmesh->ExpandSolution();

  return cmesh;
}

void MomentoExtremidadeTetra(TPZVec<REAL> &pto, TPZVec<REAL> &force){
  if ( fabs(pto[0] - 1.) > 1.e-8 ) std::cout << "\nTA TUDO ERRADO!!\n";
  REAL Y = pto[1];
  force.Resize(3);
  force.Fill(0.);
  REAL a = 240.;
  REAL H = 0.5;
  force[0] = 2. * a * Y / H;
}

TPZCompMesh * VigaEngastadaForcaVolumeTetra(int h, int p){
  SetInterpOrder(p);

  const REAL Length = 10.;
  const REAL Base   = 0.2;
  const REAL Height = 0.5; 
//  const REAL STRESS = 100.;
  const REAL EYoung = 205000.;
  const REAL Poisson= 0.0;
  
  REAL co[9][3] = {{0.,0.,0.},    {0.,0., Base},    {0.,Height, Base},    {0., Height, 0.}, 
                   {Length,0.,0.},{Length,0., Base},{Length,Height, Base},{Length, Height, 0.}, 
                   {Length/2., Height/2., Base/2.} };
  int indices[12][4] = {{0,1,3,8}, {1,2,3,8}, {3,6,7,8}, {3,2,6,8}, {0,4,3,8}, {3,4,7,8}, {0,1,4,8}, {1,5,4,8}, {1,2,5,8}, {2,5,6,8},
                        {4,5,7,8}, {5,6,7,8}};
  TPZGeoEl *elvec[1];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 9;
  int nelem = 12;
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
    TPZManVector<int,4> nodind(4);
    for(nod=0; nod < 4; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(ETetraedro,nodind,1,index);
  }

  int index;
  TPZManVector<int,3> bcincid(3);
  bcincid[0] = 0;
  bcincid[1] = 1;
  bcincid[2] = 3;
  gmesh->CreateGeoElement(ETriangle,bcincid,-1,index);
  
  bcincid[0] = 1;
  bcincid[1] = 2;
  bcincid[2] = 3;
  gmesh->CreateGeoElement(ETriangle,bcincid,-1,index);  
  
  bcincid[0] = 4;
  bcincid[1] = 5;
  bcincid[2] = 7;
  gmesh->CreateGeoElement(ETriangle,bcincid,-2,index);
  
  bcincid[0] = 5;
  bcincid[1] = 6;
  bcincid[2] = 7;
  gmesh->CreateGeoElement(ETriangle,bcincid,-2,index);
  
  gmesh->BuildConnectivity();
  
  TPZVec<TPZGeoEl*> filhos;
  for(int i = 0; i < h; i++){
    int n = gmesh->NElements();
    for(int j = 0; j < n; j++){
      if (gmesh->ElementVec()[j]->Dimension() == 3) gmesh->ElementVec()[j]->Divide(filhos);
    }
  }    
  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(3);
  
//   TPZManVector<REAL,3> NullForce(3); 
//   NullForce.Fill(0.);
  TPZManVector<REAL,3> gravity(3); 
  gravity[0] = 0.; gravity[1] = 0.; gravity[2] = -98.;
  TPZElasticity3D * mat = new TPZElasticity3D(1, EYoung, Poisson, /*NullForce*/gravity);

  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  val2.Zero();
  TPZBndCond * bcD = mat->CreateBC(-1, 0,val1,val2);
  
  val2(1,0) = 0.; //*STRESS;
  TPZBndCond * bcN = mat->CreateBC(-2, 1,val1,val2);
  
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bcD);
  cmesh->InsertMaterialObject(bcN);
  
  //cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetAllCreateFunctionsContinuous();
  
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
//   cmesh->ExpandSolution();

  return cmesh;
}
