//$Id: meshes.cpp,v 1.2 2006-03-16 13:26:01 tiago Exp $

#include "meshes.h"

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

void SetPOrder(int p){
  TPZCompEl::gOrder = p;
}

TPZCompMesh *BarraTracionada(int h, int p){
  SetPOrder(p);

  const REAL Length = 10.;
  const REAL Base   = 0.3;
  const REAL Height = 0.5; 
  const REAL STRESS = 100.;
  const REAL EYoung = 205000.;
  const REAL Poisson= 0.0;
  
  REAL co[8][3] = {{0.,0.,0.},    {0.,0., Base},    {0.,Height, Base},    {0., Height, 0.}, 
                   {Length,0.,0.},{Length,0., Base},{Length,Height, Base},{Length, Height, 0.}};
  int indices[1][8] = {{0,1,2,3,4,5,6,7}};
  TPZGeoEl *elvec[1];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 8;
  int nelem = 1;
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
      if (gmesh->ElementVec()[j]->Dimension() == 3) gmesh->ElementVec()[j]->Divide(filhos);
    }
  }    
  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(3);
  
  TPZManVector<REAL,3> NullForce(3); 
  NullForce.Fill(0.);
  TPZElasticity3D * mat = new TPZElasticity3D(1, EYoung, Poisson, NullForce);

  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  val2.Zero();
  TPZBndCond * bcD = mat->CreateBC(-1, 0,val1,val2);
  
  val2(0,0) = STRESS;
  TPZBndCond * bcN = mat->CreateBC(-2, 1,val1,val2);
  
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bcD);
  cmesh->InsertMaterialObject(bcN);
  
  //cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetAllCreateFunctionsContinuous();
  
  cmesh->AutoBuild();
//   cmesh->AdjustBoundaryElements();
//   cmesh->CleanUpUnconnectedNodes();
//   cmesh->ExpandSolution();

  return cmesh;
}

TPZCompMesh *BarraTracionadaGirada(int h, int p){
  SetPOrder(p);

  const REAL Length = 10.;
  const REAL Base   = 0.3;
  const REAL Height = 0.5; 
  const REAL STRESS = 100.;
  const REAL EYoung = 205000.;
  const REAL Poisson= 0.0;

  MElementType tipovolume;
#define CUBO
#ifdef CUBO    
  tipovolume = ECube;
  REAL co[8][3] = {{0.,1.92711,1.13412},    {0.,-0.248971, -2.22216},    {0., -1.92711, -1.13412 },    {0., 0.248971, 2.22216 }, 
                  {Length,1.92711,1.13412},    {Length,-0.248971, -2.22216},    {Length, -1.92711, -1.13412 },    {Length, 0.248971, 2.22216 } };
  int indices[1][8] = {{0,1,2,3,4,5,6,7}};
  int nnode = 8;
  int nnodesperel = 8;
#endif

//#define PRISMA
#ifdef PRISMA
  tipovolume = EPrisma;
  REAL co[6][3] = {{0.,1.92711,1.13412},    {0.,-0.248971, -2.22216},    {0., -1.92711, -1.13412 },
                  {Length,1.92711,1.13412},    {Length,-0.248971, -2.22216},    {Length, -1.92711, -1.13412 } };
  int indices[1][6] = {{0,1,2,3,4,5}};
  int nnode = 6;
  int nnodesperel = 6;
#endif
      
  TPZGeoEl *elvec[1];
  int nelem = 1;
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
    TPZManVector<int,8> nodind(nnodesperel);
    for(nod=0; nod < nnodesperel; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(tipovolume,nodind,1,index);
  }

  int index;
#ifdef CUBO  
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
#endif

#ifdef PRISMA
  TPZManVector<int,4> bcincid(3);
  bcincid[0] = 0;
  bcincid[1] = 1;
  bcincid[2] = 2;
  gmesh->CreateGeoElement(ETriangle,bcincid,-1,index);
  
  bcincid[0] = 3;
  bcincid[1] = 4;
  bcincid[2] = 5;
  gmesh->CreateGeoElement(ETriangle,bcincid,-2,index);
#endif
  
  gmesh->BuildConnectivity();

#define GeoRefine
#ifdef GeoRefine    
  TPZVec<TPZGeoEl*> filhos;
  for(int i = 0; i < h; i++){
    int n = gmesh->NElements();
    for(int j = 0; j < n; j++){
      if (gmesh->ElementVec()[j]->Dimension() == 3) gmesh->ElementVec()[j]->Divide(filhos);
    }
  }
#endif
  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(3);
  
  TPZManVector<REAL,3> NullForce(3); 
  NullForce.Fill(0.);
  TPZElasticity3D * mat = new TPZElasticity3D(1, EYoung, Poisson, NullForce);

  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  val2.Zero();
  TPZBndCond * bcD = mat->CreateBC(-1, 0,val1,val2);
  
  val2(0,0) = STRESS;
  TPZBndCond * bcN = mat->CreateBC(-2, 1,val1,val2);
  
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bcD);
  cmesh->InsertMaterialObject(bcN);
  
  //cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetAllCreateFunctionsContinuous();
  
  cmesh->AutoBuild();
  
//#define CompRefine
#ifdef CompRefine
  //Refinamento uniforme
  for(int i = 0; i < h; i++){
    int n = cmesh->ElementVec().NElements();
    TPZManVector<int> filhos;
    for(int j = 0; j < n; j++){
      TPZCompEl * cel = cmesh->ElementVec()[j];
      if (!cel) continue;
      if (cel->Dimension() == 3){
        cel->Divide(cel->Index(), filhos);
      }//if
    }//for j
  }//uniforme
  cmesh->ExpandSolution();
#endif  
  
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
//   cmesh->ExpandSolution();

  return cmesh;
}

TPZCompMesh * VigaEngastada(int h, int p){
  SetPOrder(p);

  const REAL Length = 1.;
  const REAL Base   = 0.2;
  const REAL Height = 0.5; 
  const REAL STRESS = 100.;
  const REAL EYoung = 205000.;
  const REAL Poisson= 0.0;
  
  REAL co[8][3] = {{0.,0.,0.},    {0.,0., Base},    {0.,Height, Base},    {0., Height, 0.}, 
                   {Length,0.,0.},{Length,0., Base},{Length,Height, Base},{Length, Height, 0.}};
  int indices[1][8] = {{0,1,2,3,4,5,6,7}};
  TPZGeoEl *elvec[1];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 8;
  int nelem = 1;
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
  bcN->SetForcingFunction(MomentoExtremidade);
  
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

void MomentoExtremidade(TPZVec<REAL> &pto, TPZVec<REAL> &force){
  if ( fabs(pto[0] - 1.) > 1.e-8 ) std::cout << "\nTA TUDO ERRADO!!\n";
  REAL Y = pto[1];
  force.Resize(3);
  force.Fill(0.);
  REAL a = 240.;
  REAL H = 0.5;
  force[0] = 2. * a * Y / H;
}

TPZCompMesh * VigaEngastadaForcaVolume(int h, int p){
  SetPOrder(p);

  const REAL Length = 10.;
  const REAL Base   = 0.2;
  const REAL Height = 0.5; 
  const REAL STRESS = 100.;
  const REAL EYoung = 205000.;
  const REAL Poisson= 0.0;
  
  REAL co[8][3] = {{0.,0.,0.},    {0.,0., Base},    {0.,Height, Base},    {0., Height, 0.}, 
                   {Length,0.,0.},{Length,0., Base},{Length,Height, Base},{Length, Height, 0.}};
  int indices[1][8] = {{0,1,2,3,4,5,6,7}};
  TPZGeoEl *elvec[1];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 8;
  int nelem = 1;
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

TPZCompMesh *BarraTracionadaNeumann(int h, int p){
  SetPOrder(p);

  const REAL Length = 10.;
  const REAL Base   = 0.3;
  const REAL Height = 0.5; 
  const REAL STRESS = 100.;
  const REAL EYoung = 205000.;
  const REAL Poisson= 0.0;
  
  REAL co[8][3] = {{0.,0.,0.},    {0.,0., Base},    {0.,Height, Base},    {0., Height, 0.}, 
                   {Length,0.,0.},{Length,0., Base},{Length,Height, Base},{Length, Height, 0.}};
  int indices[1][8] = {{0,1,2,3,4,5,6,7}};
  TPZGeoEl *elvec[1];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 8;
  int nelem = 1;
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
      if (gmesh->ElementVec()[j]->Dimension() == 3) gmesh->ElementVec()[j]->Divide(filhos);
    }
  }    
  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(3);
  
  TPZManVector<REAL,3> NullForce(3); 
  NullForce.Fill(0.);
  TPZElasticity3D * mat = new TPZElasticity3D(1, EYoung, Poisson, NullForce);

  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  val2.Zero();
  val2(0,0) = -1.*STRESS;
  TPZBndCond * bcD = mat->CreateBC(-1, 1,val1,val2);
  
  val2(0,0) = STRESS;
  TPZBndCond * bcN = mat->CreateBC(-2, 1,val1,val2);
  
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bcD);
  cmesh->InsertMaterialObject(bcN);
  
  //cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetAllCreateFunctionsContinuous();
  
  cmesh->AutoBuild();
//   cmesh->AdjustBoundaryElements();
//   cmesh->CleanUpUnconnectedNodes();
//   cmesh->ExpandSolution();

  return cmesh;
}

#include "pzpoisson3d.h"
TPZCompMesh * Small(int h, int p){

  SetPOrder(p);
  
  REAL co[9][2] = {{0.5,0.5},{0.5,0.},{1.,0.},{1.,0.5},{1.,1.},{0.5,1.},{0.,1.},{0.,0.5},{0.,0.}};
  int indices[4][4] = {{1,2,3,0},{0,3,4,5},{7,0,5,6},{8,1,0,7}};
  TPZGeoEl *elvec[4];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 9;
  int nelem = 4;
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el=0; el<nelem; el++) {
    TPZManVector<int,4> nodind(4);
    for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
  }

  gmesh->BuildConnectivity();
  
  TPZVec<TPZGeoEl*> filhos;
  for(int i = 0; i < h; i++){
    int n = gmesh->NElements();
    for(int j = 0; j < n; j++){
      if (gmesh->ElementVec()[j]->Dimension() == 2) gmesh->ElementVec()[j]->Divide(filhos);
    }
  }      

  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],5,-3,*gmesh); // right
  TPZGeoElBC gbc3(elvec[1],5,-3,*gmesh); // right
  TPZGeoElBC gbc4(elvec[1],6,-3,*gmesh); // top
  TPZGeoElBC gbc5(elvec[2],6,-3,*gmesh); // top
  TPZGeoElBC gbc6(elvec[2],7,-2,*gmesh); // left
  TPZGeoElBC gbc7(elvec[3],7,-2,*gmesh); // left
  TPZGeoElBC gbc8(elvec[3],4,-1,*gmesh); // bottom

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);

  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);

  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;

  REAL beta = 0.0;
  mat->SetParameters(0.01, beta, convdir);
  mat->SetInternalFlux(10.0);
  int nstate = 1;
  
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc[3];
  val2.Zero();
  bc[0] = mat->CreateBC(-3, 0,val1,val2);
  
  bc[1] = mat->CreateBC(-1, 0, val1, val2);
  
  bc[2] = mat->CreateBC(-2, 0, val1, val2);
  
  cmesh->InsertMaterialObject(mat);
  for(int ii = 0; ii < 3; ii++) cmesh->InsertMaterialObject(bc[ii]);

  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  //  cmesh->CleanUpUnconnectedNodes();
  //  cmesh->ExpandSolution();
  
  return cmesh;
}//Small
