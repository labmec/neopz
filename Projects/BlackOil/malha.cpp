//$Id: malha.cpp,v 1.1 2008-11-12 12:47:09 fortiago Exp $

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
  const REAL PressaoInjecao = 30.e6;
  const REAL PressaoProducao = 25.e6;
  const REAL SoInj = 0.;
  const REAL SoProd = 666; ///a saturacao na producao nao eh prescrita por ser um outflow


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

  TPZAutoPointer<TPZMaterial> mat = new TPZBlackOil2P3D(1, deltaT);

  TPZFMatrix val1(2,2,0.),val2(2,1,0.);
  val2(0,0) = PressaoInjecao;
  val2(1,0) = SoInj;
  TPZAutoPointer<TPZMaterial> bcI = mat->CreateBC(mat,-1, 1,val1,val2);

  val2(0,0) = PressaoProducao;
  val2(1,0) = SoProd;
  TPZAutoPointer<TPZMaterial> bcP = mat->CreateBC(mat,-2, 1,val1,val2);

  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bcI);
  cmesh->InsertMaterialObject(bcP);

  cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetDefaultOrder(0);
  cmesh->AutoBuild();
//   cmesh->AdjustBoundaryElements();
//   cmesh->CleanUpUnconnectedNodes();
//   cmesh->ExpandSolution();

  return cmesh;
}
