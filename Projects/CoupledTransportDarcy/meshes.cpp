//$Id: meshes.cpp,v 1.1 2006-01-10 19:40:31 tiago Exp $

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
#include "TPZInterfaceEl.h"
#include "pzintel.h"
#include "pzbndcond.h"
#include "pzpoisson3d.h"
#include "pzcoupledtransportdarcy.h"
#include "pzcoupledtransportdarcyBC.h"
#include "TPZShapeDisc.h"
#include "pzgmesh.h"
#include "pzshapelinear.h"

using namespace pzshape;

void SetPOrder(int p){
  TPZCompEl::gOrder = p;
//   std::cout << "POrder = " << p << std::endl;
}

TPZCompMesh * CreateMesh_TriangularDomain_ComoPhilippeQuer(int h, int SingH, int p){
  REAL L = 1.;
  const int numero = 4;
  const int nInf = pow(2,numero);
  std::cout << "\n\n**************************************" << std::endl;
  std::cout << "nInf = " << nInf << std::endl;
  std::cout << "**************************************\n\n" << std::endl;
  REAL Infinity = nInf * L;
  REAL K = 1.;
  REAL pfracture = 100.;
  std::cout << "log(nInf) = " << log(nInf) << std::endl;
  REAL KSUPORTE = 1.*(K/L)*log(nInf);
  std::cout << "Ksuporte = " << KSUPORTE << std::endl;
    
  SetPOrder(p);
  
  TPZVec< TPZVec<REAL> > co;
  TPZManVector< TPZVec< int> > indices;
  co.Resize(6);
  for(int i = 0; i < co.NElements(); i++) co[i].Resize(2);
  co[0][0] = 0.;          co[0][1] = 0.;
  co[1][0] = Infinity/2.; co[1][1] = 0.;
  co[2][0] = Infinity/2.; co[2][1] = Infinity/2.;
  co[3][0] = 0.;          co[3][1] = Infinity/2.;
  co[4][0] = Infinity;    co[4][1] = 0.;
  co[5][0] = 0.;          co[5][1] = Infinity;
  indices.Resize( 3 );
  indices[0].Resize(4);
  indices[0][0] = 0; indices[0][1] = 1; indices[0][2] = 2; indices[0][3] = 3;
  indices[1].Resize(3);
  indices[1][0] = 1; indices[1][1] = 4; indices[1][2] = 2;
  indices[2].Resize(3);
  indices[2][0] = 3; indices[2][1] = 2; indices[2][2] = 5;
  
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  
  for(int nod = 0; nod < co.NElements(); nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[nodind].Initialize(nod,co[nod],*gmesh);
  }

  TPZVec<TPZGeoEl *> elvec(indices.NElements(), NULL);
  for(int el = 0; el < indices.NElements(); el++) {
    int index;
    if (indices[el].NElements() == 4){
      elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,indices[el],1,index);
      continue;
    }
    if (indices[el].NElements() == 3){
      elvec[el] = gmesh->CreateGeoElement(ETriangle,indices[el],1,index);
      continue;
    }
    std::cout << "ERROR - " << __PRETTY_FUNCTION__;
  }//for el

  TPZManVector<int, 2> faceind(2); faceind[0] = 0; faceind[1] = 1;
  int faceindex = -1;
  TPZGeoEl * BigFace = gmesh->CreateGeoElement(EOned, faceind, -9, faceindex);
  gmesh->BuildConnectivity();  
  
  TPZManVector<TPZGeoEl *> filhos;
  TPZGeoEl * seed = BigFace;
  for(int i = 0; i < numero - 1; i++){
    seed->Divide(filhos);
    seed = filhos[0];
  }//for i
  seed->SetMaterialId(-1);

  /** Dividing elvec[0] to be consistent to refined face */
  seed = elvec[0];
  for(int i = 0; i < numero - 1; i++){
    seed->Divide(filhos);
    seed = filhos[0];
  }//for
  TPZGeoEl * ElFracture = seed;

//   std::ofstream gmeshfile("gmesh.txt");
//   gmesh->Print(gmeshfile);

  TPZGeoElBC gbc1(elvec[1],4,-2,*gmesh);
  TPZGeoElBC gbc2(elvec[2],4,-2,*gmesh);

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetNonSymmetric();
  TPZManVector<REAL,2> convdir(2,0.);
  REAL beta = 0.0;
  mat->SetParameters(K, beta, convdir);
  int nstate = 1;
  
  
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc[3];

  /** Condicao Dirichlet */
  //Dirichlet p = pfracture
  //val2(0,0) = pfracture; bc[0] = mat->CreateBC(-1, 0,val1,val2);
    
  /** Condicao mista */
  val1(0,0) = KSUPORTE; val2(0,0) = pfracture * KSUPORTE; bc[0] = mat->CreateBC(-1, 2,val1,val2);

  /** Dirichlet nulo na simetria */
  val2(0,0) = 0.;
  bc[1] = mat->CreateBC(-2, 0, val1, val2);

  /** Neumann nulo */
  val2(0,0) = 0.;
  bc[2] = mat->CreateBC(-9, 1, val1, val2);
  
  cmesh->InsertMaterialObject(mat);
  for(int ii = 0; ii < 3; ii++) cmesh->InsertMaterialObject(bc[ii]);

  cmesh->SetAllCreateFunctionsContinuous();
  
  cmesh->AutoBuild();
  
  //Refinamento uniforme
  for(int i = 0; i < h; i++){
    int n = cmesh->ElementVec().NElements();
    TPZManVector<int> filhos;
    for(int j = 0; j < n; j++){
      TPZCompEl * cel = cmesh->ElementVec()[j];
      if (!cel) continue;
      if (cel->Dimension() == 2){
        cel->Divide(cel->Index(), filhos);
      }//if
    }//for j
  }//uniforme
  cmesh->ExpandSolution();

  for(int i = 0; i < SingH; i++){
    cmesh->ExpandSolution();
    TPZStack<TPZGeoEl*> stack;
    TPZGeoElSide RefEl(ElFracture, 1);
    TPZGeoElSide neigh = RefEl.Neighbour();
    while (neigh != RefEl){
      if (neigh.Element()){
        if (neigh.Element()->Dimension() == 2){
          stack.Push(neigh.Element());
        }//if dim = 2
        neigh = neigh.Neighbour();
      }
      else {
        std::cout << "Deu errado, neigh.Element() nao podia ser NULL" << std::endl;
        break;
      }
    }//while
    SetPOrder(p+i+1);
    int n = stack.NElements();
    TPZManVector<int> filhos;
    for(int j = 0; j < n; j++){
      TPZCompEl * cel = stack[j]->Reference();
      if (!cel) continue;
      if (cel->Reference()->Dimension() == 2){
        cel->Divide(cel->Index(), filhos);
        int nsons = filhos.NElements();
        for(int ison = 0; ison < nsons; ison++){
          TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement*>(cmesh->ElementVec()[ filhos[ison] ]);
          if(intel) if (p+i+1 < 10) intel->SetPreferredSideOrder(p+i+1);
          TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cmesh->ElementVec()[ filhos[ison] ]);
          if (disc) if (p+i+1 < 10) disc->SetDegree(p+i+1);
        }//for son
      }//if
    }//j
    cmesh->ExpandSolution();
  }//i

  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();

//   std::ofstream cmeshfile("cmesh.txt");
//   cmesh->Print(cmeshfile);
  return cmesh;
}



/*************************/////

TPZCompMesh * CreateMesh_TriangularDomain_ComoPhilippeQuer_Adimensional(int h, int SingH, int p){
  REAL L = 1.;
  const int numero = 4;
  const int nInf = pow(2,numero);
  std::cout << "\n\n**************************************" << std::endl;
  std::cout << "nInf = " << nInf << std::endl;
  std::cout << "**************************************\n\n" << std::endl;
  REAL Infinity = nInf * L;
  REAL K = 1.;
  REAL pfracture = 1./2.;
  std::cout << "log(nInf) = " << log(nInf) << std::endl;
  REAL KSUPORTE = 0.; //valor variado, sera imposto via funcao
  std::cout << "Ksuporte = " << KSUPORTE << std::endl;
    
  SetPOrder(p);
  
  TPZVec< TPZVec<REAL> > co;
  TPZManVector< TPZVec< int> > indices;
  co.Resize(6);
  for(int i = 0; i < co.NElements(); i++) co[i].Resize(2);
  co[0][0] = 0.;          co[0][1] = 0.;
  co[1][0] = Infinity/2.; co[1][1] = 0.;
  co[2][0] = Infinity/2.; co[2][1] = Infinity/2.;
  co[3][0] = 0.;          co[3][1] = Infinity/2.;
  co[4][0] = Infinity;    co[4][1] = 0.;
  co[5][0] = 0.;          co[5][1] = Infinity;
  indices.Resize( 3 );
  indices[0].Resize(4);
  indices[0][0] = 0; indices[0][1] = 1; indices[0][2] = 2; indices[0][3] = 3;
  indices[1].Resize(3);
  indices[1][0] = 1; indices[1][1] = 4; indices[1][2] = 2;
  indices[2].Resize(3);
  indices[2][0] = 3; indices[2][1] = 2; indices[2][2] = 5;
  
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  
  for(int nod = 0; nod < co.NElements(); nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[nodind].Initialize(nod,co[nod],*gmesh);
  }

  TPZVec<TPZGeoEl *> elvec(indices.NElements(), NULL);
  for(int el = 0; el < indices.NElements(); el++) {
    int index;
    if (indices[el].NElements() == 4){
      elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,indices[el],1,index);
      continue;
    }
    if (indices[el].NElements() == 3){
      elvec[el] = gmesh->CreateGeoElement(ETriangle,indices[el],1,index);
      continue;
    }
    std::cout << "ERROR - " << __PRETTY_FUNCTION__;
  }//for el

  TPZManVector<int, 2> faceind(2); faceind[0] = 0; faceind[1] = 1;
  int faceindex = -1;
  TPZGeoEl * BigFace = gmesh->CreateGeoElement(EOned, faceind, -9, faceindex);
  gmesh->BuildConnectivity();  
  
  TPZManVector<TPZGeoEl *> filhos;
  TPZGeoEl * seed = BigFace;
  for(int i = 0; i < numero - 1; i++){
    seed->Divide(filhos);
    seed = filhos[0];
  }//for i
  seed->SetMaterialId(-1);

  /** Dividing elvec[0] to be consistent to refined face */
  seed = elvec[0];
  for(int i = 0; i < numero - 1; i++){
    seed->Divide(filhos);
    seed = filhos[0];
  }//for
  TPZGeoEl * ElFracture = seed;

//   std::ofstream gmeshfile("gmesh.txt");
//   gmesh->Print(gmeshfile);

  TPZGeoElBC gbc1(elvec[1],4,-2,*gmesh);
  TPZGeoElBC gbc2(elvec[2],4,-2,*gmesh);

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetNonSymmetric();
  TPZManVector<REAL,2> convdir(2,0.);
  REAL beta = 0.0;
  mat->SetParameters(K, beta, convdir);
  
  int nstate = 1;
  
  
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc[3];

  /** Condicao Dirichlet */
  //Dirichlet p = pfracture
  //val2(0,0) = pfracture; bc[0] = mat->CreateBC(-1, 0,val1,val2);
    
  /** Condicao mista */
  val1(0,0) = KSUPORTE; val2(0,0) = pfracture * KSUPORTE; bc[0] = mat->CreateBC(-1, 2,val1,val2);
  bc[0]->SetForcingFunction(KRebocoVal2);
  bc[0]->SetVal1Function( KRebocoVal1 );

  /** Dirichlet nulo na simetria */
  val2(0,0) = 0.;
  bc[1] = mat->CreateBC(-2, 0, val1, val2);

  /** Neumann nulo */
  val2(0,0) = 0.;
  bc[2] = mat->CreateBC(-9, 1, val1, val2);
  
  cmesh->InsertMaterialObject(mat);
  for(int ii = 0; ii < 3; ii++) cmesh->InsertMaterialObject(bc[ii]);

  cmesh->SetAllCreateFunctionsContinuous();
  
  cmesh->AutoBuild();
  
  //Refinamento uniforme
  for(int i = 0; i < h; i++){
    int n = cmesh->ElementVec().NElements();
    TPZManVector<int> filhos;
    for(int j = 0; j < n; j++){
      TPZCompEl * cel = cmesh->ElementVec()[j];
      if (!cel) continue;
      if (cel->Dimension() == 2){
        cel->Divide(cel->Index(), filhos);
      }//if
    }//for j
  }//uniforme
  cmesh->ExpandSolution();

  for(int i = 0; i < SingH; i++){
    cmesh->ExpandSolution();
    TPZStack<TPZGeoEl*> stack;
    TPZGeoElSide RefEl(ElFracture, 1);
    TPZGeoElSide neigh = RefEl.Neighbour();
    while (neigh != RefEl){
      if (neigh.Element()){
        if (neigh.Element()->Dimension() == 2){
          stack.Push(neigh.Element());
        }//if dim = 2
        neigh = neigh.Neighbour();
      }
      else {
        std::cout << "Deu errado, neigh.Element() nao podia ser NULL" << std::endl;
        break;
      }
    }//while
    SetPOrder(p+i+1);
    int n = stack.NElements();
    TPZManVector<int> filhos;
    for(int j = 0; j < n; j++){
      TPZCompEl * cel = stack[j]->Reference();
      if (!cel) continue;
      if (cel->Reference()->Dimension() == 2){
        cel->Divide(cel->Index(), filhos);
        int nsons = filhos.NElements();
        for(int ison = 0; ison < nsons; ison++){
          TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement*>(cmesh->ElementVec()[ filhos[ison] ]);
          if(intel) if (p+i+1 < 10) intel->SetPreferredSideOrder(p+i+1);
          TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cmesh->ElementVec()[ filhos[ison] ]);
          if (disc) if (p+i+1 < 10) disc->SetDegree(p+i+1);
        }//for son
      }//if
    }//j
    cmesh->ExpandSolution();
  }//i

  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();

//   std::ofstream cmeshfile("cmesh.txt");
//   cmesh->Print(cmeshfile);
  return cmesh;
}

TPZCompMesh * CreateMesh_ComoPhilippeQuer_Adimensional_Sem_Simetria(int h, int SingH, int p){
  REAL L = 1.;
  const int numero = 4;
  const int nInf = pow(2,numero);
  std::cout << "\n\n**************************************" << std::endl;
  std::cout << "nInf = " << nInf << std::endl;
  std::cout << "**************************************\n\n" << std::endl;
  REAL Infinity = nInf * L;
  REAL K = 1.;
  REAL pfracture = 1.;
  std::cout << "log(nInf) = " << log(nInf) << std::endl;
  REAL KSUPORTE = 0.; //valor variado, sera imposto via funcao
  std::cout << "Ksuporte = " << KSUPORTE << std::endl;
    
  SetPOrder(p);
  
  TPZVec< TPZVec<REAL> > co;
  TPZManVector< TPZVec< int> > indices;
  co.Resize(4);
  for(int i = 0; i < co.NElements(); i++) co[i].Resize(2);
  co[0][0] = 0.;          co[0][1] = 0.;
  co[1][0] = Infinity;    co[1][1] = 0.;
  co[2][0] = Infinity;    co[2][1] = Infinity;
  co[3][0] = 0.;          co[3][1] = Infinity;

  indices.Resize( 1 );
  indices[0].Resize(4);
  indices[0][0] = 0; indices[0][1] = 1; indices[0][2] = 2; indices[0][3] = 3;
  
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  
  for(int nod = 0; nod < co.NElements(); nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[nodind].Initialize(nod,co[nod],*gmesh);
  }

  TPZVec<TPZGeoEl *> elvec(indices.NElements(), NULL);
  for(int el = 0; el < indices.NElements(); el++) {
    int index;
    if (indices[el].NElements() == 4){
      elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,indices[el],1,index,0);
      continue;
    }
    if (indices[el].NElements() == 3){
      elvec[el] = gmesh->CreateGeoElement(ETriangle,indices[el],1,index,0);
      continue;
    }
    std::cout << "ERROR - " << __PRETTY_FUNCTION__;
  }//for el

  TPZManVector<int, 2> faceind(2); faceind[0] = 0; faceind[1] = 1;
  int faceindex = -1;
  TPZGeoEl * BigFace = gmesh->CreateGeoElement(EOned, faceind, -9, faceindex,0);
  TPZManVector<int, 1> nodeind(1); nodeind[0] = 2;
  
  TPZGeoElBC gbc1(elvec[0],2,-2,*gmesh);
//   int nodeindex = -1;
//   TPZGeoEl * ProductionNode = gmesh->CreateGeoElement(EPoint, nodeind, -2, nodeindex);

  gmesh->BuildConnectivity();  
  
  TPZManVector<TPZGeoEl *> filhos;
  TPZGeoEl * seed = BigFace;
  for(int i = 0; i < numero; i++){
    seed->Divide(filhos);
    seed = filhos[0];
  }//for i
  seed->SetMaterialId(-1);
  TPZGeoEl * ElFracture = seed;
  
  /** Dividing elvec[0] to be consistent to refined face */
  seed = elvec[0];
  for(int i = 0; i < numero; i++){
    seed->Divide(filhos);
    seed = filhos[0];
  }//for


  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZCoupledTransportDarcy *mat;
  mat = new TPZCoupledTransportDarcy(1, 2, 3, 2);
  TPZMatPoisson3d * FirstEq = mat->GetMaterial(0);
  TPZMatPoisson3d * SecondEq = mat->GetMaterial(1);
  FirstEq->SetNonSymmetric();
  SecondEq->SetNonSymmetric();
  
  TPZManVector<REAL,2> convdir(2,0.);
  REAL beta = 0.0;
  FirstEq->SetParameters(K, beta, convdir);
  SecondEq->SetParameters(1., beta, convdir);
  
  int nstate = 1;
  
  
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc[3];

  /** Condicao Dirichlet */
  //Dirichlet p = pfracture
  //val2(0,0) = pfracture; bc[0] = mat->CreateBC(-1, 0,val1,val2);
    
  /** Condicao mista */
  val1(0,0) = KSUPORTE; val2(0,0) = pfracture * KSUPORTE; bc[0] = mat->CreateBC(-1, 2,val1,val2);
  bc[0]->SetForcingFunction(KRebocoVal2);
  bc[0]->SetVal1Function( KRebocoVal1 );

  /** Dirichlet nulo no ponto de producao */
  val2(0,0) = 0.;
  bc[1] = mat->CreateBC(-2, 0, val1, val2);
  
  /** Neumann nulo */
  val2(0,0) = 0.;
  bc[2] = mat->CreateBC(-9, 1, val1, val2);  
 
  cmesh->InsertMaterialObject(mat);
  for(int ii = 0; ii < 3; ii++) cmesh->InsertMaterialObject(bc[ii]);

  cmesh->SetAllCreateFunctionsContinuous();
  //cmesh->SetAllCreateFunctionsDiscontinuous();
  
  cmesh->AutoBuild();
  
  //Refinamento uniforme
  for(int i = 0; i < h; i++){
    int n = cmesh->ElementVec().NElements();
    TPZManVector<int> filhos;
    for(int j = 0; j < n; j++){
      TPZCompEl * cel = cmesh->ElementVec()[j];
      if (!cel) continue;
      if (cel->Dimension() == 2){
        cel->Divide(cel->Index(), filhos);
      }//if
    }//for j
  }//uniforme
  cmesh->ExpandSolution();

  for(int i = 0; i < SingH; i++){
    cmesh->ExpandSolution();
    TPZStack<TPZGeoEl*> stack;
    TPZGeoElSide RefEl(ElFracture, 1);
    TPZGeoElSide neigh = RefEl.Neighbour();
    while (neigh != RefEl){
      if (neigh.Element()){
        if (neigh.Element()->Dimension() == 2){
          stack.Push(neigh.Element());
        }//if dim = 2
        neigh = neigh.Neighbour();
      }
      else {
        std::cout << "Deu errado, neigh.Element() nao podia ser NULL" << std::endl;
        break;
      }
    }//while
    SetPOrder(p+i+1);
    int n = stack.NElements();
    TPZManVector<int> filhos;
    for(int j = 0; j < n; j++){
      TPZCompEl * cel = stack[j]->Reference();
      if (!cel) continue;
      if (cel->Reference()->Dimension() == 2){
        cel->Divide(cel->Index(), filhos);
        int nsons = filhos.NElements();
        for(int ison = 0; ison < nsons; ison++){
          TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement*>(cmesh->ElementVec()[ filhos[ison] ]);
          if(intel) if (p+i+1 < 10) intel->SetPreferredSideOrder(p+i+1);
          TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cmesh->ElementVec()[ filhos[ison] ]);
          if (disc) if (p+i+1 < 10) disc->SetDegree(p+i+1);
        }//for son
      }//if
    }//j
    cmesh->ExpandSolution();
  }//i

  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();

//   std::ofstream cmeshfile("cmesh.txt");
//   cmesh->Print(cmeshfile);
  return cmesh;
}


void KRebocoVal1(TPZVec<REAL> &loc, TPZFMatrix &result){
  //Kreboco = Gamma * L / (3.29*10^-10 * (pfratura - pinf) )
  //Gamma = CL / ( Sqrt[(L/C)^2-(X/L)^2] * (pfratura - pinf) )
  REAL x = loc[0]; //Assume-se que a fratura está alinhada com X e começa na origem (0,0)
  
  REAL CL = 0.01;
//  REAL L = 150.; //metros
  REAL C = 0.0488;
  REAL pfratura = 34.e+6;
  REAL pinf = 27.58e+6;
  REAL DeltaP = pfratura - pinf;
  REAL K = 3.29e-10;

  REAL Kreboco = CL * C / ( DeltaP * K * sqrt(1. - x*x) );
  
//   if (x < 0.5){
//     std::cout << "Para ai";
//   }
  
  result.Redim(1,1);
  result(0,0) = Kreboco;
    
}//method

void KRebocoVal2(TPZVec<REAL> &loc, TPZVec<REAL> &result){
  TPZFMatrix val1(1,1,0.);
  KRebocoVal1(loc, val1);
  result.Resize(1);
  result[0] = val1(0,0) * 1.0;
}

// TPZCompMesh * CreateSimpleMeshWithExactSolution(int h, int p){
// 
//   TPZCompEl::gOrder = p;
// 
//   REAL co[4][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.}};
//   int indices[1][4] = {{0,1,2,3}};
//   TPZGeoEl *elvec[1];
//   TPZGeoMesh *gmesh = new TPZGeoMesh();
//   int nnode = 4;
//   int nelem = 1;
//   int nod;
//   for(nod=0; nod<nnode; nod++) {
//     int nodind = gmesh->NodeVec().AllocateNewElement();
//     TPZManVector<REAL,2> coord(2);
//     coord[0] = co[nod][0];
//     coord[1] = co[nod][1];
//     gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
//   }
// 
//   int el;
//   for(el=0; el<nelem; el++) {
//     TPZManVector<int,4> nodind(4);
//     for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
//     int index;
//     elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
//   }
// 
//   gmesh->BuildConnectivity();
//   
//   TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom
//   TPZGeoElBC gbc2(elvec[0],5,-1,*gmesh); // right
//   TPZGeoElBC gbc3(elvec[0],6,-1,*gmesh); // right
//   TPZGeoElBC gbc4(elvec[0],7,-1,*gmesh); // top
//   
//   //Refinamento uniforme
//   for(int ref = 0; ref < h; ref++){
//     TPZManVector<TPZGeoEl *> filhos;
//     int n = gmesh->NElements();    
//     for(int i = 0; i < n; i++){
//       TPZGeoEl * gel = gmesh->ElementVec()[i];
//       if (gel->Dimension() == 2) gel->Divide(filhos);
//     }//for i
//   }//ref
// 
//   
//   TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
//   cmesh->SetDimModel(2);
//   
//   TPZCoupledTransportDarcy * mat = new TPZCoupledTransportDarcy(1, 2, 3, 2);
//   TPZMatPoisson3d * FirstEq  = mat->GetMaterial(0);
//   TPZMatPoisson3d * SecondEq = mat->GetMaterial(1);
//   FirstEq->SetForcingFunction(Forcing1);
//   SecondEq->SetForcingFunction(Forcing2);
// 
//   TPZManVector<REAL,2> convdir(2,0.);
//   FirstEq->SetParameters(1., 0., convdir);
//   SecondEq->SetParameters(1./2., 0., convdir);
//   mat->SetAlpha(1./3.);
//   
//   int nstate = 1;
//   TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
//   TPZBndCond *bc;  
//   val2.Zero();
//   TPZCoupledTransportDarcyBC * BC = mat->CreateBC(-1);
//   
//   bc = FirstEq->CreateBC(-2, 0,val1,val2);
//   BC->SetMaterial(0, bc);
//   
//   val2(0,0) = 1.;
//   bc = SecondEq->CreateBC(-3, 0, val1, val2);
//   BC->SetMaterial(1, bc);
// 
//   cmesh->InsertMaterialObject(mat);
//   cmesh->InsertMaterialObject(BC);
// 
//   cmesh->SetAllCreateFunctionsContinuous();  
//   cmesh->AutoBuild();
//   cmesh->AdjustBoundaryElements();
//   cmesh->CleanUpUnconnectedNodes();
//   cmesh->ExpandSolution();
//   
//   return cmesh;
// }//end of method

TPZCompMesh * CreateSimpleMeshWithExactSolution(int h, int p){

  TPZCompEl::gOrder = p;

  REAL co[4][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.}};
  int indices[1][4] = {{0,1,2,3}};
  TPZGeoEl *elvec[1];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 4;
  int nelem = 1;
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
  
  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],5,-1,*gmesh); // right
  TPZGeoElBC gbc3(elvec[0],6,-1,*gmesh); // right
  TPZGeoElBC gbc4(elvec[0],7,-1,*gmesh); // top
  
  //Refinamento uniforme
  for(int ref = 0; ref < h; ref++){
    TPZManVector<TPZGeoEl *> filhos;
    int n = gmesh->NElements();    
    for(int i = 0; i < n; i++){
      TPZGeoEl * gel = gmesh->ElementVec()[i];
      if (gel->Dimension() == 2) gel->Divide(filhos);
    }//for i
  }//ref

  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZCoupledTransportDarcy * mat = new TPZCoupledTransportDarcy(1, 2, 3, 2);
  TPZMatPoisson3d * FirstEq  = mat->GetMaterial(0);
  TPZMatPoisson3d * SecondEq = mat->GetMaterial(1);
  FirstEq->SetForcingFunction(Forcing1);
  SecondEq->SetForcingFunction(Forcing2);

  TPZManVector<REAL,2> convdir(2,0.);
  FirstEq->SetParameters(1., 0., convdir);
  SecondEq->SetParameters(1./2., 0., convdir);
  mat->SetAlpha(1./3.);
  
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc;  
  val2.Zero();
  TPZCoupledTransportDarcyBC * BC = mat->CreateBC(-1);
  
  bc = FirstEq->CreateBC(-2, 0,val1,val2);
  bc->SetForcingFunction(Dirichlet1);
  BC->SetMaterial(0, bc);
  
  val2(0,0) = 1.;
  bc = SecondEq->CreateBC(-3, 0, val1, val2);
  BC->SetMaterial(1, bc);

  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(BC);

  cmesh->SetAllCreateFunctionsContinuous();  
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  cmesh->ExpandSolution();
  
  return cmesh;
}//end of method

void ExactSol_p(TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix &deriv){
  REAL x = pt[0];
  REAL y = pt[1];
  sol.Resize(1);
  sol[0] = x*y*(x-1.)*(y-1.);
  deriv.Resize(2,1);
  deriv(0,0) = (-1.+x)*(-1.+y)*y+x*(-1.+y)*y;
  deriv(1,0) = (-1.+x)*x*(-1.+y)+(-1.+x)*x*y;
  
  /************/
  sol[0] = x*y;
  deriv(0,0) = y;
  deriv(1,0) = x;  
}

void ExactSol_u(TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix &deriv){
  REAL x = pt[0];
  REAL y = pt[1];
  sol.Resize(1);
  sol[0] = x*y*(x-1.)*(y-1.)+1.;
  deriv.Resize(2,1);
  deriv(0,0) = (-1.+x)*(-1.+y)*y+x*(-1.+y)*y;
  deriv(1,0) = (-1.+x)*x*(-1.+y)+(-1.+x)*x*y;
  
  /*********/
  REAL val = x*y*(x - 1.)*(y - 1.) + 1.;
  sol[0] = val;
  val = (-1.+x)*(-1.+y)*y+(-1.+y)*x*y;
  deriv(0,0) = val;
  val = (-1.+x)*(-1.+y)*x+(-1.+x)*x*y;
  deriv(1,0) = val;  
}

void Dirichlet1(TPZVec<REAL> &pt, TPZVec<REAL> &force){
  REAL x = pt[0];
  REAL y = pt[1];
  force.Resize(1);
  force[0] = x*y;
}
  
void Forcing1(TPZVec<REAL> &pt, TPZVec<REAL> &force){
  REAL x = pt[0];
  REAL y = pt[1];
  force.Resize(1);
  force[0] = -2.*(-1.+x)*x-2.*(-1.+y)*y;
  force[0] *= -1.;

  /**************/  
  force[0] = 0.;
}

void Forcing2(TPZVec<REAL> &pt, TPZVec<REAL> &force){
  REAL x = pt[0];
  REAL y = pt[1];
  force.Resize(1);
  force[0] = 0.5*(-2.*(-1.+x)*x-2.*(-1.+y)*y) 
             -(1./3.)*((-1.+x)*(-1.+y)*y+x*(-1.+y)*y) * ((-1.+x)*(-1.+y)*y+x*(-1.+y)*y);
  force[0] *= -1.;
  
  /**************/
  REAL val = 1./2.*(-2.*(-1.+x)*x-2.*(-1.+y)*y)
            -1./3.*y*((-1.+x)*(-1.+y)*y+x*(-1.+y)*y);
  force[0] = -1. * val;
}

