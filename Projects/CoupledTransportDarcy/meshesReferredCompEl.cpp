//$Id: meshesReferredCompEl.cpp,v 1.1 2006-06-29 12:28:13 tiago Exp $

#include "meshes.h"
#include "meshesReferredCompEl.h"

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
#include "pzpoisson3dreferred.h"
#include "pztransientmat.h"
#include "pznonlinearpoisson3d.h"

#include "TPZShapeDisc.h"
#include "pzgmesh.h"
#include "pzshapelinear.h"

using namespace pzshape;

void CreateSimpleMeshesWithExactSolutionToReferredCompEl(TPZVec< TPZCompMesh * > & CompMeshes, int h, int p){

  CompMeshes.Resize(2);
  CompMeshes[0] = NULL;
  CompMeshes[1] = NULL;

  TPZCompEl::gOrder = p;

  int nnode = 4;
  REAL co[4][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.}};
  
  int nodesperel = 4;
  int nelem = 1;
  TPZGeoEl *elvec[1]; 
  int indices[1][4] = {{0,1,2,3}};
  
//   int nodesperel = 3;
//   int nelem = 2;
//   TPZGeoEl *elvec[2];
//   int indices[2][3] = {{0,1,2},{0,2,3}};
  
  TPZGeoMesh *gmesh = new TPZGeoMesh();

  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(3);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    coord[2] = co[nod][2];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el=0; el<nelem; el++) {
    TPZManVector<int,4> nodind(nodesperel);
    for(nod=0; nod<nodesperel; nod++) nodind[nod]=indices[el][nod];
    int index;
    if (nodesperel == 3) elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
    if (nodesperel == 4) elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
  }

  TPZManVector<int,2> nodind(2);
  int index;
  nodind[0] = 0; nodind[1] = 1;
  gmesh->CreateGeoElement(EOned, nodind, -1, index);
  nodind[0] = 1; nodind[1] = 2;
  gmesh->CreateGeoElement(EOned, nodind, -1, index);
  nodind[0] = 2; nodind[1] = 3;
  gmesh->CreateGeoElement(EOned, nodind, -1, index);
  nodind[0] = 3; nodind[1] = 0;
  gmesh->CreateGeoElement(EOned, nodind, -1, index);

  gmesh->BuildConnectivity();
  
  //Refinamento uniforme
  for(int ref = 0; ref < h; ref++){
    TPZManVector<TPZGeoEl *> filhos;
    int n = gmesh->NElements();    
    for(int i = 0; i < n; i++){
      TPZGeoEl * gel = gmesh->ElementVec()[i];
      if (gel->Dimension() == 2) gel->Divide(filhos);
    }//for i
  }//ref

  
  TPZCompMesh * firstmesh = new TPZCompMesh(gmesh);
  firstmesh->SetDimModel(2);
  
  TPZMatPoisson3d * FirstMat = new TPZMatPoisson3d(1, 2);
  FirstMat->SetInternalFlux(0.);
  TPZManVector<REAL,2> convdir(2,0.);
  FirstMat->SetParameters(1., 0., convdir);
  
  TPZCompMesh * secondmesh = new TPZCompMesh(gmesh);
  secondmesh->SetDimModel(2);
  
//   TPZMatPoisson3dReferred * SecondMat = new TPZMatPoisson3dReferred(1,2);
  TPZMatPoisson3dReferred * SecondMat = new TPZTransientMaterial< TPZNonLinearPoisson3d > (1,2, 0.);
  SecondMat->SetForcingFunction(Forcing2);
  //Values 0. and convdir are not used by TPZMatPoisson3dReferred. The convection term will
  //be obtained by the previous solution of FirstMat and firstmesh
  SecondMat->SetParameters(1./2., 0., convdir);
  SecondMat->SetAlpha(1./3.);
  SecondMat->SetSD(0.0);
  
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  val2.Zero();
  TPZBndCond *bc1st = FirstMat->CreateBC(-1, 0, val1, val2);
  bc1st->SetForcingFunction(Dirichlet1);
  firstmesh->InsertMaterialObject(FirstMat);
  firstmesh->InsertMaterialObject(bc1st);
  
  val2(0,0) = 1.;
  TPZBndCond *bc2nd = SecondMat->CreateBC(-1, 0, val1, val2);
  secondmesh->InsertMaterialObject(SecondMat);
  secondmesh->InsertMaterialObject(bc2nd);
  
//   gmesh->ResetReference();
//   firstmesh->LoadReferences();
  firstmesh->SetAllCreateFunctionsContinuous();//Referred();
  firstmesh->AutoBuild();
  firstmesh->AdjustBoundaryElements();
  firstmesh->CleanUpUnconnectedNodes();
  firstmesh->ExpandSolution();
  
   gmesh->ResetReference();
//   secondmesh->LoadReferences();
  secondmesh->SetAllCreateFunctionsContinuousReferred();
  secondmesh->AutoBuild();
  secondmesh->AdjustBoundaryElements();
  secondmesh->CleanUpUnconnectedNodes();
  secondmesh->ExpandSolution();  
  
  CompMeshes[0] = firstmesh;
  CompMeshes[1] = secondmesh;
  
  if(0){
    std::ofstream gmeshfile("gmesh.txt");
    gmesh->Print(gmeshfile);
    std::ofstream cmesh1file("cmesh1.txt");
    firstmesh->Print(cmesh1file);
    std::ofstream cmesh2file("cmesh2.txt");
    secondmesh->Print(cmesh2file);
  }
  
  
}//CreateSimpleMeshWithExactSolutionToReferredCompEl

void CreateMesh_ComoPhilippeQuer_Adimensional_Sem_Simetria(TPZVec< TPZCompMesh * > & CompMeshes, int h, int SingH, int p){
  REAL L = 1.;
  const int numero = 4;
  const int nInf = (int)pow(2,numero);
  
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
  
  TPZGeoElBC gbc1(elvec[0],2,-2,*gmesh);
  TPZGeoElBC gbc22(elvec[0],5,-3,*gmesh);
  TPZGeoElBC gbc23(elvec[0],6,-3,*gmesh);
  
//   TPZManVector<int, 1> nodeind(1); nodeind[0] = 2;
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
  
  /** Uniform refinement */
  for(int i = 0; i < h; i++){
    int n = gmesh->NElements();
    for(int j = 0; j < n; j++){
      TPZGeoEl * gel = gmesh->ElementVec()[j];
      if (gel->HasSubElement()) continue;
      if (gel->Dimension() == 2) gel->Divide(filhos);
    }
  }

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetNonSymmetric();
  TPZManVector<REAL,2> convdir(2,0.);
  REAL beta = 0.0;
  mat->SetParameters(K, beta, convdir);
  mat->SetSD(0.);
  
  int nstate = 1;  
  
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc[4];

  /** Condicao Dirichlet */
  //Dirichlet p = pfracture
//   val2(0,0) = pfracture; bc[0] = mat->CreateBC(-1, 0,val1,val2);
    
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
  
   /** Neumann nulo */
  bc[3] = mat->CreateBC(-3, 1, val1, val2);  
 
  cmesh->InsertMaterialObject(mat);
  for(int ii = 0; ii < 4; ii++) cmesh->InsertMaterialObject(bc[ii]);

  cmesh->SetAllCreateFunctionsContinuous();
//   cmesh->SetAllCreateFunctionsDiscontinuous();
  
  cmesh->AutoBuild();
  
  //Refinamento uniforme
//   for(int i = 0; i < h; i++){
//     int n = cmesh->ElementVec().NElements();
//     TPZManVector<int> filhos;
//     for(int j = 0; j < n; j++){
//       TPZCompEl * cel = cmesh->ElementVec()[j];
//       if (!cel) continue;
//       if (cel->Dimension() == 2){
//         cel->Divide(cel->Index(), filhos);
//       }//if
//     }//for j
//   }//uniforme
//   cmesh->ExpandSolution();

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
  
  if (0){
    std::ofstream mesh1file("Acmesh1.txt");
    cmesh->Print(mesh1file);
    std::ofstream gmeshfile("Agmesh.txt");
    gmesh->Print(gmeshfile);
  }
  cmesh->CleanUpUnconnectedNodes();

  if (0){
    std::ofstream mesh1file("Dcmesh1.txt");
    cmesh->Print(mesh1file);
    std::ofstream gmeshfile("Dgmesh.txt");
    gmesh->Print(gmeshfile);
  }


  TPZCompMesh * secondmesh = new TPZCompMesh(gmesh);
  secondmesh->SetDimModel(2);
  
  TPZMatPoisson3dReferred * SecondMat = new TPZTransientMaterial< TPZNonLinearPoisson3d >(1, 2, 0.2);
  //Values 0. and convdir are not used by TPZMatPoisson3dReferred. The convection term will
  //be obtained by the previous solution of FirstMat and firstmesh
  convdir[0] = 1.; convdir[1] = 1.;
  SecondMat->SetParameters(0.*1./2., 1.e-3, convdir);
  SecondMat->SetAlpha(1.);
  SecondMat->SetSD(1.);

  /** Outflow na producao */
  val2.Zero();
  TPZBndCond * bc2nd[4];
  bc2nd[0] = SecondMat->CreateBC(-3, 3,val1,val2);
  /** Temperatura na fratura */
  REAL TEMP = 1.;
  val2(0,0) = TEMP; 
  bc2nd[1] = SecondMat->CreateBC(-1, 0,val1,val2);  
  /** Neumann nulo */
  val2(0,0) = 0.;
  bc2nd[2] = SecondMat->CreateBC(-9, 1, val1, val2);  
  bc2nd[3] = SecondMat->CreateBC(-2, 1, val1, val2);
  
  secondmesh->InsertMaterialObject( SecondMat );
  secondmesh->InsertMaterialObject( bc2nd[0] );
  secondmesh->InsertMaterialObject( bc2nd[1] );
  secondmesh->InsertMaterialObject( bc2nd[2] );
  secondmesh->InsertMaterialObject( bc2nd[3] );
  gmesh->ResetReference();
  secondmesh->SetAllCreateFunctionsContinuousReferred();
  secondmesh->AutoBuild();  
  secondmesh->AdjustBoundaryElements();
  secondmesh->CleanUpUnconnectedNodes();
  secondmesh->ExpandSolution();  
  
  CompMeshes.Resize(2);
  CompMeshes[0] = cmesh;
  CompMeshes[1] = secondmesh;
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

TPZCompMesh * TesteConvectivoPuro(int h, int SingH, int p){
  REAL L = 1.;
  const int numero = 1;
  const int nInf = (int)pow(2,numero);
  
  std::cout << "\n\n**************************************" << std::endl;
  std::cout << "nInf = " << nInf << std::endl;
  std::cout << "**************************************\n\n" << std::endl;
  REAL Infinity = nInf * L;

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
  
  faceind[0] = 0; faceind[1] = 3;
  TPZGeoEl * GelEntrada = gmesh->CreateGeoElement(EOned, faceind, -12, faceindex, 0);
  faceind[0] = 1; faceind[1] = 2;
  TPZGeoEl * GelSaida = gmesh->CreateGeoElement(EOned, faceind, -13, faceindex, 0);

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
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetNonSymmetric();
  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  REAL beta = 1.0;
  mat->SetParameters(0., beta, convdir);
  mat->SetSD(1.);
  
  int nstate = 1;  
  
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc[5];

  val1(0,0) = 0.; val2(0,0) = 0.; 
  bc[0] = mat->CreateBC(-1, 1, val1,val2);
  bc[1] = mat->CreateBC(-2, 1, val1, val2);
  bc[2] = mat->CreateBC(-9, 1, val1, val2);
  val2(0,0) = 3.0;
  bc[3] = mat->CreateBC(-12, 0, val1, val2);
  val2(0,0) = 3.0;
  bc[4] = mat->CreateBC(-13, 3, val1, val2);
 
  cmesh->InsertMaterialObject(mat);
  for(int ii = 0; ii < 5; ii++) cmesh->InsertMaterialObject(bc[ii]);

  cmesh->SetAllCreateFunctionsContinuous();
//   cmesh->SetAllCreateFunctionsDiscontinuous();
  
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
  
  if (0){
    std::ofstream mesh1file("Acmesh1.txt");
    cmesh->Print(mesh1file);
    std::ofstream gmeshfile("Agmesh.txt");
    gmesh->Print(gmeshfile);
  }
  cmesh->CleanUpUnconnectedNodes();

  if (0){
    std::ofstream mesh1file("Dcmesh1.txt");
    cmesh->Print(mesh1file);
    std::ofstream gmeshfile("Dgmesh.txt");
    gmesh->Print(gmeshfile);
  }

  return cmesh;
}

TPZCompMesh * TesteConvectivoPuro2(int h, int p){

  SetPOrder(p);
  
  TPZVec< TPZVec<REAL> > co;
  TPZManVector< TPZVec< int> > indices;
  co.Resize(4);
  for(int i = 0; i < co.NElements(); i++) co[i].Resize(2);
  co[0][0] = 0.;          co[0][1] = 0.;
  co[1][0] = 1.;          co[1][1] = 0.;
  co[2][0] = 1.;          co[2][1] = 1.;
  co[3][0] = 0.;          co[3][1] = 1.;

//1 quadrado
//   indices.Resize( 1 );
//   indices[0].Resize(4);
//   indices[0][0] = 0; indices[0][1] = 1; indices[0][2] = 2; indices[0][3] = 3;

//2 triangulos
  indices.Resize( 2 );
  indices[0].Resize(3);
  indices[0][0] = 0; indices[0][1] = 1; indices[0][2] = 2;
  indices[1].Resize(3);
  indices[1][0] = 0; indices[1][1] = 2; indices[1][2] = 3;
  
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  
  for(int nod = 0; nod < co.NElements(); nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[nodind].Initialize(nod,co[nod],*gmesh);
  }

  TPZVec<TPZGeoEl *> elvec(indices.NElements(), NULL);
  for(int el = 0; el < indices.NElements(); el++) {
    int index;
    if (indices[el].NElements() == 4){
      elvec[el] = gmesh->CreateGeoElement(EQuadrilateral, indices[el],1,index,0);
    }
    if (indices[el].NElements() == 3){
      elvec[el] = gmesh->CreateGeoElement(ETriangle     , indices[el],1,index,0);
    }
  }//for el

  TPZManVector<int, 2> faceind(2); faceind[0] = 0; faceind[1] = 1;
  int faceindex = -1;
  faceind[0] = 0; faceind[1] = 3;
  TPZGeoEl * GelEntrada = gmesh->CreateGeoElement(EOned, faceind, -12, faceindex, 0);
  faceind[0] = 1; faceind[1] = 2;
  TPZGeoEl * GelSaida = gmesh->CreateGeoElement(EOned, faceind, -13, faceindex, 0);

  gmesh->BuildConnectivity();  
  
  TPZVec<TPZGeoEl*> filhos;
  for(int i = 0; i < h; i++){
    int n = gmesh->NElements();
    for(int j = 0; j < n; j++){
      if (gmesh->ElementVec()[j]->Dimension() == 2) gmesh->ElementVec()[j]->Divide(filhos);
    }
  }

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZTransientMaterial< TPZNonLinearPoisson3d >(1, 2, 0.) ;
  dynamic_cast< TPZNonLinearPoisson3d *>(mat)->SetReferred(false);
  mat->SetNonSymmetric();
  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  REAL beta = 1.0;
  mat->SetParameters(0., beta, convdir);
  mat->SetSD(1.);
  
  int nstate = 1;    
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc[2];

  val1(0,0) = 0.;
  val2(0,0) = 3.0;
  bc[0] = mat->CreateBC(-12, 0, val1, val2);
  val2(0,0) = 0.0;
  bc[1] = mat->CreateBC(-13, 3, val1, val2);
 
  cmesh->InsertMaterialObject(mat);
  for(int ii = 0; ii < 2; ii++) cmesh->InsertMaterialObject(bc[ii]);

  cmesh->SetAllCreateFunctionsContinuous();
//   cmesh->SetAllCreateFunctionsDiscontinuous();
  
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  return cmesh;
}
