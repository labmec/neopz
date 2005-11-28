//$Id: meshes.cpp,v 1.1 2005-11-28 13:49:47 tiago Exp $

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
#include "TPZShapeDisc.h"
#include "pzgmesh.h"
#include "pzshapelinear.h"

using namespace pzshape;

void SetPOrder(int p){
  TPZCompEl::gOrder = p;
//   std::cout << "POrder = " << p << std::endl;
}

TPZCompMesh *CreateMesh(int h, int p) {
  REAL L = 1.;
  REAL Infinity = 10. * L;
  REAL K = 1000.;
  REAL pfracture = 100.;
  REAL pinfinity = 1.;
  
  SetPOrder(p);
  REAL co[9][2] = {{0., 0.},{L, 0.},{L, L},{0., L},{Infinity, 0.},{Infinity, L},{Infinity, Infinity},{L, Infinity},{0., Infinity}};
  int indices[4][4] = {{0,1,2,3},{1,4,5,2},{2,5,6,7},{3,2,7,8}};
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

  //TPZGeoElBC gbc;

  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom p = pfracture
  TPZGeoElBC gbc2(elvec[1],5,-2,*gmesh); // right
  TPZGeoElBC gbc3(elvec[2],5,-2,*gmesh); // right
  TPZGeoElBC gbc4(elvec[2],6,-2,*gmesh); // top
  TPZGeoElBC gbc5(elvec[3],6,-2,*gmesh); // top

    
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
  TPZBndCond *bc[2];

  //Dirichlet p = pfracture
  val2.Zero();
  val2(0,0) = pfracture;
  bc[0] = mat->CreateBC(-1, 0,val1,val2);
  
  val2(0,0) = pinfinity;
  bc[1] = mat->CreateBC(-2, 0, val1, val2);
  
  cmesh->InsertMaterialObject(mat);
  for(int ii = 0; ii < 2; ii++) cmesh->InsertMaterialObject(bc[ii]);

//   cmesh->SetAllCreateFunctionsDiscontinuous(); 
//   TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);
  cmesh->SetAllCreateFunctionsContinuous();
  
  cmesh->AutoBuild();

  cmesh->AdjustBoundaryElements();
  //  cmesh->CleanUpUnconnectedNodes();
  //  cmesh->ExpandSolution();
  
  return cmesh;
}

TPZCompMesh *CreateMeshSingularity(int h, int SingH, int p) {
  REAL L = 1.;
  REAL Infinity = 10. * L;
  REAL K = 1.;//0.001;
  REAL pfracture = 100.;
  REAL KSUPORTE = 1.e-3;
  REAL pinfinity = 1.;
  
  SetPOrder(p);
  REAL co[9][2] = {{0., 0.},{L, 0.},{L, L},{0., L},{Infinity, 0.},{Infinity, L},{Infinity, Infinity},{L, Infinity},{0., Infinity}};
  int indices[4][4] = {{0,1,2,3},{1,4,5,2},{2,5,6,7},{3,2,7,8}};
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
//   TPZVec<TPZGeoEl*> filhos; 
//   for(int i = 0; i < h; i++){
//      int n = gmesh->NElements();
//      for(int j = 0; j < n; j++){
//        if (gmesh->ElementVec()[j]->Dimension() == 2) gmesh->ElementVec()[j]->Divide(filhos);
//      }
//   }
  //TPZGeoElBC gbc;
  
  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom p = pfracture
  TPZGeoElBC gbc2(elvec[1],5,-2,*gmesh); // right
  TPZGeoElBC gbc3(elvec[2],5,-2,*gmesh); // right
  TPZGeoElBC gbc4(elvec[2],6,-2,*gmesh); // top
  TPZGeoElBC gbc5(elvec[3],6,-2,*gmesh); // top

    
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
  TPZBndCond *bc[2];

  /** Condicao Dirichlet */
  //Dirichlet p = pfracture
  val2(0,0) = pfracture; bc[0] = mat->CreateBC(-1, 0,val1,val2);
    
  /** Condicao mista */
  //val1(0,0) = KSUPORTE; val2(0,0) = pfracture * KSUPORTE; bc[0] = mat->CreateBC(-1, 2,val1,val2);

  
  val2(0,0) = pinfinity;
  bc[1] = mat->CreateBC(-2, 0, val1, val2);
  
  cmesh->InsertMaterialObject(mat);
  for(int ii = 0; ii < 2; ii++) cmesh->InsertMaterialObject(bc[ii]);

  cmesh->SetAllCreateFunctionsDiscontinuous();
//   TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);
//  cmesh->SetAllCreateFunctionsContinuous();
  
  cmesh->AutoBuild();
  
  //Refinamento uniforme
  for(int i = 0; i < h; i++){
    int n = cmesh->ElementVec().NElements();
    TPZManVector<int> filhos;    
    for(int j = 0; j < n; j++){
      TPZCompEl * cel = cmesh->ElementVec()[j];
      if (cel->Dimension() == 2){
        int index = cel->Index();
        cel->Divide(cel->Index(), filhos);
        cmesh->ElementVec()[j] = NULL;
      }//if
    }//for j
  }//uniforme

  for(int i = 0; i < SingH; i++){
    TPZStack<TPZGeoEl*> stack;
    TPZGeoElSide RefEl(elvec[0], 1);
    TPZGeoElSide neigh = RefEl.Neighbour();
    while (neigh != RefEl){
      if (neigh.Element()){
        stack.Push(neigh.Element());
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

      if (!cel->Reference()){
        TPZInterfaceElement *face = dynamic_cast<TPZInterfaceElement*> (cel);
        if (face){
          std::cout << "INTERFACE\n";
        }
        TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc*> (cel);
        if (disc){
          std::cout << "DESCONTINUO\n";
        }

      }
      
      if (cel->Reference()->Dimension() == 2){
        int index = cel->Index();
        cel->Divide(cel->Index(), filhos);
        cmesh->ElementVec()[index] = NULL;
        int nsons = filhos.NElements();
        for(int ison = 0; ison < nsons; ison++){
          TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement*>(cmesh->ElementVec()[ filhos[ison] ]);
          if(intel) intel->SetPreferredSideOrder(p+i+1);
          TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cmesh->ElementVec()[ filhos[ison] ]);
          if (disc) disc->SetDegree(p+i+1);
        }//for son
        cmesh->ExpandSolution();
        
//         {
//           std::cout << "NShapeF = ";
//           TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement*>(cmesh->ElementVec()[ filhos[0] ]);
//           if (intel) std::cout << intel->NShapeF() << std::endl;
//           TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cmesh->ElementVec()[ filhos[0] ]);
//           if (disc) std::cout << disc->NShapeF() << std::endl;
//         }

      }//if
    }//j
  }//i
  
  cmesh->AdjustBoundaryElements();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  return cmesh;
}


TPZCompMesh *CreateMeshSingularityKeepingAspectRatio(int h, int SingH, int p) {
  REAL L = 1.;
  const int nInf = 11;
  REAL Infinity = nInf * L;
  REAL K = 0.001;//1.;
  REAL pfracture = 100.;
  REAL KSUPORTE = 1.e-3;
  REAL pinfinity = 1.;
  
  SetPOrder(p);
  
  TPZVec< TPZVec<REAL> > co;
  TPZVec< TPZVec< int> > indices;
  int nnode = nInf*nInf;
  co.Resize(nInf*nInf);
  int nelem = (nInf-1) * (nInf-1);
  indices.Resize( nelem );
  for(int i = 0; i < nInf; i++){
    for(int j = 0; j < nInf; j++){
      int index = j * nInf + i;
      co[index].Resize(2);
      co[index][0] = L * i;
      co[index][1] = L * j;
    }//j
  }//i
  
  for(int i = 0; i < nelem; i++){
    int j = (int)(i/(nInf-1));
    int col = i%(nInf-1);
    int first = j * (nInf) + col;
    indices[i].Resize(4);
    indices[i][0] = first;
    indices[i][1] = first + 1;
    indices[i][2] = first + nInf + 1;
    indices[i][3] = first + nInf;
    std::cout << "j= " << j << " i = " << i << "nos= " << indices[i][0] << "\t" << indices[i][1] << "\t" << indices[i][2] << "\t" << indices[i][3] << "\n";
  }//i
  
  TPZGeoEl *elvec[nelem];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
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
  
  //Creating Boundary elements
  

  gmesh->BuildConnectivity();
//   TPZVec<TPZGeoEl*> filhos; 
//   for(int i = 0; i < h; i++){
//      int n = gmesh->NElements();
//      for(int j = 0; j < n; j++){
//        if (gmesh->ElementVec()[j]->Dimension() == 2) gmesh->ElementVec()[j]->Divide(filhos);
//      }
//   }
  //TPZGeoElBC gbc;
  
  
  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom p = pfracture
  for(int i = 0; i < nInf - 1; i++){ // right
    int index = i * (nInf - 1) + (nInf - 1 - 1);
    TPZGeoElBC gbc2(elvec[index],5,-2,*gmesh);
  }
  
  for(int i = 0; i < nInf - 1; i++){ // top
    int index = (nInf - 1) * (nInf - 2) + i;
    TPZGeoElBC gbc3(elvec[index],6,-2,*gmesh);
  }
  
  
//   TPZGeoElBC gbc2(elvec[1],5,-2,*gmesh); // right
//   TPZGeoElBC gbc3(elvec[2],5,-2,*gmesh); // right
//   TPZGeoElBC gbc4(elvec[2],6,-2,*gmesh); // top
//   TPZGeoElBC gbc5(elvec[3],6,-2,*gmesh); // top

    
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
  TPZBndCond *bc[2];

  /** Condicao Dirichlet */
  //Dirichlet p = pfracture
  //val2(0,0) = pfracture; bc[0] = mat->CreateBC(-1, 0,val1,val2);
    
  /** Condicao mista */
  val1(0,0) = KSUPORTE; val2(0,0) = pfracture * KSUPORTE; bc[0] = mat->CreateBC(-1, 2,val1,val2);

  
  val2(0,0) = pinfinity;
  bc[1] = mat->CreateBC(-2, 0, val1, val2);
  
  cmesh->InsertMaterialObject(mat);
  for(int ii = 0; ii < 2; ii++) cmesh->InsertMaterialObject(bc[ii]);

//   cmesh->SetAllCreateFunctionsDiscontinuous(); 
//   TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);
  cmesh->SetAllCreateFunctionsContinuous();
  
  cmesh->AutoBuild();
  
  //Refinamento uniforme
  for(int i = 0; i < h; i++){
    int n = cmesh->ElementVec().NElements();
    TPZManVector<int> filhos;    
    for(int j = 0; j < n; j++){
      TPZCompEl * cel = cmesh->ElementVec()[j];
      if (cel->Dimension() == 2){
        cel->Divide(cel->Index(), filhos);
      }//if
    }//for j
  }//uniforme

  for(int i = 0; i < SingH; i++){
    TPZStack<TPZGeoEl*> stack;
    TPZGeoElSide RefEl(elvec[0], 1);
    TPZGeoElSide neigh = RefEl.Neighbour();
    while (neigh != RefEl){
      if (neigh.Element()){
        stack.Push(neigh.Element());
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
          if(intel) intel->SetPreferredSideOrder(p+i+1);
          TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cmesh->ElementVec()[ filhos[ison] ]);
          if (disc) disc->SetDegree(p+i+1);
        }//for son
        
//         {
//           std::cout << "NShapeF = ";
//           TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement*>(cmesh->ElementVec()[ filhos[0] ]);
//           if (intel) std::cout << intel->NShapeF() << std::endl;
//           TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cmesh->ElementVec()[ filhos[0] ]);
//           if (disc) std::cout << disc->NShapeF() << std::endl;
//         }

      }//if
    }//j
  }//i
  
  cmesh->AdjustBoundaryElements();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  return cmesh;
}

TPZCompMesh * CreateMeshSingularityKeepingAspectRatio_TriangularDomain(int h, int SingH, int p) {
  REAL L = 1.;
  const int nInf = 10;
  REAL Infinity = nInf * L;
  REAL K = 0.001;//1.;
  REAL pfracture = 100.;
  REAL KSUPORTE = 1.e-3;
  REAL pinfinity = 0.;
  
  SetPOrder(p);
  
  TPZVec< TPZVec<REAL> > co;
  TPZManVector< TPZVec< int> > indices;
  int nnode = nInf*nInf;
  co.Resize(nInf*nInf);
  int nMaxElem = 2 * (nInf-1) * (nInf-1);
  indices.Resize( nMaxElem );
  for(int i = 0; i < nInf; i++){
    for(int j = 0; j < nInf; j++){
      int index = j * nInf + i;
      co[index].Resize(2);
      co[index][0] = L * i;
      co[index][1] = L * j;
    }//j
  }//i
  
  int nrow = nInf-1;
  int ncol = nInf-1;
  int iel = 0;
  indices.Resize(iel+1);
  TPZVec< int > BCElements(nrow,-1);
  int ibc = 0;
  for(int j = 0; j < nrow; j++){
    for(int i = 0; i < ncol; i++){
      int first = j * nInf + i;
      if (i == ncol - 1){
        indices[iel].Resize(3);
        indices[iel][0] = first;
        indices[iel][1] = first + 1;
        indices[iel][2] = first + nInf;
        BCElements[ibc] = iel;
        ibc++;
        iel++;
      }//last element = triangular
      else{
//         indices[iel].Resize(4);
//         indices[iel][0] = first;
//         indices[iel][1] = first + 1;
//         indices[iel][2] = first + nInf + 1;
//         indices[iel][3] = first + nInf;      
//         iel++;
        indices[iel].Resize(3);
        indices[iel][0] = first;
        indices[iel][1] = first + 1;
        indices[iel][2] = first + nInf;      
        iel++;
        indices[iel].Resize(3);
        indices[iel][0] = first + 1;
        indices[iel][1] = first + nInf + 1;
        indices[iel][2] = first + nInf;      
        iel++;        
      }//other elements
    }//for j
    ncol--;
  }//for i
  const int nelem = iel;
  TPZGeoEl *elvec[nelem];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
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
  
  //Creating Boundary elements
  

  gmesh->BuildConnectivity();

  TPZGeoElBC gbc1(elvec[0],3,-1,*gmesh); // bottom p = pfracture
  for(int ii = 0; ii < BCElements.NElements(); ii++){
    if (elvec[BCElements[ii]]){
      TPZGeoElBC gbc2(elvec[BCElements[ii]],4,-2,*gmesh);
    }
    else{
      std::cout << "ERROR " << __PRETTY_FUNCTION__;
    }
  }//for

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
  TPZBndCond *bc[2];

  /** Condicao Dirichlet */
  //Dirichlet p = pfracture
  //val2(0,0) = pfracture; bc[0] = mat->CreateBC(-1, 0,val1,val2);
    
  /** Condicao mista */
  val1(0,0) = KSUPORTE; val2(0,0) = pfracture * KSUPORTE; bc[0] = mat->CreateBC(-1, 2,val1,val2);

  
  val2(0,0) = pinfinity;
  bc[1] = mat->CreateBC(-2, 0, val1, val2);
  
  cmesh->InsertMaterialObject(mat);
  for(int ii = 0; ii < 2; ii++) cmesh->InsertMaterialObject(bc[ii]);

//   cmesh->SetAllCreateFunctionsDiscontinuous();
//   TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);
  cmesh->SetAllCreateFunctionsContinuous();
  
  cmesh->AutoBuild();
  
  //Refinamento uniforme
  for(int i = 0; i < h; i++){
    int n = cmesh->ElementVec().NElements();
    TPZManVector<int> filhos;    
    for(int j = 0; j < n; j++){
      TPZCompEl * cel = cmesh->ElementVec()[j];
      if (cel->Dimension() == 2){
        cel->Divide(cel->Index(), filhos);
      }//if
    }//for j
  }//uniforme

  for(int i = 0; i < SingH; i++){
    TPZStack<TPZGeoEl*> stack;
    TPZGeoElSide RefEl(elvec[0], 1);
    TPZGeoElSide neigh = RefEl.Neighbour();
    while (neigh != RefEl){
      if (neigh.Element()){
        stack.Push(neigh.Element());
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
  }//i
  
  cmesh->AdjustBoundaryElements();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  return cmesh;
}


TPZCompMesh * CreateMeshSingularityKeepingAspectRatio_TriangularDomain2(int h, int SingH, int p) {
  REAL L = 1.;
  const int nInf = 10;
  REAL Infinity = nInf * L;
  REAL K = 1.;
  REAL pfracture = 100.;
  REAL KSUPORTE = 1.e-3;
  
  SetPOrder(p);
  
  TPZVec< TPZVec<REAL> > co;
  TPZManVector< TPZVec< int> > indices;
  int nnode = nInf*nInf;
  co.Resize(nInf*nInf);
  int nMaxElem = 2 * (nInf-1) * (nInf-1);
  indices.Resize( nMaxElem );
  for(int i = 0; i < nInf; i++){
    for(int j = 0; j < nInf; j++){
      int index = j * nInf + i;
      co[index].Resize(2);
      co[index][0] = L * i;
      co[index][1] = L * j;
    }//j
  }//i
  
  int nrow = nInf-1;
  int ncol = nInf-1;
  int iel = 0;
  TPZVec< TPZGeoEl * > TriangularElements(nrow,NULL);
  for(int j = 0; j < nrow; j++){
    for(int i = 0; i < ncol; i++){
      indices.Resize(iel+1);
      int first = j * nInf + i;
      if (i == ncol - 1){
        indices[iel].Resize(3);
        indices[iel][0] = first;
        indices[iel][1] = first + 1;
        indices[iel][2] = first + nInf;
        iel++;
      }//last element = triangular
      else{
        indices[iel].Resize(4);
        indices[iel][0] = first;
        indices[iel][1] = first + 1;
        indices[iel][2] = first + nInf + 1;
        indices[iel][3] = first + nInf;      
        iel++;
      }//square elements
    }//for j
    ncol--;
  }//for i
  const int nelem = iel;

  TPZGeoEl *elvec[nelem];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  int itriang;
  for(el=0; el<nelem; el++) {
    int index;
    if (indices[el].NElements() == 4){
      elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,indices[el],1,index);
      continue;
    }
    if (indices[el].NElements() == 3){
      elvec[el] = gmesh->CreateGeoElement(ETriangle,indices[el],1,index);
      TriangularElements[itriang] = elvec[el];
      itriang++;
      continue;
    }
    std::cout << "ERROR - " << __PRETTY_FUNCTION__;
  }//for el

  gmesh->BuildConnectivity();

  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom p = pfracture
  for(int ii = 0; ii < TriangularElements.NElements(); ii++){
    if (TriangularElements[ii]){
      TPZGeoElBC gbc2(TriangularElements[ii],4,-2,*gmesh);
    }
    else{
      std::cout << "ERROR " << __PRETTY_FUNCTION__;
    }
  }//for

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
  TPZBndCond *bc[2];

  /** Condicao Dirichlet */
  //Dirichlet p = pfracture
  val2(0,0) = pfracture; bc[0] = mat->CreateBC(-1, 0,val1,val2);
    
  /** Condicao mista */
  //val1(0,0) = KSUPORTE; val2(0,0) = pfracture * KSUPORTE; bc[0] = mat->CreateBC(-1, 2,val1,val2);

  
  val2(0,0) = 0.;
  bc[1] = mat->CreateBC(-2, 0, val1, val2);
  
  cmesh->InsertMaterialObject(mat);
  for(int ii = 0; ii < 2; ii++) cmesh->InsertMaterialObject(bc[ii]);

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
    TPZGeoElSide RefEl(elvec[0], 1);
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
  return cmesh;
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

