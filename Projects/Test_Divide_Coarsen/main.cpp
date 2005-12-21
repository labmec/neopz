//$Id: main.cpp,v 1.2 2005-12-21 12:00:31 tiago Exp $

/**
 * This program tests the methods TPZCompEl::Coarsen and TPZCompEl::Divide acting on TPZInterpolatedElement and acting on TPZCompElDisc
 * December 07, 2005
 */


#include "pzvec.h"

#include "pzcmesh.h"

#include "pzdebug.h"
#include "pzcheckgeom.h"
//#include "pzerror.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"

//#include "pzintel.h"
#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"

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
#include "pzonedref.h"

#include "pzadmchunk.h"


#include "pzbndcond.h"
#include "pzpoisson3d.h"

#include "pzvisualmatrix.h"

#include <time.h>
#include <stdio.h>

#include "gmres.h"
#include "TPZTimer.h"
using namespace std;

enum TIPO{ Continuous = 0, Discontinuous = 1, Both = 2};

void SetMaterial(TPZCompMesh & cmesh){

  cmesh.SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetNonSymmetric();
  TPZManVector<REAL,2> convdir(2,0.);
  REAL beta = 0.0;
  mat->SetParameters(1., beta, convdir);
  int nstate = 1;  
  
  /** Condicao Dirichlet */  
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc = mat->CreateBC(-1, 0, val1, val2);
  
  cmesh.InsertMaterialObject(mat);
  cmesh.InsertMaterialObject(bc);
}//method


TPZCompMesh * mesh1Elemento(TIPO tipo, const int ndivision, const int ncoarsen){

  TPZCompEl::gOrder = 2;

  const int nnodes = 4;
  const int nelem  = 1;
  TPZVec < TPZVec< REAL > > co(nnodes);
  for(int i = 0; i < co.NElements(); i++) co[i].Resize(2);
  
  co[0][0] = 0.; co[0][1] = 0.;
  co[1][0] = 1.; co[1][1] = 0.;
  co[2][0] = 1.; co[2][1] = 1.;
  co[3][0] = 0.; co[3][1] = 1.;

  
  TPZVec< TPZVec<int> > indices(nelem);
  for(int i = 0; i < indices.NElements(); i++) indices[i].Resize(4);
  
  indices[0][0] = 0; indices[0][1] = 1; indices[0][2] = 2; indices[0][3] = 3;
  
  TPZVec< TPZGeoEl* > elvec(nelem);
  
  TPZGeoMesh *gmesh = new TPZGeoMesh();

  int nod;
  for(nod = 0; nod < nnodes; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el = 0; el < nelem; el++) {
    TPZManVector<int,4> nodind(4);
    for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
  }

   gmesh->BuildConnectivity();
   TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh);

//FAZENDO CONTINUO DESCONTINUO
  TPZVec<TPZGeoEl*> continuous, discontinuous;
  if (tipo == Continuous){
    continuous.Resize(1);
    discontinuous.Resize(0);
    continuous[0] = elvec[0];
  }//if Continuous
  
  if (tipo == Discontinuous){
    continuous.Resize(0);
    discontinuous.Resize(1);
    discontinuous[0] = elvec[0];
  }//if Continuous

  if (tipo == Both){
    continuous.Resize(1);
    discontinuous.Resize(1);
    continuous[0] = elvec[0];
    discontinuous[0] = elvec[1];
  }//if Continuous  
  
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  SetMaterial(*cmesh);
  
  cmesh->AutoBuildContDisc(continuous, discontinuous);
cmesh->SetAllCreateFunctionsDiscontinuous();
//cmesh->AutoBuild();
//   TPZInterfaceElement::SetCalcStiffContDisc();
  
  if (0) {
    ofstream cmeshout("Original_cmesh.txt");
    cmesh->Print(cmeshout);
    ofstream gmeshout("Original_gmesh.txt");
    gmesh->Print(gmeshout);
  }
  
  /** Testing divide method */  
{//Continuous elements  
  set<int> newindices, oldindices;
  for(int i = 0; i < continuous.NElements(); i++) oldindices.insert(continuous[i]->Index());
  
  for(int j = 0; j < ndivision; j++){
    newindices.clear();
    set<int>::iterator e, w;
    e = oldindices.end();
    for(w = oldindices.begin(); w != e; w++){
      const int index = *w;
      TPZManVector<int,8> subindex; 
      cmesh->Divide(index, subindex, 1);
      for(int subi = 0; subi < subindex.NElements(); subi++){
        newindices.insert( subindex[subi] );
      }//for subi
    }//for w
    oldindices.clear();
    e = newindices.end();
    for(w = newindices.begin(); w != e; w++) oldindices.insert(*w);
    cmesh->ExpandSolution();
//     cmesh->AdjustInterfaceElements();
//     cmesh->AdjustBoundaryElements();
//     cmesh->CleanUpUnconnectedNodes();
  }//for j
}//continuous

{//Discontinuous elements  
  set<int> newindices, oldindices;
  for(int i = 0; i < discontinuous.NElements(); i++) oldindices.insert(discontinuous[i]->Index());
  
  for(int j = 0; j < ndivision; j++){
    newindices.clear();
    set<int>::iterator e, w;
    e = oldindices.end();
    for(w = oldindices.begin(); w != e; w++){
      const int index = *w;
      TPZManVector<int,8> subindex; 
      cmesh->Divide(index, subindex, 1);
      for(int subi = 0; subi < subindex.NElements(); subi++){
        newindices.insert( subindex[subi] );
      }//for subi
    }//for w
    oldindices.clear();
    e = newindices.end();
    for(w = newindices.begin(); w != e; w++) oldindices.insert(*w);
    cmesh->ExpandSolution();
    cmesh->RemakeAllInterfaceElements();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();    
  }//for j
}//discontinuous
  
  if (0) {
    ofstream cmeshout("Refined_cmesh.txt");
    cmesh->Print(cmeshout);
    ofstream gmeshout("Refined_gmesh.txt");
    gmesh->Print(gmeshout);
  }

/** Testing Coarsen method */
  for(int icoarsen = 0; icoarsen < ncoarsen; icoarsen++){
    const int n = cmesh->NElements();
    
    /** counting number of non NULL elements */
    int nactualel = 0;
    for(int iactual = 0; iactual < n; iactual++){
      if (cmesh->ElementVec()[iactual]) nactualel++;
    }
    
    
    for(int iel = 0; iel < n; iel++){
      if (!cmesh->ElementVec()[iel]) continue;
      TPZGeoEl * father = cmesh->ElementVec()[iel]->Reference()->Father();
      if (!father) continue;
      int nsub = father->NSubElements();
      TPZManVector<int,4> subindices(nsub);
      bool CreateDiscontinuous = false;
      for(int isub = 0; isub < nsub; isub++){
        TPZGeoEl * subel = father->SubElement(isub);
        if (!subel) continue;
        TPZCompEl * cel = subel->Reference();
        if (!cel) continue;
        int index = cel->Index();
        subindices[isub] = index;        
        TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc *>(cel);
        if (disc) CreateDiscontinuous = true;
      }//for isub
      int coarseindex = -1;      
      cmesh->Coarsen(subindices, coarseindex, CreateDiscontinuous);
    }//for iel  

    cmesh->RemakeAllInterfaceElements();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

  }//for icoarsen
  
    /** counting number of non NULL elements */
    int nactualel = 0;
    for(int iactual = 0; iactual < cmesh->NElements(); iactual++){
      if (cmesh->ElementVec()[iactual]) nactualel++;
    }  
  
  if (0) {
    ofstream cmeshout("Coarsen_cmesh.txt");
    cmesh->Print(cmeshout);
    ofstream gmeshout("Coarsen_gmesh.txt");
    gmesh->Print(gmeshout);
  }  

// /** Testing Coarsen method */
//   for(int icoarsen = 0; icoarsen < ncoarsen; icoarsen++){
//     const int n = gmesh->NElements();
//     for(int iel = 0; iel < n; iel++){
//       TPZGeoEl * father = gmesh->ElementVec()[iel];
//       if (!father) continue;
//       if (!father->HasSubElement()) continue; //if father has no subelements it is not a father
//       if ( father->SubElement(0)->HasSubElement() ) continue; //I am looking for a father of elements that are not father theirselves (I dont want a grandfather)
//       int nsub = father->NSubElements();
//       TPZManVector<int> subindices(nsub);
//       bool CreateDiscontinuous = false;
//       for(int isub = 0; isub < nsub; isub++){
//         TPZGeoEl * subel = father->SubElement(isub);
//         if (!subel) continue;
//         TPZCompEl * cel = subel->Reference();
//         if (!cel) continue;
//         int index = cel->Index();
//         subindices[isub] = index;        
//         TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc *>(cel);
//         if (disc) CreateDiscontinuous = true;
//       }//for isub
//       int coarseindex = -1;      
//       cmesh->Coarsen(subindices, coarseindex, CreateDiscontinuous);
//     }//for iel  
//   }//for icoarsen


  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  return cmesh;
}//mesh 1 elemento


TPZCompMesh * mesh(TIPO tipo, const int ndivision, const int ncoarsen){

  TPZCompEl::gOrder = 2;

  const int nnodes = 6;
  const int nelem  = 2;
  TPZVec < TPZVec< REAL > > co(nnodes);
  for(int i = 0; i < co.NElements(); i++) co[i].Resize(2);
  
  co[0][0] = 0.; co[0][1] = 0.;
  co[1][0] = 1.; co[1][1] = 0.;
  co[2][0] = 1.; co[2][1] = 1.;
  co[3][0] = 0.; co[3][1] = 1.;
  co[4][0] = 2.; co[4][1] = 0.;
  co[5][0] = 2.; co[5][1] = 1.;
  
  TPZVec< TPZVec<int> > indices(nelem);
  for(int i = 0; i < indices.NElements(); i++) indices[i].Resize(4);
  
  indices[0][0] = 0; indices[0][1] = 1; indices[0][2] = 2; indices[0][3] = 3;
  indices[1][0] = 1; indices[1][1] = 4; indices[1][2] = 5; indices[1][3] = 2;
  
  TPZVec< TPZGeoEl* > elvec(nelem);
  
  TPZGeoMesh *gmesh = new TPZGeoMesh();

  int nod;
  for(nod = 0; nod < nnodes; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el = 0; el < nelem; el++) {
    TPZManVector<int,4> nodind(4);
    for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
  }

   gmesh->BuildConnectivity();
   TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh);
   TPZGeoElBC gbc2(elvec[1],4,-1,*gmesh);
//   TPZVec<int> nodind;
//   int index = -1;
//   nodind.Resize(2);
//   nodind[0] = 0; nodind[1] = 1;
//   gmesh->CreateGeoElement(EOned,nodind,-1,index);
//   nodind[0] = 1; nodind[1] = 4;
//   gmesh->CreateGeoElement(EOned,nodind,-1,index);
//   gmesh->BuildConnectivity();

//FAZENDO CONTINUO DESCONTINUO
  TPZVec<TPZGeoEl*> continuous, discontinuous;
  if (tipo == Continuous){
    continuous.Resize(2);
    discontinuous.Resize(0);
    continuous[0] = elvec[0];
    continuous[1] = elvec[1];
  }//if Continuous
  
  if (tipo == Discontinuous){
    continuous.Resize(0);
    discontinuous.Resize(2);
    discontinuous[0] = elvec[0];
    discontinuous[1] = elvec[1];
  }//if Continuous

  if (tipo == Both){
    continuous.Resize(1);
    discontinuous.Resize(1);
    continuous[0] = elvec[0];
    discontinuous[0] = elvec[1];
  }//if Continuous  
  
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  SetMaterial(*cmesh);
  
  cmesh->AutoBuildContDisc(continuous, discontinuous);
// cmesh->SetAllCreateFunctionsDiscontinuous();
// cmesh->AutoBuild();
//   TPZInterfaceElement::SetCalcStiffContDisc();
  
  if (0) {
    ofstream cmeshout("Original_cmesh.txt");
    cmesh->Print(cmeshout);
    ofstream gmeshout("Original_gmesh.txt");
    gmesh->Print(gmeshout);
  }
  
  /** Testing divide method */  
{//Continuous elements  
  set<int> newindices, oldindices;
  for(int i = 0; i < continuous.NElements(); i++) oldindices.insert(continuous[i]->Index());
  
  for(int j = 0; j < ndivision; j++){
    newindices.clear();
    set<int>::iterator e, w;
    e = oldindices.end();
    for(w = oldindices.begin(); w != e; w++){
      const int index = *w;
      TPZManVector<int,8> subindex; 
      cmesh->Divide(index, subindex, 1);
      for(int subi = 0; subi < subindex.NElements(); subi++){
        newindices.insert( subindex[subi] );
      }//for subi
    }//for w
    oldindices.clear();
    e = newindices.end();
    for(w = newindices.begin(); w != e; w++) oldindices.insert(*w);
    cmesh->ExpandSolution();
//     cmesh->AdjustInterfaceElements();
//     cmesh->AdjustBoundaryElements();
//     cmesh->CleanUpUnconnectedNodes();
  }//for j
}//continuous

{//Discontinuous elements  
  set<int> newindices, oldindices;
  for(int i = 0; i < discontinuous.NElements(); i++) oldindices.insert(discontinuous[i]->Index());
  
  for(int j = 0; j < ndivision; j++){
    newindices.clear();
    set<int>::iterator e, w;
    e = oldindices.end();
    for(w = oldindices.begin(); w != e; w++){
      const int index = *w;
      TPZManVector<int,8> subindex; 
      cmesh->Divide(index, subindex, 1);
      for(int subi = 0; subi < subindex.NElements(); subi++){
        newindices.insert( subindex[subi] );
      }//for subi
    }//for w
    oldindices.clear();
    e = newindices.end();
    for(w = newindices.begin(); w != e; w++) oldindices.insert(*w);
    cmesh->ExpandSolution();
    cmesh->RemakeAllInterfaceElements();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();    
  }//for j
}//discontinuous
  
  if (0) {
    ofstream cmeshout("Refined_cmesh.txt");
    cmesh->Print(cmeshout);
    ofstream gmeshout("Refined_gmesh.txt");
    gmesh->Print(gmeshout);
  }

/** Testing Coarsen method */
  for(int icoarsen = 0; icoarsen < ncoarsen; icoarsen++){
    const int n = cmesh->NElements();
    
    /** counting number of non NULL elements */
    int nactualel = 0;
    for(int iactual = 0; iactual < n; iactual++){
      if (cmesh->ElementVec()[iactual]){
        if (cmesh->ElementVec()[iactual]->Reference()->Dimension() == 2) nactualel++;
      }
    }
    
    TPZVec< TPZCompEl * > CopyElVec( cmesh->ElementVec().NElements() );
    for(int iel = 0; iel < n; iel++){
      TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement *>(cmesh->ElementVec()[iel]);
      if (face) CopyElVec[iel] = NULL;
      else CopyElVec[iel] = cmesh->ElementVec()[iel];
    }

        
    for(int iel = 0; iel < n; iel++){
      TPZCompEl * cel = CopyElVec[iel];
      if (!cel) continue;      
      MElementType tipoel = cel->Type();
      int index = cel->Index();
      TPZGeoEl * father = cel->Reference()->Father();
      if (!father) continue;
      int ID = father->Id();
      int nsub = father->NSubElements();
      TPZManVector<int,4> subindices(nsub);
      bool CreateDiscontinuous = false;
      for(int isub = 0; isub < nsub; isub++){
        TPZGeoEl * subel = father->SubElement(isub);
        if (!subel) continue;
        TPZCompEl * cel = subel->Reference();
        if (!cel) continue;
        int index = cel->Index();
        subindices[isub] = index;
        TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc *>(cel);
        if (disc) CreateDiscontinuous = true;
             
      }//for isub
      int coarseindex = -1;      
      
      for(int isub = 0; isub < subindices.NElements(); isub++) CopyElVec[ subindices[isub] ] = NULL;
      cmesh->Coarsen(subindices, coarseindex, CreateDiscontinuous);
    }//for iel  

    cmesh->RemakeAllInterfaceElements();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

  }//for icoarsen
  
    /** counting number of non NULL elements */
    int nactualel = 0;
    for(int iactual = 0; iactual < cmesh->NElements(); iactual++){
      if (cmesh->ElementVec()[iactual]) nactualel++;
    }  
  
  if (0) {
    ofstream cmeshout("Coarsen_cmesh.txt");
    cmesh->Print(cmeshout);
    ofstream gmeshout("Coarsen_gmesh.txt");
    gmesh->Print(gmeshout);
  }  

// /** Testing Coarsen method */
//   for(int icoarsen = 0; icoarsen < ncoarsen; icoarsen++){
//     const int n = gmesh->NElements();
//     for(int iel = 0; iel < n; iel++){
//       TPZGeoEl * father = gmesh->ElementVec()[iel];
//       if (!father) continue;
//       if (!father->HasSubElement()) continue; //if father has no subelements it is not a father
//       if ( father->SubElement(0)->HasSubElement() ) continue; //I am looking for a father of elements that are not father theirselves (I dont want a grandfather)
//       int nsub = father->NSubElements();
//       TPZManVector<int> subindices(nsub);
//       bool CreateDiscontinuous = false;
//       for(int isub = 0; isub < nsub; isub++){
//         TPZGeoEl * subel = father->SubElement(isub);
//         if (!subel) continue;
//         TPZCompEl * cel = subel->Reference();
//         if (!cel) continue;
//         int index = cel->Index();
//         subindices[isub] = index;        
//         TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc *>(cel);
//         if (disc) CreateDiscontinuous = true;
//       }//for isub
//       int coarseindex = -1;      
//       cmesh->Coarsen(subindices, coarseindex, CreateDiscontinuous);
//     }//for iel  
//   }//for icoarsen


  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  return cmesh;
}//mesh

using namespace std;

int main(){

int nmaxeq = 5000;//numero maximo de equacoes do problema

//  TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);


  TPZCompMesh * cmesh = mesh(Both/*Discontinuous*/, 3, 1); 
    
  TPZGeoMesh *gmesh = cmesh->Reference();

  TPZAnalysis an(cmesh);

#define direct
#ifdef direct
  //TPZParFrontStructMatrix <TPZFrontNonSym> full(cmesh);
  TPZFStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect( ELU /*ECholesky*/ /*ELDLt*/ );
  an.SetSolver(step);
#endif

//#define iter
#ifdef iter
  cout << "ITER_SOLVER" << endl;  
  /*TPZFStructMatrix */TPZSpStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);  
  TPZStepSolver step( full.Create() );
  an.SetSolver(step);  
  REAL tol = 1.e-14;
// Sem pre-condicionador 
  TPZCopySolve precond( full.Create() );
  step.ShareMatrix( precond );  
  step.SetGMRES( 2000, 20, precond, tol, 0 ); 
  an.SetSolver(step);
#endif

   std::cout << "\nNumero de equacoes: " << an.Mesh()->NEquations() << std::endl;
   if (an.Mesh()->NEquations() > nmaxeq) { 
     cout << "skipping simulation...\n" << endl;
     delete cmesh;
     delete gmesh;
     return -1; //break;
   }
  
  char filedx[20];
  sprintf(filedx,"sol.dx"/*,p,h*/);  

  if (0) {
    ofstream cmeshout("cmesh.txt");
    cmesh->Print(cmeshout);
    ofstream gmeshout("gmesh.txt");
    gmesh->Print(gmeshout);
  }
  
  an.Run();

/**** Aqui faz DX ****/
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(1);
  scalnames[0] = "Solution";
  vecnames[0] = "Derivate";
  an.DefineGraphMesh(2,scalnames,vecnames,filedx);
  an.PostProcess(2);
  
  delete cmesh;
  delete gmesh;

  return 0;
}

