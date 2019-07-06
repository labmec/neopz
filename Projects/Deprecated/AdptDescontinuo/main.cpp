#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include "pzvec.h"
//#include "pzerror.h"
//#include <pzmat1dlin.h>
//#include <pzmattest3d.h>
//#include <pzmattest.h>
//#include "TPZRefPattern.h"
//#include "tpzgeoelrefpattern.h"
#include <pzpoisson3d.h>

//Shape Classes
#include <pzshapepoint.h>
#include <pzshapelinear.h>
#include <pzshapetriang.h>
#include <pzshapequad.h>
#include <pzshapetetra.h>
#include <pzshapepiram.h>
#include <pzshapeprism.h>
#include <pzshapecube.h>

//Refine Classes

#include <TPZRefCube.h>
#include <TPZRefLinear.h>
#include <pzrefpoint.h>
#include <pzrefprism.h>
#include <pzrefpyram.h>
#include <pzrefquad.h>
#include <pzreftetrahedra.h>
#include <pzreftriangle.h>

//Geometry Classes
#include <TPZGeoElement.h>
#include <pzgeopoint.h>
#include <TPZGeoLinear.h>
#include <pzgeotriangle.h>
#include <pzgeoquad.h>
#include <pzgeotetrahedra.h>
#include <pzgeopyramid.h>
#include <pzgeoprism.h>
#include <TPZGeoCube.h>

//Computational classes
#include <TPZCompElDisc.h>

//matrix and analysis classes
//#include "pzfmatrix.h"
//#include "pzanalysis.h"
//#include "pzskylstrmatrix.h"

//solvers
//#include "pzsolve.h"
//#include "pzstepsolver.h"



//IO
#include <iostream>
using namespace std;
int gDebug = 0;

void WriteElement (TPZGeoEl *el,int elindex, ofstream &arq,TPZVec<int> &elementtype);
void WriteMesh(TPZGeoMesh *mesh,ofstream &arq);

int main(){
  
  REAL coordinates [8][3]={{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.},
                           {0.,0.,1.},{1.,0.,1.},{1.,1.,1.},{0.,1.,1.}};
  TPZVec <REAL> coord (3,0.);
  int i,j,index;
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  
  for (i=0;i<8;i++){
    for (j=0;j<3;j++){
      coord[j] = coordinates[i][j];
    }
    index = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[index] = TPZGeoNode(index,coord,*gmesh);
  }

  TPZVec <int> elconnect (8,0);
  for (i=0;i<8;i++) elconnect[i] = i;

  gmesh->CreateGeoElement(ECube,elconnect,1,index);
  
 
  gmesh->BuildConnectivity();
   
  TPZVec <TPZGeoEl *> subelindex (8,0);
  TPZGeoEl *gel = gmesh->ElementVec()[0];
  gel->Divide(subelindex);
  
  elconnect.Resize(4);
  elconnect [0] = 0;
  elconnect [1] = 4;
  elconnect [2] = 7;
  elconnect [3] = 3;
  
  gmesh->CreateGeoElement(EQuadrilateral,elconnect,2,index);
  
  gmesh->BuildConnectivity();
  
  TPZGeoElement< TPZShapeTetra, TPZGeoTetrahedra, TPZRefTetrahedra >::
      SetCreateFunction( TPZCompElDisc::CreateDisc );

   TPZGeoElement< TPZShapeTriang, TPZGeoTriangle, TPZRefTriangle >::
      SetCreateFunction( TPZCompElDisc::CreateDisc );

   TPZGeoElement< TPZShapeCube, TPZGeoCube, TPZRefCube >::
      SetCreateFunction( TPZCompElDisc::CreateDisc );

   TPZGeoElement< TPZShapeQuad, TPZGeoQuad, TPZRefQuad >::
      SetCreateFunction( TPZCompElDisc::CreateDisc );

   TPZGeoElement< TPZShapePrism, TPZGeoPrism, TPZRefPrism >::
      SetCreateFunction( TPZCompElDisc::CreateDisc );

   TPZGeoElement< TPZShapePiram, TPZGeoPyramid, TPZRefPyramid >::
      SetCreateFunction( TPZCompElDisc::CreateDisc );

   TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
   TPZMatPoisson3d *mat = new TPZMatPoisson3d (1,3);
   TPZMatPoisson3d *mat1 = new TPZMatPoisson3d (2,3);
   cmesh->InsertMaterialObject(mat);
   cmesh->InsertMaterialObject(mat1);
   
   cmesh->AutoBuild();
   
   TPZVec<int> subindex;
   
   int nelem = cmesh->NElements();
   //cmesh->Print(cout);
   //cmesh->Reference()->Print(cout);
   
   for (i=0;i<nelem;i++){
     TPZCompEl *cel = cmesh->ElementVec()[i];
     if (!cel || cel->Type() != EDiscontinuous || cel->Reference()->Dimension() != 3) continue;
     cel->Divide(i,subindex,0);
   }
   
   gmesh = cmesh->Reference();
   gmesh->ResetReference();
   
//    for (i=0;i<cmesh->NElements();i++){
//       TPZCompEl *cel = cmesh->ElementVec()[i];
//       if (!cel || cel->Type() != EDiscontinuous || cel->Reference()->Dimension() != 3) continue;
//       delete cel;
//    }
   ofstream malha1 ("malha1.txt");
   cmesh->Reference()->Print(malha1);
   cmesh->Print(cout);
   /*
   for (i=0;i<gmesh->NElements();i++){
     TPZGeoEl *gel = gmesh->ElementVec()[i];
     if (!gel) continue;
     gel->RemoveConnectivities();
   }
   */
   cout << "nel antes " << gmesh->ReallyNEl() << endl;
   delete cmesh;
//    cout << "nel depois " << gmesh->ReallyNEl() << endl;
//    
//    gmesh->Print(cout);
//    
//    ofstream arq ("teste.dx");
//    WriteMesh(gmesh,arq);
// 
//    gmesh->Print(cout);
//    
//    gmesh->BuildConnectivity2();
//    
//    gmesh->Print(cout);
//    
//    TPZCompMesh *cmesh2 = new TPZCompMesh(gmesh);
//    TPZMatPoisson3d *mat2 = new TPZMatPoisson3d (1,3);
//    cmesh2->InsertMaterialObject(mat2);
//    cmesh2->AutoBuild();
//    ofstream malha2 ("malha2.txt");
//    cmesh2->Print(cout);
//    cmesh2->Reference()->Print(malha2);
//    
//    delete cmesh2;   
//    gmesh->ResetReference();
//    ofstream arq2 ("teste2.dx");
//    WriteMesh (gmesh,arq2);
   delete gmesh;
}    

void WriteMesh(TPZGeoMesh *mesh,ofstream &arq){
  arq << "object 1 class array type float rank 1 shape 3 items ";
  arq << mesh->NodeVec().NElements() << " data follows" << endl;
  int i;
  //Print Nodes
  for (i=0;i<mesh->NodeVec().NElements(); i++){
    TPZGeoNode *node = &mesh->NodeVec()[i];
    arq /*<< node->Id() << "\t"*/
      << node->Coord(0) << "\t"
      << node->Coord(1) << "\t"
      << node->Coord(2) << endl;
  }
  arq << "object 2 class array type integer rank 1 shape 8 items ";
  arq << mesh->ElementVec().NElements() << " data follows" << endl;
  TPZVec<int> elementtype(mesh->ElementVec().NElements(),0);
  for (i=0;i<mesh->ElementVec().NElements();i++){
    TPZGeoEl *el = mesh->ElementVec()[i];
    if ( !el ) continue;
    WriteElement (el,i,arq,elementtype);
  }
  arq << "attribute \"element type\" string \"cubes\"" << endl
    << "attribute \"ref\" string \"positions\"" << endl;
  arq << "object 3 class array type integer rank 0 items ";
  arq << mesh->ElementVec().NElements() << " data follows" << endl;
  for (i=0;i<mesh->ElementVec().NElements();i++){
    TPZGeoEl *el = mesh->ElementVec()[i];
    if ( !el ) continue;
    arq << elementtype[i] << endl;
  }
  arq << "attribute \"dep\" string \"connections\"" << endl;
  arq << "object 4 class field" << endl
    << "component \"positions\" value 1" << endl
    << "component \"connections\" value 2" << endl
    << "component \"data\" value 3" << endl;
}

void WriteElement (TPZGeoEl *el,int elindex, ofstream &arq,TPZVec<int> &elementtype){
  int ncon = el->NNodes();
  elementtype[elindex] = ncon;
  switch (ncon)  {
  case (2) : {
    //rib
    int ni = el->NodePtr(0)->Id();
    int nf = el->NodePtr(1)->Id();
    arq << ni << "\t" << nf << "\t"
        << ni << "\t" << nf << "\t"
        << ni << "\t" << nf << "\t"
        << ni << "\t" << nf << endl;
    break;
  }
  case (3) : {
    //triangle
    int n0 = el->NodePtr(0)->Id();
    int n1 = el->NodePtr(1)->Id();
    int n2 = el->NodePtr(2)->Id();
    arq << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n2 << "\t"
        << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n2 << endl;
    break;
  }
  case (4) : {
    if (el->Dimension() == 2){
      //quad
      int n0 = el->NodePtr(0)->Id();
      int n1 = el->NodePtr(1)->Id();
      int n2 = el->NodePtr(3)->Id();
      int n3 = el->NodePtr(2)->Id();
      arq << n0 << "\t" << n1 << "\t"
          << n2 << "\t" << n3 << "\t"
          << n0 << "\t" << n1 << "\t"
          << n2 << "\t" << n3 << endl;
    }else{
      //tetrahedre
      int n0 = el->NodePtr(0)->Id();
      int n1 = el->NodePtr(1)->Id();
      int n2 = el->NodePtr(2)->Id();
      int n3 = el->NodePtr(3)->Id();
      arq << n0 << "\t" << n1 << "\t"
          << n2 << "\t" << n2 << "\t"
          << n3 << "\t" << n3 << "\t"
          << n3 << "\t" << n3 << endl;
    }
    break;
  }
  case (5) : {
    //pyramid
    int n0 = el->NodePtr(0)->Id();
    int n1 = el->NodePtr(1)->Id();
    int n2 = el->NodePtr(3)->Id();
    int n3 = el->NodePtr(2)->Id();
    int n4 = el->NodePtr(4)->Id();
    arq << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n3 << "\t"
        << n4 << "\t" << n4 << "\t"
        << n4 << "\t" << n4 << endl;
    break;
  }
  case (6) : {
    //pyramid
    int n0 = el->NodePtr(0)->Id();
    int n1 = el->NodePtr(1)->Id();
    int n2 = el->NodePtr(2)->Id();
    int n3 = el->NodePtr(3)->Id();
    int n4 = el->NodePtr(4)->Id();
    int n5 = el->NodePtr(5)->Id();
    arq << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n2 << "\t"
        << n3 << "\t" << n4 << "\t"
        << n5 << "\t" << n5 << endl;
    break;
  }
  case (8) : {
    int n0 = el->NodePtr(0)->Id();
    int n1 = el->NodePtr(1)->Id();
    int n2 = el->NodePtr(3)->Id();
    int n3 = el->NodePtr(2)->Id();
    int n4 = el->NodePtr(4)->Id();
    int n5 = el->NodePtr(5)->Id();
    int n6 = el->NodePtr(7)->Id();
    int n7 = el->NodePtr(6)->Id();
    arq << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n3 << "\t"
        << n4 << "\t" << n5 << "\t"
        << n6 << "\t" << n7 << endl;
    break;
  }
  default:
    cout << "Erro..." << endl;
  }
  return;
}


