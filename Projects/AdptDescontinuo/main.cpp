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
   cmesh->InsertMaterialObject(mat);
   
   cmesh->AutoBuild();
   
   TPZVec<int> subindex;
   
   int nelem = cmesh->NElements();
   cmesh->Print(cout);
   cmesh->Reference()->Print(cout);
   
   for (i=0;i<nelem;i++){
     TPZCompEl *cel = cmesh->ElementVec()[i];
     if (!cel || cel->Type() != EDiscontinuous || cel->Reference()->Dimension() != 3) continue;
     cel->Divide(i,subindex,0);
   }
   
   gmesh = cmesh->Reference();
   gmesh->ResetReference();
   delete cmesh;
   gmesh->Print(cout);
}
