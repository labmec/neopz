#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include "pzvec.h"
#include "pzerror.h"
#include <pzmat1dlin.h>
#include <pzmattest3d.h>
#include <pzmattest.h>
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
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
#include "pzfmatrix.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"

//solvers
#include "pzsolve.h"
#include "pzstepsolver.h"

//IO
#include <iostream>
#include <fstream> 
using namespace std;
int gDebug = 0;

void error(char *msg)
{
  cout << msg << endl;
}

int ReadNodes(istream &arq,TPZGeoMesh *gmesh){
   int nnodes = 0;
   TPZVec<REAL> coord (3,0.);
   arq >> nnodes;
   if(nnodes <= 0){
      cout << "TMBGeometricMesh::ReadNodes. Error: invalid number of nodes.\n";
      return 0;
   }
   int i,j;
   gmesh->NodeVec().Resize(nnodes);
   for (i = 0; i < nnodes; i++){
      for (j = 0; j < 3; j++){
         arq >> coord[j];
      }
      //Node Id = node index begin in 0!
      gmesh->NodeVec()[i].Initialize(i, coord, *gmesh);
   }
   cout << "Nodes sucessfully read: ";
   cout << nnodes;
   cout << " nodes\n";
   return 1;
}


int ReadElements(istream &arq,TPZGeoMesh *gmesh){
   int nelements = 0;
   TPZVec<int> connectids;
   int nelnodes = 0;
   arq >> nelements;
   if(nelements <= 0){
      cout << "TMBGeometricMesh::ReadElements. Error: invalid number of elements.\n";
      return 0;
   }
   int i = 0, j;
   int npoints = 0;
   int nlinears = 0;
   int ntriangs = 0;
   int ntetras = 0;
   int nprisms = 0;
   int npyrams = 0;
   int ncubes = 0;
   int index;
   for (i = 0; i < nelements; i++){
      arq >> nelnodes;
      if (nelnodes < 1){
         cout << "ReadElements. Error: invalid number of element nodes.\n";
         return 0;
      }
      connectids.Resize(nelnodes);
      int aux;
      for (j = 0; j < nelnodes; j++){
         arq >> aux;
         // our node index and node ids begin from 0 - Azevedo from 1
         connectids[j] = aux - 1;
         if (aux < 1){
            cout << "ReadElements. Error: invalid node id in connectivity.\n";
            return 0;
         }
      }
      switch (nelnodes){
         case 4:{
            gmesh->CreateGeoElement( ETetraedro, connectids, 1, index );
            ntetras++;
            break;
         }
         case 5:{
            gmesh->CreateGeoElement( EPiramide, connectids, 1, index );
            npyrams++;
            break;
         }
         case 6:{
            gmesh->CreateGeoElement( EPrisma, connectids, 1, index );
            nprisms++;
            break;
         }
         case 8:{
            gmesh->CreateGeoElement( ECube, connectids, 1, index );
            ncubes++;
            break;
         }
         default:{
            cout << "TMBGeometricMesh::ReadElements. Error: element ";
            cout << i+1;
            cout << " with ";
            cout << nelnodes;
            cout << " nodes. Unknwon element type!\n";
            return 0;
         }
      }
   }
   cout << "Elements sucessfully read: ";
   cout << nelements << " elements, being\n";
   cout << "\t"<< npoints  << " points\n";
   cout << "\t"<< nlinears << " linears\n";
   cout << "\t"<< ntriangs << " triangles\n";
   cout << "\t"<< ntetras  << " tetrahedra\n";
   cout << "\t"<< npyrams  << " pyramids\n";
   cout << "\t"<< nprisms  << " prims\n";
   cout << "\t"<< ncubes   << " cubes\n";
   return 1;
}


int ReadBoundaries(istream &arq,TPZGeoMesh *gmesh){
   int nbound = 0;
   TPZVec<int> face_connects (4,-1);
   arq >> nbound;
   if((nbound < 1)){
      cout << "ReadBoundaries. Error: invalid number of boundary conditions.\n";
      return 0;
   }
   int i, j;
   int type, index;
   for (i = 0; i < nbound; i++){
      int nodeid =  -1;
      arq >> type;
      if((type < 1) || (type > 8)){
         cout << "TMBGeometricMesh::ReadBoundaries. Error: invalid boundary condition.\n";
         return 0;
      }
      for (j = 0; j < 4; j++){
         arq >> nodeid;
         face_connects[j] = nodeid - 1;
      }
      if (face_connects[3] == -1){
         TPZVec<int>  aux(3,0);
         for (j = 0; j < 3; j++){
            aux[j] = face_connects[j];
         }

         //Note: type+1
         gmesh->CreateGeoElement( ETriangle, aux, type+1, index );
      }
      else{
         //Note: type+1
         gmesh->CreateGeoElement( EQuadrilateral, face_connects, type+1, index );
      }
   }
   cout << "Boundary conditions sucessfully read: ";
   cout << nbound;
   cout << " boundaries\n";
   return 1;
}

int ReadFathers(istream &arq,TPZGeoMesh *gmesh){
  int nelements = 0;
  TPZVec<int> connectids;
  arq >> nelements;
  TPZVec<int> RegSubElement (nelements,0);
  if(nelements <= 0) {
    cout << "TMBGeometricMesh::ReadFathers. Error: invalid number of elements.\n";
    return 0;
  }
  int i = 0;
  int index;
  for (i = 0; i < nelements; i++){
    arq >> index;
    if (index < 0){
      //The element hasn't father
      continue;
    }
    TPZGeoEl *gel = gmesh->ElementVec()[i];
    if (!gel) {
      cout << "TMBGeometricMesh::ReadFathers. Error: invalid element.\n";
      return 0;
    }
    TPZGeoEl *father = gmesh->ElementVec()[index];
    if (!father) {
      cout << "TMBGeometricMesh::ReadFathers. Error: invalid father.\n";
      return 0;
    }
    //The subelements are read and write in asceding sequence!!!
    gel->SetFather(father);
    // THIS WILL WORK ONLY FOR QUAD AND TRIANGLE BEWARE !!!!
    int ncor = gel->NCornerNodes();
    int ic;
    for(ic=0 ;ic<ncor; ic++)
    {
      TPZGeoElSide small(gel,ic),large(father,ic);
      if(small.NeighbourExists(large)) 
      {
        father->SetSubElement(ic,gel);
        break;
      }
      if(ic == ncor-1)
      {
        cout << "Father Son non implemented element type?\n";
      }
    }
//    father->SetSubElement(RegSubElement[index],gel);
    RegSubElement[index] ++;      
  }
  return 1;  
}


int main(){
  TPZGeoMesh * gmesh = new TPZGeoMesh;
  ifstream nodes ("adapt_nodes.2");
  ifstream incid ("adapt_incid.3");
  ifstream father ("adapt.0");
  ifstream bcs ("funil_ref9.5");
  
  ReadNodes(nodes,gmesh);
  //ReadBoundaries (bcs,gmesh);  
  ReadElements(incid,gmesh);
  ReadBoundaries (bcs,gmesh);
  
  cout << endl << endl << "Antes Build Connectivities" << endl << endl;
  
  gmesh->BuildConnectivity();
  
  cout << endl << endl << "Depois Build Connectivities" << endl << endl;
  
  ReadFathers(father,gmesh);
  

 // gmesh->Print(cout);
  
//  gmesh->BuildConnectivity();
  
  ofstream out("output.txt");
  gmesh->Print(out);
  
  
//     TPZGeoElement< TPZShapeTetra, TPZGeoTetrahedra, TPZRefTetrahedra >::
//       SetCreateFunction( TPZCompElDisc::CreateDisc );
// 
//    TPZGeoElement< TPZShapeTriang, TPZGeoTriangle, TPZRefTriangle >::
//       SetCreateFunction( TPZCompElDisc::CreateDisc );

   TPZGeoElement< TPZShapeCube, TPZGeoCube, TPZRefCube >::
      SetCreateFunction( TPZCompElDisc::CreateDisc );

   TPZGeoElement< TPZShapeQuad, TPZGeoQuad, TPZRefQuad >::
      SetCreateFunction( TPZCompElDisc::CreateDisc );

//    TPZGeoElement< TPZShapePrism, TPZGeoPrism, TPZRefPrism >::
//       SetCreateFunction( TPZCompElDisc::CreateDisc );
// 
//    TPZGeoElement< TPZShapePiram, TPZGeoPyramid, TPZRefPyramid >::
//       SetCreateFunction( TPZCompElDisc::CreateDisc );

   TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
   TPZMatPoisson3d *mat = new TPZMatPoisson3d (1,3);
   cmesh->InsertMaterialObject(mat);
   
   TPZMatPoisson3d *mat_bc[10];
   for (int i=2;i<=9;i++){
     mat_bc[i] = new TPZMatPoisson3d (i-1,3);
     cmesh->InsertMaterialObject(mat_bc[i]);
   }

  cmesh->AutoBuild();
  cmesh->Print(cout);
  
  delete cmesh;
  delete gmesh;
  cout << "fim" << endl;
  
  return 0;
}

