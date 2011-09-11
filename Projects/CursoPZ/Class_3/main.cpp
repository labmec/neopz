/**
 * @file
 * @brief Implements the use of the jacobian method as tutorial example of the matrix NeoPZ module
 */
#include "pzvec.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include <pzskylmat.h>
#include <pzstepsolver.h>
#include <pzgmesh.h>
#include <pzgeoel.h>
//IO
#include <iostream>
using namespace std;


int main(){
  REAL coordinates [6][3] = {	{0.,0.,0.},{1.,0.,1.},{2.,0.,0.},
                                {0.,1.,0.},{1.,1.,1.},{2.,1.,0.}};
                   //0 , 1  , 2 , 3 ,  4 , 5
 int nodeids [6] = {100,1100,200,300,1200,400};
 int elconnect[2][4] = {{0,1,4,3},{4,1,2,5}};
 int elid [2] = {10,20};
 int i,j;
 TPZGeoMesh *gmesh = new TPZGeoMesh;
 int index = 0;
 TPZVec<REAL> coord(3,0.);
 //Creates the nodes like specified
 for (i=0;i<6;i++){
   for (j=0;j<3;j++) coord[j] = coordinates[i][j];
   index = gmesh->NodeVec().AllocateNewElement();
   gmesh->NodeVec()[index] = TPZGeoNode(nodeids[i],coord,*gmesh);
 }

 TPZGeoEl *elvec[2];
 TPZVec<int> connect(4,0);
 //Creates the elements like specified
 for (i=0;i<2;i++){
   for (j=0;j<4;j++){
     connect[j] = elconnect[i][j];
   }
   elvec[i] = gmesh->CreateGeoElement(EQuadrilateral,connect,1,elid[i]);
 }

 gmesh->BuildConnectivity();
 gmesh->Print();

 //Given point-->> (1.50,0.65,0.50)
 int dim = elvec[1]->Dimension();
 TPZVec<REAL> gpt(3,0.);
 TPZVec<REAL> mpt(dim,0.);
 gpt[0] = 1.5; gpt[1] = 0.65; gpt[2] = 0.5;
 elvec[1]->ComputeXInverse (gpt, mpt);
 
 TPZFMatrix jac(dim,dim,0.);
 TPZFMatrix axes(3,3,0.);
 REAL detjac;
 TPZFMatrix jacinv(dim,dim,0.);
 elvec[1]->Jacobian(mpt, jac, axes, detjac, jacinv);
 cout << "Jacobiana\n" << jac << endl << "Axes\n" << axes
       << "Det = " << detjac << "Inversa Jacobiana\n" << jacinv << endl;			

 return 0;
}

