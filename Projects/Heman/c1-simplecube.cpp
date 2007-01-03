#include "pzvec.h"

#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzgnode.h"

#include "pzmatrix.h"

#include "pzpoisson3d.h"
#include "pzmat2dlin.h"
#include "pzbndcond.h"

using namespace std;
static TPZCompMesh * CreateCubeMesh();

//*************************************
//************Option 0*****************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateCubeMesh(){
  //malha 2 cubos
  const int nelem = 2;
  //número de nós
  const int ncoord = 12;
  //TPZVec<REAL> coord(ncoord,0.);
  REAL Coord[ncoord][3] = { { 0.,0.,0.} , 
		       { 1.,0.,0.} ,
		       { 2.,0.,0.} ,
		       { 0.,1.,0.} , 
		       { 1.,1.,0.} ,
		       { 2.,1.,0.} ,
                       { 0.,0.,1.} , 
		       { 1.,0.,1.} ,
		       { 2.,0.,1.} ,
		       { 0.,1.,1.} , 
		       { 1.,1.,1.} ,
		       { 2.,1.,1.}                
                                            
                                   };
  int Connect[nelem][8] = { {0,1,4,3,6,7,10,9},
			    {1,2,5,4,7,8,11,10} };
  int nConnect[nelem] = {8,8};
  
  // criar um objeto tipo malha geometrica
  TPZGeoMesh *geomesh = new TPZGeoMesh();
  
  // criar nos
  int i,j;
  for(i=0; i<(ncoord); i++) {
    int nodind = geomesh->NodeVec().AllocateNewElement();
    TPZVec<REAL> coord(3);
    for (j=0; j<3; j++) {
      coord[j] = Coord[i][j];
    }
    geomesh->NodeVec()[nodind] = TPZGeoNode(i,coord,*geomesh);
  }

  // criação dos elementos
  TPZGeoEl *gel[nelem];

  for(i=0;i<nelem;i++) {  
    TPZVec<int> indices(nConnect[i]);
    for(j=0;j<nConnect[i];j++) {
      indices[j] = Connect[i][j];
    }
    int index;
    switch (nConnect[i]){
    case (4): 
      gel[i] = geomesh->CreateGeoElement(EQuadrilateral,indices,1,index);
      break;
    case(3):
      gel[i] = geomesh->CreateGeoElement(ETriangle,indices,1,index);
      break;
    case (8) :
      gel[i] = geomesh->CreateGeoElement(ECube,indices,1,index);  
      break; 
    default:
      cout << "Erro : elemento não implementado" << endl;
    }
  }
 
  // Descomentar o trecho abaixo para habilitar a
  // divisão dos elementos geométricos criados 
                      
  geomesh->BuildConnectivity();
  //  geomesh->Print(cout);

  //Divisão dos elementos
  // TPZVec<TPZGeoEl *> sub,subsub;
  //  gel[0]->Divide(sub);
  //  sub[0]->Divide(subsub);
  //  subsub[2]->Divide(sub);
  
  //  for (i=0;i< (sub.NElements()-1) ;i++){
  //    sub[i]->Divide(subsub);
  //  }
  
  // Criação das condições de contorno geométricas
  TPZGeoElBC heman_1(gel[0],20,-1,*geomesh);
  TPZGeoElBC heman_2(gel[1],20,-1,*geomesh); 
  //  geomesh->BuildConnectivity2();
  //geomesh->Print(cout);

  // Criação da malha computacional
  TPZCompMesh *comp = new TPZCompMesh(geomesh);

  // Criar e inserir os materiais na malha
  TPZAutoPointer<TPZMaterial> mat = new TPZMatPoisson3d(1,3);
  comp->InsertMaterialObject(mat);
 
  TPZAutoPointer<TPZMaterial> meumat = mat;

  // Condições de contorno
  // Dirichlet
  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  TPZAutoPointer<TPZMaterial> bnd = meumat->CreateBC (meumat,-1,0,val1,val2);
  comp->InsertMaterialObject(bnd);
  bnd = meumat->CreateBC (meumat,-1,0,val1,val2);
  
  // comp->Print(cout);

  // Ajuste da estrutura de dados computacional
  comp->AutoBuild();
  return comp;
}
