#include "pzvec.h"

#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzgnode.h"

#include "pzmatrix.h"

#include "pzelasmat.h"
#include "pzmat2dlin.h"
#include "pzbndcond.h"

using namespace std;
static TPZCompMesh * CreateSillyMesh();

//*************************************
//************Option 0*****************
//*******L Shape Quadrilateral*********
//*************************************

TPZCompMesh *CreateSillyMesh(){
  //malha 1 quadrado e 2 triangulos
  const int nelem = 1;
  //nmero de n�
  const int ncoord = 4;
  TPZVec<REAL> coord(ncoord,0.);
  REAL Coord[ncoord][3] = { { 0.,0.,0.} , 
		       { 1.,0.,0.} ,
                       { 0.,1.,0.} ,
                       { 0.,0.,1.} };
  int Connect[nelem][4] = { {0,1,2,3} };
  int nConnect[nelem] = {4};
  
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

  // cria�o dos elementos
  TPZGeoEl *gel[nelem];

  for(i=0;i<nelem;i++) {  
    TPZVec<int> indices(nConnect[i]);
    for(j=0;j<nConnect[i];j++) {
      indices[j] = Connect[i][j];
    }
    int index;
    switch (nConnect[i]){
    case (4): 
      gel[i] = geomesh->CreateGeoElement(ETetraedro,indices,1,index,1);
      break;
    case(3):
      gel[i] = geomesh->CreateGeoElement(ETriangle,indices,1,index,1);
      break;
    case(2):
      gel[i] = geomesh->CreateGeoElement(EOned,indices,1,index,1);
      break;
    
      default:
      cout << "Erro : elemento n� implementado" << endl;
    }
  }
 
  // Descomentar o trecho abaixo para habilitar a
  // divis� dos elementos geom�ricos criados 
                      
  geomesh->BuildConnectivity();
  
  //geomesh->Print(cout);

  //Divis� dos elementos
  // TPZVec<TPZGeoEl *> sub,subsub;
  //  gel[0]->Divide(sub);
  //  sub[0]->Divide(subsub);
  //  subsub[2]->Divide(sub);
  
  //  for (i=0;i< (sub.NElements()-1) ;i++){
  //    sub[i]->Divide(subsub);
  //  }
  
  // Cria�o das condi�es de contorno geom�ricas
  TPZGeoElBC heman_1(gel[0],5,-1,*geomesh);
  TPZGeoElBC heman_2(gel[0],9,-1,*geomesh);
  //  geomesh->BuildConnectivity2();
  //geomesh->Print(cout);

  // Cria�o da malha computacional
  TPZCompMesh *comp = new TPZCompMesh(geomesh);

  // Criar e inserir os materiais na malha
  TPZAutoPointer<TPZMaterial> mat = new TPZElasticityMaterial(1,1.e5,0.2,0,0);
  comp->InsertMaterialObject(mat);
 
  TPZAutoPointer<TPZMaterial> meumat = mat;

  // Condi�es de contorno
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
