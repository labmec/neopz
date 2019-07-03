#include "pzvec.h"

#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzgnode.h"

#include "pzmatrix.h"

#include "pzelasmat.h"
//#include "pzplaca.h"
#include "pzmat2dlin.h"


static TPZCompMesh *CreateTriangularMesh();


//*************************************
//************Option 2*****************
//******Simple Triangular Mesh*********
//*************************************
TPZCompMesh *CreateTriangularMesh(){

  const int nelem = 1;
  const int nnode = 3;
  
  REAL co[nnode][3] = {{0.,0.,0.},{1.,0.,0.},{0.,1.,0.}};
  int indices[nelem][nnode] = {{0,1,2}};
  
  TPZGeoEl *elvec[nelem];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZVec<REAL> coord(3);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    coord[2] = co[nod][2];
    gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
  }

  int el;
  for(el=0; el<nelem; el++) {
    TPZVec<int> nodind(3);
    for(nod=0; nod<3; nod++) nodind[nod]=indices[el][nod];
    //    elvec[el] = new TPZGeoElT2d(el,nodind,1,*gmesh);
    int index;
    elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
  }
  
  gmesh->BuildConnectivity2();
  //  TPZStack<TPZGeoEl*> subel;
  //  elvec[0]->Divide(subel);
  
  //  TPZGeoElBC gbc;
  
  // bc -1 -> Dirichlet
  TPZGeoElBC gbc1(elvec[0],3,-1,*gmesh);

  // bc -2 -> Neumann at the right x==1
  TPZGeoElBC gbc2(elvec[0],5,-2,*gmesh);
  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  
 // inserir os materiais
 TPZMaterial *mat = new TPZElasticityMaterial(1,1.e5,0.2,0,0);
 cmesh->InsertMaterialObject(mat);
 
 TPZMaterial *meumat = mat;

 // inserir a condicao de contorno
 TPZFMatrix val1(3,3,0.),val2(3,1,0.);

 TPZMaterial *bnd = meumat->CreateBC (-1,0,val1,val2);
 cmesh->InsertMaterialObject(bnd);
 bnd = meumat->CreateBC (-2,0,val1,val2);

 val2(0,0)=1.;
 bnd = meumat->CreateBC (-2,1,val1,val2);
 cmesh->InsertMaterialObject(bnd);

 cmesh->AutoBuild();
 cmesh->AdjustBoundaryElements();
 cmesh->CleanUpUnconnectedNodes();
  
 return cmesh;
}
