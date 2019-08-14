
#include "pzvec.h"

#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzgnode.h"

#include "pzmatrix.h"

#include "pzelasmat.h"
#include "pzbndcond.h"
//#include "pzplaca.h"
#include "pzmat2dlin.h"

using namespace std;

static TPZCompMesh *CreatePlanMesh();

//*************************************
//************Option 3*****************
//*****Plane Quad & Triang Mesh********
//*************************************
TPZCompMesh *CreatePlanMesh() {
  
  REAL co[5][2] = {
    {0.,0.},
    {1.,0.},
    {2.,0.},
    {0.,1.},
    {1.,1.}
  };

  int indices[2][4] = {
    {0,1,4,3},
    {1,2,4,-1}
  };
  
  TPZGeoEl *elvec[2];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  
  int nnode = 5;
  int nelem = 2;

  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZVec<REAL> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    TPZGeoNode pznode(nod,coord,*gmesh);
    gmesh->NodeVec()[nodind] = pznode;
  }

  int el;

  for(el=0; el<nelem; el++) {
    int ncnodes = el > 0 ? 3 : 4;
    TPZVec<int> nodind(ncnodes);
    for(nod=0; nod<ncnodes; nod++) nodind[nod]=indices[el][nod];
    int index;
    if (el == 0){
      //      elvec[el] = new TPZGeoElQ2d(el,nodind,1,*gmesh);
      elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }else{
      //      elvec[el] = new TPZGeoElT2d(el,nodind,1,*gmesh);
      elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
    }
  }

  gmesh->BuildConnectivity2();
  
  // bc -1 -> Dirichlet
  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh);
  TPZGeoElBC gbc2(elvec[1],3,-2,*gmesh);

  // bc -4 -> Neumann at the top y==1
  TPZGeoElBC gbc3(elvec[0],6,-3,*gmesh);

  gmesh->SetName ("Original Geometric Mesh");
  gmesh->Print(cout);

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetName ("Original Computational Mesh");

  // Criar e inserir os materiais na malha
  TPZAutoPointer<TPZMaterial> mat = new TPZElasticityMaterial(1,1.e5,0.2,0,0);
  cmesh->InsertMaterialObject(mat);
 
  TPZAutoPointer<TPZMaterial> meumat = mat;

  // Condições de contorno
  // Dirichlet
  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  TPZAutoPointer<TPZMaterial> bnd = meumat->CreateBC (meumat,-1,0,val1,val2);
  cmesh->InsertMaterialObject(bnd);

  bnd = meumat->CreateBC (meumat,-2,0,val1,val2);
  cmesh->InsertMaterialObject(bnd);

  // Neumann
  val2(0,0) = 1.;
  bnd = meumat->CreateBC (meumat,-3,1,val1,val2);
  cmesh->InsertMaterialObject(bnd);

  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  //cmesh->ExpandSolution();

  //cmesh->Print(cout);
  return cmesh;
}
