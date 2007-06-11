#include "includes.h"
#include "pzbndcond.h"

static TPZCompMesh *CreateSimple3DMesh();


//*************************************
//************Option 4*****************
//********Simple 3D Cube Mesh**********
//*************************************
TPZCompMesh *CreateSimple3DMesh() {

  REAL co[8][3] = {
    {0.,0.,0.},
    {1.,0.,0.},
    {1.,1.,0.},
    {0.,1.,0.},
    {0.,0.,1.},
    {1.,0.,1.},
    {1.,1.,1.},
    {0.,1.,1.}
  };
  int indices[1][8] = {{0,1,2,3,4,5,6,7}};

  const int nelem = 1;
  int nnode = 8;

  TPZGeoEl *elvec[nelem];
  TPZGeoMesh *gmesh = new TPZGeoMesh();

  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZVec<REAL> coord(3);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    coord[2] = co[nod][2];
    TPZGeoNode pznode (nod,coord,*gmesh);
    gmesh->NodeVec()[nodind] = pznode;
  }

  int el;
  for(el=0; el<nelem; el++) {
    TPZVec<int> nodind(8);
    for(nod=0; nod<8; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
    //    elvec[el] = new TPZGeoElC3d(el,nodind,1,*gmesh);
  }

  gmesh->BuildConnectivity();

  //TPZVec<TPZGeoEl *> sub;
  //elvec[0]->Divide(sub);
  //   	elvec[1]->Divide(sub);
  //   	elvec[2]->Divide(sub);

//  TPZGeoElBC gbc;

  // bc -1 -> Dirichlet at the bottom face of the cube
  TPZGeoElBC gbc1(elvec[0],20,-1,*gmesh);

  // bc -2 -> Neumann at the top face of the cube
  TPZGeoElBC gbc2(elvec[0],25,-2,*gmesh);


  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

  TPZAutoPointer<TPZMaterial> mat;
    mat = new TPZMaterialTest3D(1);
    TPZFMatrix mp (3,1,0.);
    TPZMaterialTest3D * mataux = dynamic_cast<TPZMaterialTest3D *> (mat.operator ->());
    TPZMaterialTest3D::geq3=1;
    mataux->SetMaterial(mp);
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);
  TPZAutoPointer<TPZMaterial> bc[2];
  bc[0] = mat->CreateBC(mat,-1,0,val1,val2);
  int i;

  val2(0,0)=1.;
  bc[1] = mat->CreateBC(mat,-2,1,val1,val2);

  cmesh->InsertMaterialObject(mat);
  for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);

  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  //cmesh->ExpandSolution();

  return cmesh;
}
