#include "includes.h"

static TPZCompMesh *Create3DTetraMesh();

//*************************************
//************Option 6*****************
//*********3D Tetrahedra Mesh**********
//*************************************
TPZCompMesh *Create3DTetraMesh() {

  REAL co[4][3] = {
    {0.,0.,0.},
    {1.,0.,0.},
    {0.,1.,0.},
    {0.,0.,1.}
  };

  int indices[1][4] = {{0,1,2,3}};
  const int nelem = 1;
  int nnode = 4;

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
    TPZVec<int> nodind(4);
    for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(ETetraedro,nodind,1,index);

    //    elvec[el] = new TPZGeoElT3d(el,nodind,1,*gmesh);
  }

  // bc -1 -> Dirichlet
  TPZGeoElBC gbc1(elvec[0],10,-1,*gmesh);

  // bc -2 -> Neumann at the right x==1
  TPZGeoElBC gbc2(elvec[0],13,-2,*gmesh);

  gmesh->BuildConnectivity2();

  TPZStack<TPZGeoEl*> subel;
//  elvec[0]->Divide(subel);
//  subel[0]->Divide(subel);
  //  TPZGeoElBC gbc;

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

  TPZMaterial *mat;
  //  if(nstate == 3) {
    //		mat = new TPZMatHyperElastic(1,2.,400);
    mat = new TPZMaterialTest3D(1);
    TPZFMatrix mp (3,1,1.);

    TPZMaterialTest3D * mataux = dynamic_cast<TPZMaterialTest3D *> (mat);
    TPZMaterialTest3D::geq3=1;
    mataux->SetMaterial(mp);
    /*  } else {
    TPZMat2dLin *mat2d = new TPZMat2dLin(1);
    int ist,jst;
    TPZFMatrix xk(nstate,nstate,1.),xc(nstate,nstate,0.),xf(nstate,1,0.);
    for(ist=0; ist<nstate; ist++) {
      if(nstate != 1) xf(ist,0) = 1.;
      for(jst=0; jst<nstate; jst++) {
	if(ist != jst) xk(ist,jst) = 0.;
      }
    }
    mat2d->SetMaterial(xk,xc,xf);
    mat = mat2d;
    }*/
  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  TPZBndCond *bc[2];
  bc[0] = mat->CreateBC(-1,0,val1,val2);
  int i;

  val2(2,0)=-1.;
  bc[1] = mat->CreateBC(-2,1,val1,val2);

  cmesh->InsertMaterialObject(mat);
  for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);

  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  //cmesh->ExpandSolution();

  return cmesh;
}
