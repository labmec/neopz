#include "includes.h"

static TPZCompMesh *CreateAleatorioMesh();

 //*************************************
//************Option 10*****************
//**********Aleatorio Mesh**************
//*************************************
TPZCompMesh *CreateAleatorioMesh() {

  REAL co[6][3] = {
    {0.,0.,0.},
    {1.,0.,0.},
    {0.,1.,0.},
    {0.,0.,1.},
    {1.,0.,1.},
    {0.,1.,1.}
  };

  int indices[1][6] = {{0,1,2,3,4,5}};

  const int nelem = 1;
  int nnode = 6;

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
    TPZVec<int> nodind(6);
    for(nod=0; nod<6; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(EPrisma,nodind,1,index);
    //    elvec[el] = new TPZGeoElPr3d(el,nodind,1,*gmesh);
  }

  //  TPZStack<TPZGeoEl*> subel;
  //  elvec[0]->Divide(subel);

  //  TPZGeoElBC gbc;

  // bc -1 -> Dirichlet
  TPZGeoElBC gbc1(elvec[0],15,-1,*gmesh);

  // bc -2 -> Neumann at the right x==1
  TPZGeoElBC gbc2(elvec[0],19,-2,*gmesh);

  gmesh->BuildConnectivity2();
  //ofstream MALHAG("malhageometrica");
  //  gmesh->Print(MALHAG);

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

  TPZMaterial *mat;
  //  if(nstate == 3) {
    //		mat = new TPZMatHyperElastic(1,2.,400);
    mat = new TPZMaterialTest3D(1);
    TPZFMatrix mp (3,1,0.);

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
  val2(0,0)=-1.;
  bc[1] = mat->CreateBC(-2,0,val1,val2);

  cmesh->InsertMaterialObject(mat);
  for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);

  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();


  TPZVec <int> subelvec;
  cmesh->ElementVec()[0]->Divide(0,subelvec,1);

//  int  pord = 3;
  cmesh->ElementVec()[subelvec [7]]->Divide(subelvec [7],subelvec,1);
  TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cmesh->ElementVec()[subelvec [7]]);
  intel->PRefine(2);
  cmesh->ElementVec()[subelvec [7]]->Divide(subelvec [7],subelvec,1);
  TPZInterpolatedElement *intel1 = dynamic_cast<TPZInterpolatedElement *> (cmesh->ElementVec()[subelvec [7]]);
  intel1->PRefine(3);

 //  cmesh->ElementVec()[subelvec [7]]->PRefine(4);
  //  cmesh->ElementVec()[0]->Divide(,subelvec,1);


  //  cmesh->Print(cout);
  return cmesh;
}
