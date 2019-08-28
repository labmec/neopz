#include "includes.h"

static TPZCompMesh *Create3DPrismMesh();

//*************************************
//************Option 7*****************
//**********3D Prism Mesh**************
//*************************************
TPZCompMesh *Create3DPrismMesh() {

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

  TPZGeoElBC gbc;

  // bc -1 -> Neumann
  TPZGeoElBC gbc1(elvec[0],15,-1,*gmesh);

  // bc -2 -> Neumann at the right x==1
  TPZGeoElBC gbc2(elvec[0],19,-2,*gmesh);

  // bc -2 -> Dirichlet at LOWER FACE
  TPZGeoElBC gbc3(elvec[0],15,-3,*gmesh);


  gmesh->BuildConnectivity2();
  //ofstream MALHAG("malhageometrica");
  //  gmesh->Print(MALHAG);

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

  TPZMaterial *mat;
  //  if(nstate == 3) {
    //		mat = new TPZMatHyperElastic(1,2.,400);
    mat = new TPZMaterialTest3D(1);
    TPZFMatrix mp (1,1,0.);

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
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);
  TPZBndCond *bc[3];
  bc[0] = mat->CreateBC(-3,0,val1,val2);
  int i;
  val2(0,0)=-1.;
  bc[1] = mat->CreateBC(-2,1,val1,val2);
  val2(0,0)=1.;
  bc[2] = mat->CreateBC(-1,1,val1,val2);

  cmesh->InsertMaterialObject(mat);
  for(i=0; i<3; i++) cmesh->InsertMaterialObject(bc[i]);

  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();

  cmesh->Print(cout);

/*   int o; */
/*   for (o=0;o<10;o++){ */
/*     TPZIntPrism3D prismrule(2*o+2,2*o+2); */
/*     TPZVec<int> ord (3,2*o+2); */
/*     prismrule.SetOrder(ord); */
/*     int np = prismrule.NPoints(); */
/*     int p; */
/*     cout << endl << endl <<"Ordem o = " << 2*o+2 << endl; */
/*     for (p=0;p<np;p++){ */
/*       TPZVec<REAL> loc(3,0.); */
/*       REAL weight = -1.; */
/*       prismrule.Point(p,loc,weight); */
/*       cout << "Point " << p << "  (x,y,z) = " << loc[0] << " , "   */
/* 	   << loc[1] << " , "  << loc[2] << "  weight = "  << weight << endl; */
/*     } */
/*   } */

  return cmesh;
}
