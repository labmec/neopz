#include "includes.h"

static TPZCompMesh *CreateTestMesh();

//*************************************
//************Option 8*****************
//*****All element types Mesh**********
//*************************************
TPZCompMesh * CreateTestMesh() {

  REAL nodeco[12][3] = {
    {0.,0.,0.},
    {1.,0.,0.},
    {2.,0.,0.},
    {0.,1.,0.},
    {1.,1.,0.},
    {2.,1.,0.},
    {0.,0.,1.},
    {1.,0.,1.},
    {2.,0.,1.},
    {0.,1.,1.},
    {1.,1.,1.},
    {2.,1.,1.}
  };

  int nodind[7][8] = {
    {0,1,4,3,6,7,10,9},
    {2,4,10,8,5},
    {8,10,11,5},
    {2,4,1,8,10,7},
    {0,1},
    {0,1,7,6},
    {1,2,7}
  };

  int numnos[7] = {8,5,4,6,2,4,3};

  TPZGeoMesh *gmesh = new TPZGeoMesh();

  int noind[12];
  int no;
  for(no=0; no<12; no++) {
    noind[no] = gmesh->NodeVec().AllocateNewElement();
    TPZVec<REAL> coord(3);
    coord[0] = nodeco[no][0];
    coord[1] = nodeco[no][1];
    coord[2] = nodeco[no][2];
    gmesh->NodeVec()[noind[no]].Initialize(coord,*gmesh);
  }
  int matid = 1;
  TPZVec<int> nodeindex;
  int nel;
  TPZVec<TPZGeoEl *> gelvec;
  gelvec.Resize(4);
  for(nel=0; nel<4; nel++) {
    int in;
    nodeindex.Resize(numnos[nel]);
    for(in=0; in<numnos[nel]; in++) {
      nodeindex[in] = nodind[nel][in];
    }
    int index;
    switch(nel) {
    case 0:
      //      elvec[el] = gmesh->CreateGeoElement(ECube,nodeindex,1,index);
//      gelvec[nel]=new TPZGeoElC3d(nodeindex,matid,*gmesh);
      break;
    case 1:
      gelvec[nel] = gmesh->CreateGeoElement(EPiramide,nodeindex,matid,index);
      //       gelvec[nel]=new TPZGeoElPi3d(nodeindex,matid,*gmesh);
      break;
    case 2:
      gelvec[nel] = gmesh->CreateGeoElement(ETetraedro,nodeindex,matid,index);
    //       gelvec[nel]=new TPZGeoElT3d(nodeindex,matid,*gmesh);
      break;
    case 3:
//       gelvec[nel]=new TPZGeoElPr3d(nodeindex,matid,*gmesh);
//      gelvec[nel] = gmesh->CreateGeoElement(EPrisma,nodeindex,matid,index);
      break;
    case 4:
      //      gelvec[nel]=new TPZGeoEl1d(nodeindex,2,*gmesh);
      break;
    case 5:
      //      gelvec[nel]=new TPZGeoElQ2d(nodeindex,3,*gmesh);
      break;
    case 6:
      //      gelvec[nel]=new TPZGeoElT2d(nodeindex,3,*gmesh);
      break;
    default:
      break;
    }
  }
  gmesh->BuildConnectivity2();

  //TPZVec<TPZGeoEl *> sub;
  //elvec[0]->Divide(sub);
  //   	elvec[1]->Divide(sub);
  //   	elvec[2]->Divide(sub);

  TPZGeoElBC gbc;

  // bc -1 -> Dirichlet
//  TPZGeoElBC gbc1(gelvec[0],20,-1,*gmesh);
  TPZGeoElBC gbc11(gelvec[1],14,-1,*gmesh);
//  TPZGeoElBC gbc12(gelvec[3],15,-1,*gmesh);



  // bc -2 -> Neumann at the right x==1
//  TPZGeoElBC gbc2(gelvec[0],25,-2,*gmesh);
//  TPZGeoElBC gbc21(gelvec[3],19,-2,*gmesh);
  TPZGeoElBC gbc22(gelvec[2],10,-2,*gmesh);

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

  TPZMaterial *mat;
  //  if(nstate == 3) {
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
  val2(0,0) = 1.;
  bc[1] = mat->CreateBC(-2,1,val1,val2);
  cmesh->InsertMaterialObject(mat);

  int i;
  for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);

  gmesh->Print(cout);

  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();

  gmesh->Print(cout);
  return cmesh;
}
