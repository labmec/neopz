#include "includes.h"
#include "pzbndcond.h"

static TPZCompMesh *Create3DMesh();


//*************************************
//************Option 5*****************
//********Corner 3D Cube Mesh**********
//*************************************
TPZCompMesh *Create3DMesh() {

  REAL co[26][3] = {
    {0.,0.,0.},{1.,0.,0.},{2.,0.,0.},
    {0.,1.,0.},{1.,1.,0.},{2.,1.,0.},
    {0.,2.,0.},{1.,2.,0.},{2.,2.,0.},
    {0.,0.,1.},{1.,0.,1.},{2.,0.,1.},
    {0.,1.,1.},{1.,1.,1.},{2.,1.,1.},
    {0.,2.,1.},{1.,2.,1.},{2.,2.,1.},
    {0.,0.,2.},{1.,0.,2.},{2.,0.,2.},
    {0.,1.,2.},{1.,1.,2.},{2.,1.,2.},
    {0.,2.,2.},{1.,2.,2.}
  };
  
  int indices[7][8] = {
    {0,1,4,3,9,10,13,12},
    {1,2,5,4,10,11,14,13},
    {3,4,7,6,12,13,16,15},
    {4,5,8,7,13,14,17,16},
    {9,10,13,12,18,19,22,21},
    {10,11,14,13,19,20,23,22},
    {12,13,16,15,21,22,25,24}
  };
  
  const int nelem = 7;
  int nnode = 26;

  TPZGeoEl *elvec[nelem];
  TPZGeoMesh *gmesh = new TPZGeoMesh();

  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZVec<REAL> coord(3);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    coord[2] = co[nod][2];
    TPZGeoNode pznode(nod,coord,*gmesh);
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

  gmesh->BuildConnectivity2();
	
  //  	TPZVec<TPZGeoEl *> sub;
  //   	elvec[0]->Divide(sub);
  //   	elvec[1]->Divide(sub);
  //   	elvec[2]->Divide(sub);
  
//  TPZGeoElBC gbc;

  // bc -1 -> Dirichlet
  TPZGeoElBC gbc1(elvec[0],21,-1,*gmesh);
  
  TPZGeoElBC gbc2(elvec[0],22,-1,*gmesh);
  
  TPZGeoElBC gbc3(elvec[0],25,-1,*gmesh);

  // bc -3 -> Neumann at the right x==1
  TPZGeoElBC gbc4(elvec[3],25,-2,*gmesh);

  // bc -4 -> Neumann at the top y==1
  TPZGeoElBC gbc5(elvec[5],23,-3,*gmesh);

  // bc -4 -> Neumann at the top y==1
  TPZGeoElBC gbc6(elvec[6],22,-4,*gmesh);

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  int nstate = 3;
  TPZAutoPointer<TPZMaterial> mat;
  //  if(nstate == 2) {
  //    mat = new TPZMatHyperElastic(1,2.,400);
  //  } else {
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
    //  }
  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  TPZAutoPointer<TPZMaterial> bc[4];
  bc[0] = mat->CreateBC(mat,-1,0,val1,val2);
  int i;
  val2(2,0)=-1.;
  bc[1] = mat->CreateBC(mat,-2,1,val1,val2);
  val2(2,0)=0.;
  val2(1,0)=-1.;
  bc[2] = mat->CreateBC(mat,-3,1,val1,val2);
  val2(1,0)=0.;
  val2(0,0)=-1.;
  bc[3] = mat->CreateBC(mat,-4,1,val1,val2);
   
  cmesh->InsertMaterialObject(mat);
  for(i=0; i<4; i++) cmesh->InsertMaterialObject(bc[i]);

  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  //cmesh->ExpandSolution();
  
  return cmesh;
}
