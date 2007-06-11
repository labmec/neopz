
#include "includes.h"

static TPZCompMesh *Create3DDiscMesh();
//static void BCSolution(TPZVec<REAL> &x,TPZVec<REAL> &result);
//static void Solution(TPZVec<REAL> &x,TPZVec<REAL> &result,TPZFMatrix &deriv);
//static void LoadSolution(TPZFMatrix &axes,TPZVec<REAL> &X,TPZFMatrix &u,TPZFMatrix &du);


//*************************************
//************Option 12*****************
//*******3D Discontinuous Mesh**********
//*************************************
TPZCompMesh *Create3DDiscMesh() {

  REAL co[26][3] = {
    {0.,0.,0.},{0.5,0.,0.},{1,0.,0.},{0.,0.5,0.},{0.5,0.5,0.},{1,0.5,0.},{0.,1.,0.},{0.5,1.,0.},{1,1.,0.},
    {0.,0.,0.5},{0.5,0.,0.5},{1,0.,0.5},{0.,0.5,0.5},{0.5,0.5,0.5},{1,0.5,0.5},{0.,1.,0.5},{0.5,1.,0.5},{1,1.,0.5},
    {0.,0.,1.},{0.5,0.,1.},{1,0.,1.},{0.,0.5,1.},{0.5,0.5,1.},{1,0.5,1.},{0.,1.,1.},{0.5,1.,1.}
  };

  int indices[7][8] = {
    {0,1,4,3,9,10,13,12},
    {1,2,5,4,10,11,14,13},
    {3,4,7,6,12,13,16,15},
    {4,5,8,7,13,14,17,16},
    {9,10,13,12,18,18,22,21},
    {10,11,14,13,19,20,23,22},
    {12,13,16,15,21,22,25,24}
  };

  const int nelem = 7;
  int nnode = 26;
  const int nodeperel = 8;

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
    TPZVec<int> nodind(nodeperel);
    for(nod=0; nod<nodeperel; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
    //    elvec[el] = new TPZGeoElPr3d(el,nodind,1,*gmesh);
  }


//  TPZGeoElBC gbc;

  //Condicoes de Neumann
  // bc -1 -> Face inferior
  TPZGeoElBC gbc1(elvec[0],20,-1,*gmesh);
  TPZGeoElBC gbc2(elvec[1],20,-1,*gmesh);
  TPZGeoElBC gbc3(elvec[2],20,-1,*gmesh);
  TPZGeoElBC gbc4(elvec[3],20,-1,*gmesh);

  // bc -2 -> Face lateral esquerda
  TPZGeoElBC gbc5(elvec[0],24,-2,*gmesh);
  TPZGeoElBC gbc6(elvec[2],24,-2,*gmesh);
  TPZGeoElBC gbc7(elvec[4],24,-2,*gmesh);
  TPZGeoElBC gbc8(elvec[6],24,-2,*gmesh);

  // bc -3 -> Face frontal
  TPZGeoElBC gbc9(elvec[0],21,-3,*gmesh);
  TPZGeoElBC gbc10(elvec[1],21,-3,*gmesh);
  TPZGeoElBC gbc11(elvec[4],21,-3,*gmesh);
  TPZGeoElBC gbc12(elvec[5],21,-3,*gmesh);

  // bc -3 -> Face lateral direita
  TPZGeoElBC gbc13(elvec[1],22,-4,*gmesh);
  TPZGeoElBC gbc14(elvec[3],22,-4,*gmesh);
  TPZGeoElBC gbc15(elvec[5],22,-4,*gmesh);

  // bc -3 -> Face posterior
  TPZGeoElBC gbc16(elvec[2],23,-5,*gmesh);
  TPZGeoElBC gbc17(elvec[3],23,-5,*gmesh);
  TPZGeoElBC gbc18(elvec[6],23,-5,*gmesh);

  //Condicoes Dirichlet
  TPZGeoElBC gbc19(elvec[3],25,-6,*gmesh);
  TPZGeoElBC gbc20(elvec[5],23,-6,*gmesh);
  TPZGeoElBC gbc21(elvec[6],22,-6,*gmesh);

  gmesh->BuildConnectivity2();
  //ofstream MALHAG("malhageometrica");
  //gmesh->Print(MALHAG);
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

  TPZMaterialTest3D *mat3D = new TPZMaterialTest3D(1);
  TPZFMatrix mp (3,1,0.);
  TPZMaterialTest3D::geq3=1;
  mat3D->SetMaterial(mp);
  TPZAutoPointer<TPZMaterial> mat(mat3D);

  //  mat->SetForcingFunction(ForcingFunction3DExp);

  //CreateBC
  //....


  cmesh->InsertMaterialObject(mat);
  //int i;
  // for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);

  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();


  TPZVec <int> subelvec;
  cmesh->ElementVec()[0]->Divide(0,subelvec,1);

  //int  pord = 3;
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
