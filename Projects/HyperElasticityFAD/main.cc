#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgnode.h"
#include "pzmaterial.h"
//#include "pzerror.h"
#include "pzgeoel.h"
//#include "pzcosys.h"
#include "pzmatrix.h"
#include "pzelgq2d.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzintel.h"
#include "pzelcq2d.h"
#include "pzskylstrmatrix.h"
#include "pzpoisson3d.h"
#include "pzcheckgeom.h"
#include <time.h>
#include <stdio.h>
#include "pzmathyperelastic.h"
#include "pzelgc3d.h"

TPZCompMesh *CreateMesh();


void error(char * err)
{
  PZError << "FADERROR: " << err << endl;
};

int main(){

  TPZCompMesh *cmesh = CreateMesh();

      TPZAnalysis an (cmesh);
      TPZSkylineStructMatrix strskyl(cmesh);
      an.SetStructuralMatrix(strskyl);

      TPZStepSolver direct;
      direct.SetDirect(ECholesky);
      an.SetSolver(direct);

      an.Run();
      an.Rhs().Print();
      an.Solution().Print();
  TPZMatrixSolver::Diagnose();
  return 0;
}



//*************************************
//************Option 0*****************
//*******Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateMesh(){

  //malha quadrada de nr x nc
  const	int numrel = 1;
  const	int numcel = 1;
  const int numzel = 1;
  //  int numel = numrel*numcel;
  TPZVec<REAL> coord(3,0.);

  // criar um objeto tipo malha geometrica
  TPZGeoMesh *geomesh = new TPZGeoMesh();

  // criar nos
  int i,j,k;
  for(k=0;k<numzel+1;k++)
  for(i=0; i<(numrel+1); i++) {
    for (j=0; j<(numcel+1); j++) {
      int nodind = geomesh->NodeVec().AllocateNewElement();
      TPZVec<REAL> coord(3);
      coord[0] = j;
      coord[1] = i;
      coord[2] = k;
      geomesh->NodeVec()[nodind] = TPZGeoNode(k*(numzel+1)*(numcel+1)+i*(numcel+1)+j,coord,*geomesh);
    }
  }

  TPZVec<int> indices(8);
  // criação dos elementos
  TPZGeoEl *gel[1];

      indices[0] = 0;
      indices[1] = 1;
      indices[3] = 2;
      indices[2] = 3;
      indices[4] = 4;
      indices[5] = 5;
      indices[7] = 6;
      indices[6] = 7;
      // O proprio construtor vai inserir o elemento na malha
      gel[0] = new TPZGeoElC3d(0,indices,1,*geomesh);

  //Divisão dos elementos
  TPZVec<TPZGeoEl *> sub;

  geomesh->BuildConnectivity();
  //geomesh->Print(cout);

//  gel[0]->Divide(sub);

  // Criação das condições de contorno geométricas
  TPZGeoElBC t3(gel[0],4,-1,*geomesh);
  TPZGeoElBC t4(gel[0],6,-2,*geomesh);


  // Criação da malha computacional
  TPZCompMesh *comp = new TPZCompMesh(geomesh);

  // Criar e inserir os materiais na malha
  TPZMatHyperElastic *mat = new
  TPZMatHyperElastic(1,/*REAL e*/1e5,/*REAL mu*/.25/*,REAL nu=-1.,REAL lambda=-1.,REAL coef1=-1.,REAL coef2=-1.,REAL coef3=-1.*/);

  comp->InsertMaterialObject(mat);

  TPZMaterial *meumat = mat;

  // Condições de contorno
  // Dirichlet
  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  TPZMaterial *bnd = meumat->CreateBC (-1,0,val1,val2);
  comp->InsertMaterialObject(bnd);

  // Neumann
  val2(0,0)=7.;
  val2(1,0)=3.;
  val2(2,0)=5.;
  bnd = meumat->CreateBC (-2,1,val1,val2);
  comp->InsertMaterialObject(bnd);

  comp->Reference()->Print();
  comp->Print(cout);

  // Ajuste da estrutura de dados computacional
  comp->AutoBuild();


  comp->AdjustBoundaryElements();
  comp->CleanUpUnconnectedNodes();
    comp->Print(cout);
    comp->SetName("Malha Computacional Original");
    return comp;
}

