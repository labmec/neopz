#include "pzvec.h"

#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzgnode.h"

#include "pzmatrix.h"

#include "pzelasmat.h"
//#include "pzplaca.h"
#include "pzmat2dlin.h"
//#include "pzmathyperelastic.h"
#include "pzmattest.h"
//#include "pzmatplaca2.h"

static TPZCompMesh * CreateSillyMesh();

//*************************************
//************Option 0*****************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateSillyMesh(){
  //malha quadrada de nr x nc
  const	int numrel = 1;
  const	int numcel = 1;
  //  int numel = numrel*numcel;
  TPZVec<REAL> coord(2,0.);
  
  // criar um objeto tipo malha geometrica
  TPZGeoMesh *geomesh = new TPZGeoMesh();
  
  // criar nos
  int i,j;
  for(i=0; i<(numrel+1); i++) {
    for (j=0; j<(numcel+1); j++) {
      int nodind = geomesh->NodeVec().AllocateNewElement();
      TPZVec<REAL> coord(2);
      coord[0] = j;
      coord[1] = i;
      geomesh->NodeVec()[nodind] = TPZGeoNode(i*(numcel+1)+j,coord,*geomesh);
    }
  }

  // criação dos elementos
  int elc, elr;
  TPZGeoEl *gel[numrel*numcel];
  TPZVec<int> indices(4);
  for(elr=0; elr<numrel; elr++) {  
    for(elc=0; elc<numcel; elc++) {
      indices[0] = (numcel+1)*elr+elc;
      indices[1] = indices[0]+1;
      indices[3] = indices[0]+numcel+1;
      indices[2] = indices[1]+numcel+1;
      // O proprio construtor vai inserir o elemento na malha
      //      gel[elr*numcel+elc] = new TPZGeoElQ2d(elr*numcel+elc,indices,1,*geomesh);
      int index;
      gel[elr*numcel+elc] = geomesh->CreateGeoElement(EQuadrilateral,indices,1,index);
    }
  }
 
  // Descomentar o trecho abaixo para habilitar a
  // divisão dos elementos geométricos criados 
  // Criação das condições de contorno geométricas
  TPZGeoElBC t3(gel[0],4,-1,*geomesh);

  TPZGeoElBC t4(gel[0],6,-2,*geomesh);

  geomesh->BuildConnectivity2();


  //Divisão dos elementos
  TPZVec<TPZGeoEl *> sub;
  gel[0]->Divide(sub);
  /* for (i=0;i<(sub.NElements()-1);i++){
    TPZVec<TPZGeoEl *> subsub;
    sub[i]->Divide(subsub);
  }
  */

  geomesh->Print(cout);
   
  // Criação da malha computacional
  TPZCompMesh *comp = new TPZCompMesh(geomesh);

  // Criar e inserir os materiais na malha
  TPZMaterial *mat = new TPZElasticityMaterial(1,1.e5,0.2,0,0);
  //TPZMaterial *mat = new TPZMaterialTest(1,1.,1.);
  comp->InsertMaterialObject(mat);
 
  TPZMaterial *meumat = mat;

  // Condições de contorno
  // Dirichlet
  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  TPZMaterial *bnd = meumat->CreateBC (-1,0,val1,val2);
  comp->InsertMaterialObject(bnd);
  bnd = meumat->CreateBC (-2,0,val1,val2);

  // bnd->SetForcingFunction(Forcing1);
  // comp->InsertMaterialObject(bnd);

  // Neumann
  TPZFMatrix val3(3,3,1);
  val2(0,0)=1.;
  bnd = meumat->CreateBC (-2,1,val1,val2);
  comp->InsertMaterialObject(bnd);
  
  // comp->Print(cout);

  // Ajuste da estrutura de dados computacional
  comp->AutoBuild();
  //  comp->Print(cout);
  //comp->AdjustBoundaryElements();
  //  comp->Print(cout);
  //comp->CleanUpUnconnectedNodes();
  //  comp->Print(cout);

/*  //	comp->Print(output); */
/*  TPZInterpolatedElement *intel = dynamic_cast <TPZInterpolatedElement *> (comp->ElementVec()[1]); */
/*  TPZVec<int> subelindex; */
/*  intel->Divide(1,subelindex,1); */
/*  int isub; */
/*  int nsides = intel->NConnects(); */
/*  int porder = intel->PreferredSideOrder(nsides-1); */
/*   for (isub=0; isub<subelindex.NElements();isub++){ */
/*     TPZInterpolatedElement *cintel = dynamic_cast<TPZInterpolatedElement *> (comp->ElementVec()[subelindex[isub]]); */
/*     cintel->PRefine(porder+1); */
/*   } */
/*comp->ExpandSolution(); */

    comp->SetName("Malha Computacional Original");
    //   comp->Print(cout);
    //    cout << endl << "Number of equations: " << comp->NEquations() << endl;
    // cout.flush();
    return comp;
}
