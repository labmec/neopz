//$Id: main.cpp,v 1.14 2005-01-14 16:47:15 tiago Exp $

/**
 * Galerkin descontinuo: visita do professor Igor.
 * 24/11/2003
 */

//pra acertar as condicoes de contorno
//DIFUSAO_EXP eh difusao pura
//DESCONTINUIDADE eh conveccao-difusao sem vetor de carga
//DEBUGM nao pergunta parametros h e p, nem nome dos arquivos
#define DIFUSAO_EXP
//#define DESCONTINUIDADE
//#define DEBUGM


#include "pzvec.h"

#include "pzcmesh.h"

#include "pzdebug.h"
#include "pzcheckgeom.h"
//#include "pzerror.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"

//#include "pzintel.h"
#include "pzcompel.h"
#include "TPZCompElDisc.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "TPZCopySolve.h"


#include "pzbstrmatrix.h"
#include "pzstepsolver.h"    		
#include "pzonedref.h"


#include "pzbndcond.h"
#include "pzpoisson3d.h"

#include "pzvisualmatrix.h"

#include <time.h>
#include <stdio.h>

static REAL PI, gDif;

static int gMeshType;

int gPrintLevel = 0;

int gDivide[6];

/*void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp) {
  disp[0] = -2.0 * PI * PI * sin( PI * x[0]) * sin(PI * x[1]) ;
}
void ExactSolution(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  u[0] = sin(PI*x[0])*sin(PI*x[1]);
  deriv(0,0) = PI*cos(PI*x[0])*sin(PI*x[1]);
  deriv(1,0) = PI*cos(PI*x[1])*sin(PI*x[0]); 
}**/

void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp) {
  disp[0] = - exp(0.75 * (x[0] + x[1])) * (8. * (1. - x[1] * x[1]) + 12. * x[0] * (1. - x[1] * x[1]) -4.5 * (1. - x[0] * x[0])* (1. - x[1] * x[1] )+
					  8. * (1. - x[0] * x[0]) + 12. * x[1] * (1. - x[0] * x[0]));
}
void ExactSolution(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  u[0] = 4. * (1. - x[0] * x[0] ) * (1. - x[1] * x[1]) * exp(0.75 * (x[0] + x[1]));
  deriv(0,0) = (-8. * x[0] * (1. - x[1] * x[1]) + 3. * (1. - x[0] * x[0] ) * (1. - x[1] * x[1]) ) * exp(0.75 * (x[0] + x[1]));
  deriv(1,0) = (-8. * x[1] * (1. - x[0] * x[0]) + 3. * (1. - x[0] * x[0] ) * (1. - x[1] * x[1]) ) * exp(0.75 * (x[0] + x[1]));
}

void Dirichlet1(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] = exp(-(x[0]*x[0]+1));
}

void Dirichlet2(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] = exp(-(x[1]*x[1]+1));
}

/*void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp) {
  disp[0] =  PI * PI * sin(PI * (x[1]+1.)/2.) / 4.;
}
void ExactSolution(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  u[0] = sin(PI*(x[1]+1.)/2.);
  deriv(0,0) = 0.;
  deriv(1,0) = PI*cos(PI*(x[1]+1.)/2.)/2.;
}
*/
/*void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp) {
  //  disp[0] =  -2.;
  disp[0] = 0.;
}
void ExactSolution(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  u[0] = (x[1]+1.);
  deriv(0,0) = 0.;
  deriv(1,0) = 1.;
}*/

//1 quadrado
TPZCompMesh * SimpleMesh();

//4 quadrados iguais, podendo ser divididos
TPZCompMesh *CreateMesh();

//8 triangulos
TPZCompMesh *CreateMesh2();
//8 triangulos alinhados com a descontinuidade
TPZCompMesh *CreateMesh3();

//malha de triangulos do Philippe
TPZCompMesh *CreateMeshPhil();
 

////////////////////////////////////////////////////////////////////////////////


int main(){
   
//   TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);

  int p, h;
  char filename[20];
  char filedx[20];

#ifdef DEBUGM
  p = 8;
  h = 2;
#endif

  cout << "\nQuadrado = 1; Triang = 2" << endl;
  cin >> gMeshType;

#ifndef DEBUGM
  cout << "\nOrdem p" << endl;
  cin >> p;
  cout << "\nRefinamento" << endl;
  cin >> h;
#endif


#ifdef DESCONTINUIDADE
  cout << "\nDifusao" << endl;
  //  cin >> gDif; 
  gDif = 0.001;
#endif

#ifdef DIFUSAO_EXP
 gDif = 1.0;
#endif

#ifdef DEBUGM
  strcpy(filename,"teste.dat");  
  strcpy(filedx,"teste.dx");
#endif

#ifndef DEBUGM
  cout << "\nArquivo" << endl;
  cin >> filename;
  cout << "\nArquivoDX" << endl;
  cin >> filedx;
#endif

  if (h == 0) {
    gDivide[0] = 0;
    gDivide[1] = 0;
    gDivide[2] = 0;
    gDivide[3] = 0;
    gDivide[4] = 0;
    gDivide[5] = 0;
  }
  if (h == 1) {
    gDivide[0] = 1;
    gDivide[1] = 0;
    gDivide[2] = 0;
    gDivide[3] = 0;
    gDivide[4] = 0;
    gDivide[5] = 0;
  } 
  if (h == 2) {
    gDivide[0] = 1;
    gDivide[1] = 1;
    gDivide[2] = 0;
    gDivide[3] = 0;
    gDivide[4] = 0;
    gDivide[5] = 0;
  }  
  if (h == 3) {
    gDivide[0] = 1;
    gDivide[1] = 1;
    gDivide[2] = 1;
    gDivide[3] = 0;
    gDivide[4] = 0;
    gDivide[5] = 0;
  }
  if (h == 4) {
    gDivide[0] = 1;
    gDivide[1] = 1;
    gDivide[2] = 1;
    gDivide[3] = 1;
    gDivide[4] = 0;
    gDivide[5] = 0;
  }
  if (h == 5) {
    gDivide[0] = 1;
    gDivide[1] = 1;
    gDivide[2] = 1;
    gDivide[3] = 1;
    gDivide[4] = 1;
    gDivide[5] = 0;
  }
  if (h == 6) {
    gDivide[0] = 1;
    gDivide[1] = 1;
    gDivide[2] = 1;
    gDivide[3] = 1;
    gDivide[4] = 1;
    gDivide[5] = 1;
  }


  PI = 4.*atan(1.);

  ofstream out(filename);

  TPZCompEl::gOrder = p;
  TPZCompElDisc::gDegree = p;

//  gDebug = 0;

  TPZCompMesh *cmesh;
  if (gMeshType == 1) //Quadrado
    cmesh = CreateMesh();
//     cmesh = SimpleMesh();
  if (gMeshType == 2) //Triangulo
    //    cmesh = CreateMesh2();
    cmesh = CreateMesh3();

  TPZGeoMesh *gmesh = cmesh->Reference();

  TPZAnalysis an(cmesh);

  //METODOS DE RESOLUCAO: ITERATIVO OU DIRETO

#define ITER_SOLVER
//_PRECOND
//#define DIRECT_SOLVER

#ifdef ITER_SOLVER

  cout << "ITER_SOLVER" << endl;  
  TPZSpStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);  
  TPZCopySolve precond( full.Create()  );
  TPZStepSolver step;
  step.ShareMatrix( precond );  
  step.SetGMRES( 2000000, 30, precond, 1.e-13, 0 ); 
  //step.SetCG( 2000000, precond, 1e-13, 0);
  an.SetSolver(step);

#endif

#ifdef ITER_SOLVER_PRECOND
  cout << "ITER_SOLVER_PRECOND" << endl;
  TPZSpStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);
  
  TPZBlockDiagonalStructMatrix strblock(cmesh);
  TPZBlockDiagonal * block = new TPZBlockDiagonal();
  strblock.AssembleBlockDiagonal(*block);
  
  TPZStepSolver step, precond(block);
  precond.SetDirect(ELU);

  step.SetGMRES( 2000000, 30, precond, 1.e-13, 0); 
  an.SetSolver(step);

#endif

#ifdef DIRECT_SOLVER

  cout << "DIRECT_SOLVER" << endl;
  TPZParFrontStructMatrix <TPZFrontNonSym> full(cmesh);
//  full.SetNumberOfThreads(3);
//  TPZFStructMatrix full(cmesh);  
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect(ELU);
  an.SetSolver(step);
#endif

{  
  TPZFMatrix fillin;
  cmesh->ComputeFillIn(50,fillin);
  //fillin.Print("Fillin of the computable mesh");
  VisualMatrix(fillin , filedx);
}

  an.Run();

/**** Aqui faz DX ****/
/*  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(1);
  scalnames[0] = "Solution";
  vecnames[0] = "Derivate";
  an.DefineGraphMesh(2,scalnames,vecnames,filedx);

  if (gMeshType == 1)
    an.PostProcess(2);

  if (gMeshType == 2)
    an.PostProcess(0);
*/
#ifdef DIFUSAO_EXP
  an.SetExact(ExactSolution);
  TPZVec<REAL> pos;
  an.PostProcess(pos,out);
  out << "\nNumero de equacoes: " << an.Solution().Rows() << endl;  
#endif
  
  delete cmesh;
  delete gmesh;
  return 0;
}




//////////////////////// COMECAM AS MALHAS:
/////////////////////// CREATEMESH -> 4 QUADRADOS
////////////////////// CREATEMESH2 -> 8 TRIANGULOS
///////////////////// CREATEMESH3 -> 8 TRIANGULOS ALINHADOS COM A DESCONTINUIDADE
//////////////////// CREATEMESHPHIL -> MALHA DE TRIANGULOS DO PHILIPPE




#include "TPZGeoElement.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"

TPZCompMesh *SimpleMesh(){
  REAL co[4][2] = {{1.,-1.},{1.,1.},{-1.,1.},{-1,-1}};
  int indices[1][4] = {{0,1,2,3}};
  TPZGeoEl *elvec[4];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 4;
  int nelem = 1;
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el=0; el<nelem; el++) {
    TPZManVector<int,4> nodind(4);
    for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
  }

  gmesh->BuildConnectivity();

  TPZGeoElBC gbc;

  TPZGeoElBC gbc1(elvec[0],4,-3,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],5,-3,*gmesh); // right
  TPZGeoElBC gbc3(elvec[0],6,-3,*gmesh); // right
  TPZGeoElBC gbc4(elvec[0],7,-3,*gmesh); // top
  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetForcingFunction(Forcing1);

  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;
  //Difusao pura. conveccao = 0.
  mat->SetParameters(gDif, 0., convdir);
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc;
  
  val2.Zero();
  bc = mat->CreateBC(-3, 0,val1,val2);

  
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bc);

  TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  //template class TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>;
  
  cmesh->AutoBuild();
  //  cmesh->AdjustBoundaryElements();
  //  cmesh->CleanUpUnconnectedNodes();
  //  cmesh->ExpandSolution();
  
  return cmesh;
}


TPZCompMesh *CreateMesh() {

  REAL co[9][2] = {{0.,0.},{0.,-1.},{1.,-1.},{1.,0.},{1.,1.},{0.,1.},{-1.,1.},{-1.,0.},{-1,-1}};
  int indices[4][4] = {{0,1,2,3},{0,3,4,5},{0,5,6,7},{0,7,8,1}};
  TPZGeoEl *elvec[4];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 9;
  int nelem = 4;
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el=0; el<nelem; el++) {
    TPZManVector<int,4> nodind(4);
    for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
  }

  gmesh->BuildConnectivity();

  for(int i = 0; i < nelem; i++){
    TPZVec<TPZGeoEl *> children, netos, bisnetos, tata1, tata2, tata3;
    //    cout << "\ngDivide[0] = \n" << gDivide[0];
    if(gDivide[0] == 1) {
      elvec[i]->Divide(children);
      //      cout <<  "\n Primeira divisao \n" ;
      if (gDivide[1] == 1) {
        for(int j = 0; j < children.NElements(); j++) {
          children[j]->Divide(netos);
	  //          cout <<  "\n Segunda divisao \n" ;
          if(gDivide[2] == 1) {
            for(int k = 0; k < netos.NElements(); k++) {
              netos[k]->Divide(bisnetos); 
	      //              cout <<  "\n Terceira divisao \n" ;
              if(gDivide[3] == 1) {
                for(int k = 0; k < bisnetos.NElements(); k++) {
                  bisnetos[k]->Divide(tata1); 
		  //                  cout <<  "\n Quarta divisao \n" ;
                  if(gDivide[4] == 1) {
                    for(int k = 0; k < tata1.NElements(); k++) {
                      tata1[k]->Divide(tata2); 
		      //                      cout <<  "\n Quinta divisao \n" ;
                      if(gDivide[5] == 1) {
                        for(int k = 0; k < tata2.NElements(); k++) tata2[k]->Divide(tata3); 
			//                        cout <<  "\n Sexta divisao \n" ;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }  
 
  TPZGeoElBC gbc;

#ifdef DESCONTINUIDADE
  TPZGeoElBC gbc1(elvec[0],5,-2,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],6,-2,*gmesh); // right
  TPZGeoElBC gbc3(elvec[1],5,-1,*gmesh); // right
  TPZGeoElBC gbc4(elvec[1],6,-1,*gmesh); // top
  TPZGeoElBC gbc5(elvec[2],5,-2,*gmesh); // top
  TPZGeoElBC gbc6(elvec[2],6,-2,*gmesh); // left
  TPZGeoElBC gbc7(elvec[3],5,-1,*gmesh); // left
  TPZGeoElBC gbc8(elvec[3],6,-1,*gmesh); // bottom
  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);

  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;
  mat->SetParameters(gDif, 1., convdir);
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc[2];
  
  val2(0,0) = 1.;
  bc[0] = mat->CreateBC(-1,0,val1,val2);

  val2.Zero();
  bc[1] = mat->CreateBC(-2,0,val1,val2);

  
  cmesh->InsertMaterialObject(mat);
  int i;
  for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);

#endif

#ifdef DIFUSAO_EXP
  TPZGeoElBC gbc1(elvec[0],5,-3,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],6,-3,*gmesh); // right
  TPZGeoElBC gbc3(elvec[1],5,-3,*gmesh); // right
  TPZGeoElBC gbc4(elvec[1],6,-3,*gmesh); // top
  TPZGeoElBC gbc5(elvec[2],5,-3,*gmesh); // top
  TPZGeoElBC gbc6(elvec[2],6,-3,*gmesh); // left
  TPZGeoElBC gbc7(elvec[3],5,-3,*gmesh); // left
  TPZGeoElBC gbc8(elvec[3],6,-3,*gmesh); // bottom
  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetForcingFunction(Forcing1);

  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;
  //Difusao pura. conveccao = 0.
  mat->SetParameters(gDif, 0., convdir);
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc;
  
  val2.Zero();
  bc = mat->CreateBC(-3, 0,val1,val2);

  
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bc);

#endif

  
  TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::SetCreateFunction(TPZCompElDisc::CreateDisc); 
  //template class TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>;
  
  cmesh->AutoBuild();
  //  cmesh->AdjustBoundaryElements();
  //  cmesh->CleanUpUnconnectedNodes();
  //  cmesh->ExpandSolution();
  
  return cmesh;
}


TPZCompMesh *CreateMesh2() {
  REAL co[9][2] = {{0.,0.},{0.,-1.},{1.,-1.},{1.,0.},{1.,1.},{0.,1.},{-1.,1.},{-1.,0.},{-1,-1}};
  int numquad = 0;
  int indicesquad[2][4] = {{0,1,2,3},{0,5,6,7}};
  int numtriang = 8;
  int indicestriang[8][3] = {{0,2,3},{0,1,2},{0,3,5},{3,4,5},{0,7,1},{7,8,1},{0,5,6},{0,6,7}};
  TPZGeoEl *elvec[4];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 9;
  int nelem = numquad + numtriang;
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el=0; el<numquad; el++) {
    TPZManVector<int,4> nodind(4);
    for(nod=0; nod<4; nod++) nodind[nod] = indicesquad[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index); 
  }
 
  for(el=0; el<numtriang; el++) {
    TPZManVector<int,3> nodind(3);
    for(nod=0; nod<3; nod++) nodind[nod]=indicestriang[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
  }


  gmesh->BuildConnectivity();

    
  for(int i = 0; i < nelem; i++){
    TPZVec<TPZGeoEl *> children, netos, bisnetos, tata1, tata2, tata3;
    //    cout << "\ngDivide[0] = \n" << gDivide[0];
    if(gDivide[0] == 1) {
      elvec[i]->Divide(children);
      //      cout <<  "\n Primeira divisao \n" ;
      if (gDivide[1] == 1) {
        for(int j = 0; j < children.NElements(); j++) {
          children[j]->Divide(netos);
	  //          cout <<  "\n Segunda divisao \n" ;
          if(gDivide[2] == 1) {
            for(int k = 0; k < netos.NElements(); k++) {
              netos[k]->Divide(bisnetos); 
	      //              cout <<  "\n Terceira divisao \n" ;
              if(gDivide[3] == 1) {
                for(int k = 0; k < bisnetos.NElements(); k++) {
                  bisnetos[k]->Divide(tata1); 
		  //                  cout <<  "\n Quarta divisao \n" ;
                  if(gDivide[4] == 1) {
                    for(int k = 0; k < tata1.NElements(); k++) {
                      tata1[k]->Divide(tata2); 
		      //                      cout <<  "\n Quinta divisao \n" ;
                      if(gDivide[5] == 1) {
                        for(int k = 0; k < tata2.NElements(); k++) tata2[k]->Divide(tata3); 
			//                        cout <<  "\n Sexta divisao \n" ;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }  
 
  TPZGeoElBC gbc;


#ifdef DESCONTINUIDADE
  TPZGeoElBC gbc1(elvec[0],4,-2,*gmesh); 
  TPZGeoElBC gbc2(elvec[1],4,-2,*gmesh); 
  TPZGeoElBC gbc3(elvec[6],4,-2,*gmesh); 
  TPZGeoElBC gbc4(elvec[7],4,-2,*gmesh); 
  TPZGeoElBC gbc5(elvec[3],3,-1,*gmesh); 
  TPZGeoElBC gbc6(elvec[3],4,-1,*gmesh); 
  TPZGeoElBC gbc7(elvec[5],3,-1,*gmesh); 
  TPZGeoElBC gbc8(elvec[5],4,-1,*gmesh); 

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);

  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;
  mat->SetParameters(gDif, 1., convdir);
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);

  TPZBndCond *bc[2];
  
  val2(0,0) = 1.;
  bc[0] = mat->CreateBC(-1,0,val1,val2);

  val2.Zero(); 
  bc[1] = mat->CreateBC(-2,0,val1,val2);

  cmesh->InsertMaterialObject(mat);
  int i;
  for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);
#endif

#ifdef DIFUSAO_EXP
  TPZGeoElBC gbc1(elvec[1],4,-3,*gmesh); 
  TPZGeoElBC gbc2(elvec[0],4,-3,*gmesh); 
  TPZGeoElBC gbc3(elvec[3],3,-3,*gmesh); 
  TPZGeoElBC gbc4(elvec[3],4,-3,*gmesh); 
  TPZGeoElBC gbc5(elvec[6],4,-3,*gmesh); 
  TPZGeoElBC gbc6(elvec[7],4,-3,*gmesh); 
  TPZGeoElBC gbc7(elvec[5],3,-3,*gmesh); 
  TPZGeoElBC gbc8(elvec[5],4,-3,*gmesh); 

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetForcingFunction(Forcing1);

  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;
  mat->SetParameters(gDif, 0., convdir);
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);

  TPZBndCond *bc;
  val2.Zero(); val1.Zero(); bc = mat->CreateBC(-3, 0, val1, val2);
  
  cmesh->InsertMaterialObject(mat);

  cmesh->InsertMaterialObject(bc);

#endif
  
  TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  //template class TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>;
  
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  
  return cmesh;
}

TPZCompMesh *CreateMeshPhil() {
  REAL co[9][2] = {{0.,0.},{0.,-1.},{1.,-1.},{1.,0.},{1.,1.},{0.,1.},{-1.,1.},{-1.,0.},{-1,-1}};
  int indices[4][4] = {{0,1,2,3},{0,3,4,5},{0,5,6,7},{0,7,8,1}};
  TPZGeoEl *elvec[8];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 9;
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

// triangular mesh
  int nelem = 8;
  int el;
  for(el=0; el<nelem/2; el++) {
    TPZManVector<int,4> nodind(3);
    int index;
    
    for(nod=0; nod<3; nod++) nodind[nod]=indices[el][(nod+1-el%2)%4];
    elvec[2*el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
    for(nod=0; nod<3; nod++) nodind[nod]=indices[el][(nod+3-el%2)%4];
    elvec[2*el+1] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
  }
  TPZGeoElBC gbc1(elvec[0],3,-2,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],4,-2,*gmesh); // right
  TPZGeoElBC gbc3(elvec[2],4,-1,*gmesh); // right
  TPZGeoElBC gbc4(elvec[3],3,-1,*gmesh); // top
  TPZGeoElBC gbc5(elvec[4],3,-2,*gmesh); // top
  TPZGeoElBC gbc6(elvec[4],4,-2,*gmesh); // left
  TPZGeoElBC gbc7(elvec[6],4,-1,*gmesh); // left
  TPZGeoElBC gbc8(elvec[7],3,-1,*gmesh); // bottom

// quadrilateral mesh  
//  int nelem = 4;
//  int el;
//  for(el=0; el<nelem; el++) {
//    TPZManVector<int,4> nodind(3);
//    for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
//    int index;
//    elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
//  }
//
//  TPZGeoElBC gbc1(elvec[0],5,-2,*gmesh); // bottom
//  TPZGeoElBC gbc2(elvec[0],6,-2,*gmesh); // right
//  TPZGeoElBC gbc3(elvec[1],5,-1,*gmesh); // right
//  TPZGeoElBC gbc4(elvec[1],6,-1,*gmesh); // top
//  TPZGeoElBC gbc5(elvec[2],5,-2,*gmesh); // top
//  TPZGeoElBC gbc6(elvec[2],6,-2,*gmesh); // left
//  TPZGeoElBC gbc7(elvec[3],5,-1,*gmesh); // left
//  TPZGeoElBC gbc8(elvec[3],6,-1,*gmesh); // bottom

  
  gmesh->BuildConnectivity();

    
  for(int i = 0; i < nelem; i++){
    TPZVec<TPZGeoEl *> children, netos, bisnetos, tata1, tata2, tata3;
    cout << "\ngDivide[0] = \n" << gDivide[0];
    if(gDivide[0] == 1) {
      elvec[i]->Divide(children);
      cout <<  "\n Primeira divisao \n" ;
      if (gDivide[1] == 1) {
        for(int j = 0; j < children.NElements(); j++) {
          children[j]->Divide(netos);
          cout <<  "\n Segunda divisao \n" ;
          if(gDivide[2] == 1) {
            for(int k = 0; k < netos.NElements(); k++) {
              netos[k]->Divide(bisnetos); 
              cout <<  "\n Terceira divisao \n" ;
              if(gDivide[3] == 1) {
                for(int k = 0; k < bisnetos.NElements(); k++) {
                  bisnetos[k]->Divide(tata1); 
                  cout <<  "\n Quarta divisao \n" ;
                  if(gDivide[4] == 1) {
                    for(int k = 0; k < tata1.NElements(); k++) {
                      tata1[k]->Divide(tata2); 
                      cout <<  "\n Quinta divisao \n" ;
                      if(gDivide[5] == 1) {
                        for(int k = 0; k < tata2.NElements(); k++) tata2[k]->Divide(tata3); 
                        cout <<  "\n Sexta divisao \n" ;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }  
 
  TPZGeoElBC gbc;

/* Lineare em y:
  TPZGeoElBC gbc1(elvec[0],5,-1,*gmesh); // bottom
  TPZGeoElBC gbc4(elvec[1],6,-2,*gmesh); // top
  TPZGeoElBC gbc5(elvec[2],5,-2,*gmesh); // top
  TPZGeoElBC gbc8(elvec[3],6,-1,*gmesh); // bottom
*/

  // bc -1 -> Dirichlet homogeneo
  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);


  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;
  mat->SetParameters(0.,1.,convdir);
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc[2]; 
    
  val2(0,0) = 1.;
  bc[0] = mat->CreateBC(-1,0,val1,val2);

  val2.Zero();
  bc[1] = mat->CreateBC(-2,0,val1,val2);
  
  cmesh->InsertMaterialObject(mat);
  int i;
  for(i=0; i<2; i++) 
    cmesh->InsertMaterialObject(bc[i]);

  
  TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  //template class TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>;
  
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  //  cmesh->ExpandSolution();
  
  return cmesh;
}


TPZCompMesh *CreateMesh3() {
  REAL co[9][2] = {{0.,0.},{0.,-1.},{1.,-1.},{1.,0.},{1.,1.},{0.,1.},{-1.,1.},{-1.,0.},{-1,-1}};
  int numquad = 0;
  int indicesquad[2][4] = {{0,1,2,3},{0,5,6,7}};
  int numtriang = 8;
  int indicestriang[8][3] = {{1,2,3},{0,1,3},{0,3,4},{0,4,5},{0,5,7},{5,6,7},{0,7,8},{0,8,1}};
  TPZGeoEl *elvec[4];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 9;
  int nelem = numquad + numtriang;
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el=0; el<numquad; el++) {
    TPZManVector<int,4> nodind(4);
    for(nod=0; nod<4; nod++) nodind[nod] = indicesquad[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index); 
  }
 
  for(el=0; el<numtriang; el++) {
    TPZManVector<int,3> nodind(3);
    for(nod=0; nod<3; nod++) nodind[nod]=indicestriang[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
  }


  gmesh->BuildConnectivity();

    
  for(int i = 0; i < nelem; i++){
    TPZVec<TPZGeoEl *> children, netos, bisnetos, tata1, tata2, tata3;
    //    cout << "\ngDivide[0] = \n" << gDivide[0];
    if(gDivide[0] == 1) {
      elvec[i]->Divide(children);
      //      cout <<  "\n Primeira divisao \n" ;
      if (gDivide[1] == 1) {
        for(int j = 0; j < children.NElements(); j++) {
          children[j]->Divide(netos);
	  //          cout <<  "\n Segunda divisao \n" ;
          if(gDivide[2] == 1) {
            for(int k = 0; k < netos.NElements(); k++) {
              netos[k]->Divide(bisnetos); 
	      //              cout <<  "\n Terceira divisao \n" ;
              if(gDivide[3] == 1) {
                for(int k = 0; k < bisnetos.NElements(); k++) {
                  bisnetos[k]->Divide(tata1); 
		  //                  cout <<  "\n Quarta divisao \n" ;
                  if(gDivide[4] == 1) {
                    for(int k = 0; k < tata1.NElements(); k++) {
                      tata1[k]->Divide(tata2); 
		      //                      cout <<  "\n Quinta divisao \n" ;
                      if(gDivide[5] == 1) {
                        for(int k = 0; k < tata2.NElements(); k++) tata2[k]->Divide(tata3); 
			//                        cout <<  "\n Sexta divisao \n" ;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }  
 
  TPZGeoElBC gbc;


#ifdef DESCONTINUIDADE
  TPZGeoElBC gbc1(elvec[0],3,-2,*gmesh); 
  TPZGeoElBC gbc2(elvec[0],4,-2,*gmesh); 
  TPZGeoElBC gbc3(elvec[2],4,-1,*gmesh); 
  TPZGeoElBC gbc4(elvec[3],4,-1,*gmesh); 
  TPZGeoElBC gbc5(elvec[5],3,-2,*gmesh); 
  TPZGeoElBC gbc6(elvec[5],4,-2,*gmesh); 
  TPZGeoElBC gbc7(elvec[6],4,-1,*gmesh); 
  TPZGeoElBC gbc8(elvec[7],4,-1,*gmesh); 

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);

  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;
  mat->SetParameters(gDif, 1., convdir);
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);

  TPZBndCond *bc[2];
  
  val2(0,0) = 1.;
  bc[0] = mat->CreateBC(-1,0,val1,val2);

  val2.Zero(); 
  bc[1] = mat->CreateBC(-2,0,val1,val2);

  cmesh->InsertMaterialObject(mat);
  int i;
  for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);
#endif

#ifdef DIFUSAO_EXP
  TPZGeoElBC gbc1(elvec[0],3,-3,*gmesh); 
  TPZGeoElBC gbc2(elvec[0],4,-3,*gmesh); 
  TPZGeoElBC gbc3(elvec[2],4,-3,*gmesh); 
  TPZGeoElBC gbc4(elvec[3],4,-3,*gmesh); 
  TPZGeoElBC gbc5(elvec[5],3,-3,*gmesh); 
  TPZGeoElBC gbc6(elvec[5],4,-3,*gmesh); 
  TPZGeoElBC gbc7(elvec[6],4,-3,*gmesh); 
  TPZGeoElBC gbc8(elvec[7],4,-3,*gmesh); 

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetForcingFunction(Forcing1);

  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;
  mat->SetParameters(gDif, 0., convdir);
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);

  TPZBndCond *bc;
  val2.Zero(); val1.Zero(); bc = mat->CreateBC(-3, 0, val1, val2);
  
  cmesh->InsertMaterialObject(mat);

  cmesh->InsertMaterialObject(bc);

#endif
  
  TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  //template class TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>;

  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  
  return cmesh;
}
