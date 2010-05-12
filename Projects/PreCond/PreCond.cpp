//$Id: PreCond.cpp,v 1.6 2010-05-12 18:22:49 phil Exp $

/**
 * This project test some preconditioners for 3D elastic problems.
 * Aug 24, 2005
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzlog.h"

#include "pzvec.h"

#include "pzcmesh.h"

#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"

#include "TPZCopySolve.h"
#include "TPZStackEqnStorage.h"

#include "pzsbstrmatrix.h"
#include "pzstepsolver.h"

#include "pzadmchunk.h"

#include "pzbndcond.h"
#include "pzelast3d.h"
#include "pzblockdiag.h"
#include "pzvisualmatrix.h"
#include "cg.h"

#include <time.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "gmres.h"

#include "TPZTimer.h"

using namespace std;

TPZCompMesh * CreateSimpleMesh(int h);
TPZGeoMesh * CreateSimpleGeometry();
void PostProcessing(TPZAnalysis &an, char * filedx);
void VisualMatrix(TPZCompMesh &cmesh, int resolution, char *filedx);
void ReadMesh();


int main2(){
  const int neq = 5000;
  TPZFMatrix matrix(neq,neq,1.), F(neq, 1, 1.), result(neq, 1, 0.), residual(neq, 1, 0.), scratch(neq,1,0.);
  //result(0,0) = 9999999.;
  for(int i = 0; i < neq; i++) matrix(i,i) = 5;//1000.*(i+1);
  TPZCopySolve precond( &matrix );
  int niter = 1000;
  REAL tol = 1.e-10;
//  matrix.SolveCGown(niter, precond, F, result, &residual, tol, 1);
  matrix.SolveCG(niter, precond, F, result, &residual, tol, 1);
//  matrix.SolveDirect(F, ECholesky);
//  matrix.SolveJacobi(niter, F, result, &residual, scratch, tol, 1);
  TPZCounter after(TPZFlopCounter::gCount);

  std::cout << "\n\nSolucao\n";
  result.Print("Solucao", std::cout);
  std::cout << "\n\nresiduo\n";
  residual.Print("Residuo", std::cout);

  std::cout << "\n\nF\n";
  F.Print("F", std::cout);
  return 0;
}
#include "pzshapelinear.h"
using namespace pzshape;

int main(){
  InitializePZLOG("log4cxx.cfg");
TPZMaterial::gBigNumber= 1.e6;
//TPZShapeLinear::fOrthogonal = TPZShapeLinear::Legendre;
//TPZCompEl::fOrthogonal = TPZCompEl::Chebyshev;
//TPZShapeLinear::fOrthogonal = TPZCompEl::Chebyshev;

//  TPZBlockDiagonal::main();
  int nmaxeq = 500000;//numero maximo de equacoes do problema
  int h, p;
  for(p = 4; p < 5; p++){
    for(h =4; h < 5; h++){
      TPZCompEl::SetgOrder(p);

      TPZCompMesh *cmesh;
      cmesh = CreateSimpleMesh(h);
      TPZGeoMesh *gmesh = cmesh->Reference();

      TPZAnalysis an(cmesh);

      cout << "Analysis created" << endl;
      //Escolha do metodo de resolucao do sistema linear

      //metodo direto: foi escolhido o metodo frontal
//#define direct
#ifdef direct

      // Melhor caso para GEM e elementos finitos. Melhor metodo direto
      /*TPZParFrontStructMatrix <TPZFrontNonSym>*/ 
		TPZSkylineStructMatrix Matrix(cmesh);

      an.SetStructuralMatrix(Matrix);
      TPZStepSolver step;
      step.SetDirect( ECholesky );
      an.SetSolver(step);
#endif




#define iter
#ifdef iter
      //TPZFStructMatrix Matrix(cmesh); //Matriz cheia
      TPZSpStructMatrix Matrix(cmesh); //Matriz esparsa
      //TPZSkylineStructMatrix Matrix(cmesh); //Matriz skyline SimetricaTimesBetaPlusZ
      //TPZSBandStructMatrix Matrix(cmesh); //Matriz Banda Simetrica

      cout << "before SetStructuralMatrix\n";
      an.SetStructuralMatrix(Matrix);
      cout << "before step constructor\n";
      TPZStepSolver step( Matrix.Create() );
      an.SetSolver(step);

      cout << "matrix created " << step.Matrix()->Rows() <<  endl;
      REAL tol = 1.e-10;

      //Com precondicionador
      //TPZMatrixSolver * precond = an.BuildPreconditioner(TPZAnalysis::EBlockJacobi , true);
      //TPZMatrixSolver * precond = an.BuildPreconditioner(TPZAnalysis::EElement , true);
//      TPZMatrixSolver * precond = an.BuildPreconditioner(TPZAnalysis::ENodeCentered , false);
      // Sem pre-condicionador
      //TPZCopySolve * precond = new TPZCopySolve( Matrix.Create() );  step.ShareMatrix( *precond );
      //Precondicionador jacobi
		TPZStepSolver * precond = new TPZStepSolver( Matrix.Create() ); step.ShareMatrix( *precond ); precond->SetJacobi(1, 0.0, 0);
      //TPZStepSolver * precond = new TPZStepSolver( Matrix.Create() ); step.ShareMatrix( *precond ); precond->SetSOR(1, 0.0, 0);
      //TPZStepSolver * precond = new TPZStepSolver( Matrix.Create() ); step.ShareMatrix( *precond ); precond->SetSSOR(1,1., 0.0, 0);
      TPZStepSolver jac;
      jac.SetSSOR(1,1.1,0.,0);
//      jac.SetJacobi(1,0.,0);
      jac.ShareMatrix(step);

      //step.SetGMRES( 2000, 2000, jac, tol, 0 );
      //step.SetGMRES(  an.Mesh()->NEquations() + 2,  an.Mesh()->NEquations() + 1, *precond, tol, 0 );
      //step.CGown(&A, &x, &b, &max_iter, Real &tol)
      step.SetCG(2000, *precond, tol, 0);
      //step.SetCG(100, jac, tol, 0);
      //step.SetCG( an.Mesh()->NEquations() + 2, *precond, tol,0);

      an.SetSolver(step);

      delete precond;
#endif

      //Isso nao serve pra nada. Eh so pra evitar de tentar resolver problemas muito grandes.
      //nmaxeq esta definido no comeco do codigo
      cout << "\nNumero de equacoes: " << an.Mesh()->NEquations() << endl;
      if (an.Mesh()->NEquations() > nmaxeq) {
	cout << "skipping simulation...\n" << endl;
	delete cmesh;
	delete gmesh;
	break;
      }

      //Resolve o problema
      cout << "p = " << p << " h = " << h << endl;
		an.Assemble();
		// analisar o tempo para esta chamada
		an.Solve();

#define geradxTPZShapeLinear
#ifdef geradx

      //Gera arquivo de soluicao .dx
      char filedx[20];
      sprintf(filedx,"test_p%d_h%d.dx",p,h);
      //      PostProcessing(an, filedx);

      //Gera arquivo da densidade da matriz (+- Spy do MatLab)
      char filevisualdx[20];
      sprintf(filevisualdx,"Visual_p%d_h%d.dx",p,h);
      VisualMatrix(*cmesh, 100, filevisualdx);

 #endif

      delete cmesh;
      delete gmesh;
      cout << "\nacabei" << endl;
    }
  }

  return 0;
}

void PostProcessing(TPZAnalysis &an, char* filedx){
  /**** Aqui faz DX ****/
	TPZVec<std::string> scalnames(1);
	TPZVec<std::string> vecnames(1);
  scalnames[0] = "Solution";
  vecnames[0] = "Derivate";
  an.DefineGraphMesh(2,scalnames,vecnames,filedx);
  an.PostProcess(2);
}

void VisualMatrix(TPZCompMesh &cmesh, int resolution, char *filedx){
  TPZFMatrix fillin;
  cmesh.ComputeFillIn(resolution,fillin);
  VisualMatrix(fillin , filedx);
}

#include "pzmatorthotropic.h"
#include "pzpoisson3d.h"
TPZCompMesh * CreateSimpleMesh(int h){
  TPZGeoMesh * gmesh = CreateSimpleGeometry();

  TPZVec<TPZGeoEl*> filhos;
  for(int i = 0; i < h; i++){
    int n = gmesh->NElements();
    for(int j = 0; j < n; j++){
      if (gmesh->ElementVec()[j]->Dimension() == 3) gmesh->ElementVec()[j]->Divide(filhos);
    }
  }

  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);

      TPZFMatrix axes(3,3,0.0); axes(0,0) = 1.; axes(1,1) = 1., axes(2,2) = 1.;
      REAL E = 205000.;
      REAL poisson = 0.3;
//      REAL G = 79000.;
//       TPZMatOrthotropic* mat = new TPZMatOrthotropic(1, axes, E, E, E, poisson, poisson, poisson, G, G, G);
//       TPZFMatrix bodyforces(3,1,0.0);
//       bodyforces(2,0) = 0.08;
//       mat->SetMaterial(bodyforces);


  TPZVec<REAL> force(3,1.);
  TPZElasticity3D *mat = new TPZElasticity3D(1,E,poisson,force);
	TPZAutoPointer<TPZMaterial> matauto(mat);
//        TPZMatPoisson3d *mat = new TPZMatPoisson3d(1, 3);
//        TPZVec<REAL> dir(3, 0.);
//        mat->SetParameters(1., 0., dir);
//        mat->SetInternalFlux(1.0);

  int nstate = mat->NStateVariables();
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond * bc;
  bc = mat->CreateBC(matauto,-1,0,val1,val2);
  cmesh->InsertMaterialObject(matauto);
  cmesh->InsertMaterialObject(bc);
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  return cmesh;
}

TPZGeoMesh * CreateSimpleGeometry(){

//#define lermalha
//#ifdef lermalha

//  void ReadMesh();

//#endif

  REAL co[8][3] = {{0.,0.,0.},{0.,0.,0.1},{0.,0.1,0.1},{0.,0.1,0.},{1.,0.,0.},{1.,0.,0.1},{1.,0.1,0.1},{1.,0.1,0.}};
  int indices[1][8] = {{0,1,2,3,4,5,6,7}};
  TPZGeoEl *elvec[1];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 8;
  int nelem = 1;
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,3> coord(3);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    coord[2] = co[nod][2];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el, index;
  for(el=0; el<nelem; el++) {
    TPZManVector<int,8> nodind(8);
    for(nod=0; nod<8; nod++) nodind[nod]=indices[el][nod];
    elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
  }

  gmesh->BuildConnectivity();

  //Definicao de el contorno via lado 3D
  TPZGeoElBC gbc1(elvec[0],20,-1,*gmesh);

  //Definicao de el contorno via noh 2D
//   TPZManVector<int,4> nodind(4);
//   nodind[0] = 0;
//   nodind[1] = 1;
//   nodind[2] = 2;
//   nodind[3] = 3;
//   gmesh->CreateGeoElement(EQuadrilateral,nodind,-1,index);
  return gmesh;
}

// void ReadMesh(){
//
//   int dado1, dado2, dado3;
//   char texto[256];
//   ifstream entrada("break_d.1.node");
//   std::string str, chave;
//   entrada >> str;
//   //entrada >> nnod;
//   cout << "str " << str << endl;
//   //cout << "nnod " << nnod << endl;
//   std::string str2(str,0,8);
//   cout << "str2 " << str2 << endl;
//
//   chave = "Vertices";
//   cout << "chave " << chave << endl;
//   while (str2 != chave){
//     entrada.getline(texto,sizeof(texto));
//     entrada >> str2;
// //    cout << "str " << str << endl;
//     cout << "str2 " << str2 << endl;
//   }
//   entrada.getline(texto,sizeof(texto));
//   entrada >> dado1 >> dado2 >> dado3;
//   cout << "\n" << dado1 << "\n" << dado2 << "\n" << dado3;
//   entrada.close();
//
//   ofstream saida("elipse_a.1.txt");
//   saida << "\n" << dado1 << "\t" << dado2 << "\t" << dado3;
//   saida.close();
//
// }

