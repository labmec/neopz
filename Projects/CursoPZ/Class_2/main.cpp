/**
 * @file
 * @brief Implements the use of the matrices and solvers as first tutorial example of the matrix and analysis NeoPZ modules
 */
#include "pzvec.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include <pzskylmat.h>
#include <pzstepsolver.h> 

#include <iostream>

using namespace std;

int main() {
  long neq=1000;
//  long banda=50;
 
  TPZFMatrix<REAL>  *cheia = new TPZFMatrix<REAL>;
  cheia->AutoFill(neq,neq,0);
//  FillMatrix(*cheia,neq,banda);


  TPZSkylMatrix<REAL>  *skyline = new TPZSkylMatrix<REAL>;
  skyline->AutoFill(neq,neq,1);
//  FillMatrix(*skyline,neq,banda);

//  cheia->Print( "Matriz cheia ",cout);
//  skyline->Print( "Matriz skyline ",cout);

  TPZFMatrix<REAL> F;
  F.AutoFill(neq,1,0);
//  FillF(F,neq,1);
  F.Print("Vetor de Carga ",cout);

  TPZFMatrix<REAL> resultcheia;
  /*  TPZStepSolver direct(&cheia);
  direct.SetDirect(ECholesky);
  direct.Solve(F,resultcheia);*/

  TPZStepSolver<REAL> step(cheia);
  TPZStepSolver<REAL> precond(step);
  long numiterpre =2;
  long numiter = 5;
  double overrelax = 1.1;
  double tol = 1e-8;
  precond.SetSSOR(numiterpre,overrelax,tol,0);
  step.SetCG(numiter,precond,tol,0);
  step.Solve(F,resultcheia);

  
  resultcheia.Print("Solucao ",cout);
	return 0;
}

