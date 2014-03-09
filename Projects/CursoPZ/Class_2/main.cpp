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
  long banda=50;
  long i;
 
  TPZFMatrix<REAL>  *cheia = new TPZFMatrix<REAL>(neq,neq, 0.);
  cheia->AutoFill();
//  FillMatrix(*cheia,neq,banda);

  TPZVec <long> skyvec(neq,0);
  for (i=0;i<neq;i++){
    skyvec[i] = i-banda;
    if(skyvec[i] < 0) skyvec[i] = 0;
  }

  TPZSkylMatrix<REAL>  *skyline = new TPZSkylMatrix<REAL>(neq,skyvec);
  skyline->AutoFill();
//  FillMatrix(*skyline,neq,banda);

  cheia->Print( "Matriz cheia ",cout);
  skyline->Print( "Matriz skyline ",cout);

  TPZFMatrix<REAL> F(neq,1);
  F.AutoFill();
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

