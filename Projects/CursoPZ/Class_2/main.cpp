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

void FillMatrix(TPZMatrix &mat,int neq, int banda);
void FillF(TPZFMatrix &f, int neq, int nst);

int main(){
  int neq=1000;
  int banda=50;
  int i;
 
  TPZFMatrix  *cheia = new TPZFMatrix(neq,neq, 0.);
  FillMatrix(*cheia,neq,banda);

  TPZVec <int> skyvec(neq,0);
  for (i=0;i<neq;i++){
    skyvec[i] = i-banda;
    if(skyvec[i] < 0) skyvec[i] = 0;
  }

  TPZSkylMatrix  *skyline = new TPZSkylMatrix(neq,skyvec);
  FillMatrix(*skyline,neq,banda);

  cheia->Print( "Matriz cheia ",cout);
  skyline->Print( "Matriz skyline ",cout);

  TPZFMatrix F;
  FillF(F,neq,1);
  F.Print("Vetor de Carga ",cout);

  TPZFMatrix resultcheia;
  /*  TPZStepSolver direct(&cheia);
  direct.SetDirect(ECholesky);
  direct.Solve(F,resultcheia);*/

  TPZStepSolver step(cheia);
  TPZStepSolver precond(step);
  int numiterpre =2;
  int numiter = 5;
  double overrelax = 1.1;
  double tol = 1e-8;
  precond.SetSSOR(numiterpre,overrelax,tol,0);
  step.SetCG(numiter,precond,tol,0);
  step.Solve(F,resultcheia);

  
  resultcheia.Print("Solucao ",cout);
	return 0;
}

void FillMatrix(TPZMatrix &mat,int neq, int banda){
  int i,j;
  for ( i=0;i<neq;i++){
    mat(i,i)=2000.;
    int col = i+banda+1;
    if (col >= neq ) col = neq; 
    for ( j=i+1;j<col;j++){
      mat(i,j)=500./(j-i);
      mat(j,i)=500./(j-i);
    }
  }
}

void FillF(TPZFMatrix &f,int neq, int nst){
  double PI = 3.1416;
  f.Resize(neq,nst);
  int i,j;
  for (i=0;i<neq;i++){
    for (j=0;j<nst;j++){
      double val = (double)(i%10);
      f(i,j)=PI*val;
    }
  }
}
