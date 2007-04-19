#include "pznonlinanalysis.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzfmatrix.h"
#include "pztempmat.h"
#include "pzsolve.h"
#include "pzmaterial.h"
#include "pzelmat.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "checkconv.h" 
#include <fstream>
using namespace std;
#include <string.h>
#include <stdio.h>
#include "pzstrmatrix.h"

TPZNonLinearAnalysis::TPZNonLinearAnalysis() : TPZAnalysis() {
	Mesh()->Solution().Zero();
	fSolution.Zero();
}

TPZNonLinearAnalysis::TPZNonLinearAnalysis(TPZCompMesh *mesh,ostream &out) : TPZAnalysis(mesh,out) {
	Mesh()->Solution().Zero();
	fSolution.Zero();
}

TPZNonLinearAnalysis::~TPZNonLinearAnalysis(void){}

void TPZNonLinearAnalysis::IterativeProcess(ostream &out,REAL tol,int numiter) {

   int iter = 0;
   REAL error = 1.e10;
   int numeq = fCompMesh->NEquations();
   //Mesh()->Solution().Zero();
   //fSolution->Zero();

   TPZFMatrix prevsol(fSolution);
   if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
//	TPZVec<REAL> coefs(1,1.);
//	TPZFMatrix range(numeq,1,0.01);
//    CheckConvergence<TPZNonLinearAnalysis>(*this,*fSolution,range,coefs);	
	
   while(error > tol && iter < numiter) {

      fSolution.Redim(0,0);
      Run();
      REAL norm = SolutionNorm();
      cout << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << endl;
      out << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << endl;

      fSolution += prevsol;
	  prevsol = fSolution;
	  TPZAnalysis::LoadSolution();
      if(norm < tol) {
         cout << "\nTolerancia atingida na iteracao : " << (iter+1) << endl;
         cout << "\n\nNorma da solucao |Delta(Un)|  : " << norm << endl << endl;
         out << "\nTolerancia atingida na iteracao : " << (iter+1) << endl;
         out << "\n\nNorma da solucao |Delta(Un)|  : " << norm << endl << endl;

      } else
      if( (norm - error) > 1.e-9 ) {
         cout << "\nDivergent Method\n";
         out << "\nDivergent Method\n";
      }
      error = norm;
	  iter++;
   }
}

void NullForce(TPZVec<REAL> &point,TPZVec<REAL>&val,TPZFMatrix &deriv);
REAL TPZNonLinearAnalysis::SolutionNorm(){

//   REAL trueerr,L2,estimate;
//   fCompMesh->EvaluateError(NullForce,trueerr,L2,estimate);
	return Norm(fSolution);
}

void NullForce(TPZVec<REAL> &/*point*/,TPZVec<REAL> &val,TPZFMatrix &deriv) {

    int i,cap = val.NElements() ;
    deriv.Zero();
    for(i=0;i<cap;i++) val[i] = 0.;
}

void TPZNonLinearAnalysis::ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &coefs, int icase){

	int neq = fCompMesh->NEquations();
	tangent.Redim(neq,neq);
	TPZFMatrix rhs(neq,1);
	TPZStructMatrix::Assemble(tangent, rhs, *Mesh());
}

int TPZNonLinearAnalysis::NumCases(){
	return 1;
}

void TPZNonLinearAnalysis::Residual(TPZFMatrix &residual, int icase){
	int neq = fCompMesh->NEquations();
	TPZFMatrix tangent(neq,neq);
	residual.Redim(neq,1);
	TPZStructMatrix::Assemble(tangent, residual, *Mesh());
	residual *= -1.;
}

void TPZNonLinearAnalysis::LoadSolution(TPZFMatrix &state){

	Mesh()->LoadSolution(state);
}

void TPZNonLinearAnalysis::LoadState(TPZFMatrix &state){

	Mesh()->LoadSolution(state);
}
