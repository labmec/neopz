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
#include "pzstrmatrix.h"

#include <stdio.h>
#include <fstream>

using namespace std;


TPZNonLinearAnalysis::TPZNonLinearAnalysis() : TPZAnalysis() {
	Mesh()->Solution().Zero();
	fSolution.Zero();
}

TPZNonLinearAnalysis::TPZNonLinearAnalysis(TPZCompMesh *mesh,ostream &out) : TPZAnalysis(mesh,out) {
	Mesh()->Solution().Zero();
	fSolution.Zero();
}

TPZNonLinearAnalysis::~TPZNonLinearAnalysis(void){}

  ofstream alphafile("alpha.txt");
void TPZNonLinearAnalysis::LineSearch(TPZFMatrix &Wn, TPZFMatrix &DeltaW, TPZFMatrix &NextW, REAL tol, int niter){
  REAL error = 2.*tol+1.;
  TPZFMatrix ak, bk, lambdak, muk, Interval;
  REAL NormResLambda, NormResMu;
  ///ak = Wn + 0.1 * DeltaW
  ak = DeltaW;
  ak *= 0.1;
  ak += Wn;
  ///bk = Wn + DeltaW
  bk = Wn; bk += DeltaW;
  ///Interval = (bk-ak)
  Interval = bk; Interval -= ak;
  int iter = 0;
  while(error > tol && iter < niter){
    iter++;
    ///lambdak = ak + 0.382*(bk-ak)
    lambdak = Interval; lambdak *= 0.382; lambdak += ak;
    ///muk = ak + 0.618*(bk-ak)
    muk = Interval; muk *= 0.618; muk += ak;
    ///computing residuals
    this->LoadSolution(lambdak);
    this->Assemble();
//     TPZStructMatrix::Assemble(fRhs, *(this->fCompMesh));
    NormResLambda = Norm(fRhs);
    this->LoadSolution(muk);
    this->Assemble();
//     TPZStructMatrix::Assemble(fRhs, *(this->fCompMesh));
    NormResMu = Norm(fRhs);
    if (NormResLambda > NormResMu){
      ak = lambdak;
    }
    else{
      bk = muk;
    }
    ///error = Norm(bk-ak)
    Interval = bk; Interval -= ak; error = Norm(Interval);
  }///while
  
  NextW = ak;
  NextW += bk;
  NextW *= 0.5;
  
  
  ///debug: valor do alpha
  TPZFMatrix alpha;
  alpha = NextW;
  alpha -= Wn;
  REAL MeanAlpha = 0.;
  for(int i = 0; i < alpha.Rows(); i++){
    alpha(i,0) = alpha(i,0)/DeltaW(i,0);
    MeanAlpha += alpha(i,0)/alpha.Rows();
  }
  alphafile << MeanAlpha << "\t";
  alphafile.flush();

}///void

void TPZNonLinearAnalysis::IterativeProcess(ostream &out,REAL tol,int numiter, bool linesearch) {

   int iter = 0;
   REAL error = 1.e10;
   int numeq = fCompMesh->NEquations();
   //Mesh()->Solution().Zero();
   //fSolution->Zero();

   TPZFMatrix prevsol(fSolution);
   if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
   
	TPZVec<REAL> coefs(1,1.);
	TPZFMatrix range(numeq,1,0.01);
	CheckConvergence(*this,fSolution,range,coefs);	
//         CheckConvergence<TPZNonLinearAnalysis>(*this,*fSolution,range,coefs);	
	
   while(error > tol && iter < numiter) {

      fSolution.Redim(0,0);
      Run();

      if (linesearch){
        TPZFMatrix nextSol;
        REAL LineSearchTol = 1e-3 * Norm(fSolution);
        const int niter = 10;
        this->LineSearch(prevsol, fSolution, nextSol, LineSearchTol, niter);
        fSolution = nextSol;
      }
      else{
        fSolution += prevsol;
      }
      
      prevsol -= fSolution;
      REAL norm = Norm(prevsol);
      out << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << endl;

      prevsol = fSolution;
      TPZAnalysis::LoadSolution();
      if(norm < tol) {
         out << "\nTolerancia atingida na iteracao : " << (iter+1) << endl;
         out << "\n\nNorma da solucao |Delta(Un)|  : " << norm << endl << endl;

      } else
      if( (norm - error) > 1.e-9 ) {
         out << "\nDivergent Method\n";
      }
      error = norm;
	  iter++;
	  out.flush();
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
