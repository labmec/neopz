/**
 * @file
 * @brief Contains implementations of the TPZNonLinearAnalysis methods.
 */

#include "pznonlinanalysis.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzfmatrix.h"
#include "TPZMatrixSolver.h"
#include "TPZMaterial.h"
#include "pzelmat.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "checkconv.h"

#include "pzlog.h"

#include <stdio.h>
#include <fstream>

#ifdef PZ_LOG
static TPZLogger logger("pz.nonlinearanalysis");
#endif

using namespace std;

TPZNonLinearAnalysis::TPZNonLinearAnalysis() : TPZLinearAnalysis() {
	if(Mesh()) Mesh()->Solution().Zero();
	fSolution.Zero();
}

TPZNonLinearAnalysis::TPZNonLinearAnalysis(TPZCompMesh *mesh,std::ostream &out) : TPZLinearAnalysis(mesh,true,out) {
	if(Mesh()) Mesh()->Solution().Zero();
	fSolution.Zero();
}

TPZNonLinearAnalysis::~TPZNonLinearAnalysis() {
}


//#define DEBUGLINESEARCH
#ifdef PZDEBUGLINESEARCH
ofstream alphafile("c:\\Temp\\tmp\\alpha.txt");
#endif
REAL TPZNonLinearAnalysis::LineSearch(const TPZFMatrix<STATE> &Wn, TPZFMatrix<STATE> DeltaW, TPZFMatrix<STATE> &NextW, REAL tol, int niter){
	REAL error = 2.*tol+1.;
	REAL A = 0.1, B = 2., L = 0, M = 0.;
	TPZFMatrix<STATE> ak, bk, lambdak, muk, Interval;
	REAL NormResLambda = 0., NormResMu = 0.;
	//ak = Wn + 0.1 * DeltaW
	ak = DeltaW;
	ak *= A;
	ak += Wn;
	//bk = Wn + 2. DeltaW
	bk = DeltaW;
	bk *= B;
	bk += Wn;
	//Interval = (bk-ak)
	Interval = bk; Interval -= ak;
	int iter = 0;
	int KeptVal = -1; //0 means I have residual(labmda); 1 means I have residual(mu); -1 means I have nothing
	while(error > tol && iter < niter){
		iter++;
		
		if (KeptVal != 0){
			L = 0.382*(B-A)+A;
			//lambdak = ak + 0.382*(bk-ak)
			lambdak = Interval; lambdak *= 0.382; lambdak += ak;
			//computing residual
			this->LoadSolution(lambdak);
#ifdef PZ_LOG
			LOGPZ_DEBUG(logger,"After LoadSolution")
#endif
			//		LogWellSolution(*this->Mesh(), 6);
			this->AssembleResidual();
#ifdef PZ_LOG
			LOGPZ_DEBUG(logger,"After AssembleResidual")
#endif
			//		LogWellSolution(*this->Mesh(), 6);
			NormResLambda = Norm(fRhs);
		}
		
		if (KeptVal != 1){
			//muk = ak + 0.618*(bk-ak)
			M = 0.618*(B-A)+A;
			muk = Interval; muk *= 0.618; muk += ak;
			this->LoadSolution(muk);
			this->AssembleResidual();
			NormResMu = Norm(fRhs);
		}
		
		if (NormResLambda > NormResMu){
			A = L;
			L = M;
			ak = lambdak;
			lambdak = muk;
			NormResLambda = NormResMu;
			KeptVal = 0;
		}
		else{
			B = M;
			M = L;
			bk = muk;
			muk = lambdak;
			NormResMu = NormResLambda;
			KeptVal = 1;
		}
		//error = Norm(bk-ak)
		Interval = bk; Interval -= ak; error = Norm(Interval);
		
		//alpha shall be alpha <= 1
		if(A > 1. && B > 1.) break;
		
	}//while
	
	double ALPHA = 0.5*(A + B);
	NextW = ak;
	NextW += bk;
	NextW *= 0.5;
	
	
#ifdef PZDEBUGLINESEARCH
	//debug: valor do alpha
	TPZFMatrix<REAL> alpha;
	alpha = NextW;
	alpha -= Wn;
	REAL sum = 0.;
	int ncontrib = 0;
	for(int i = 0; i < alpha.Rows(); i++){
		if (DeltaW(i,0)){
			alpha(i,0) = alpha(i,0)/DeltaW(i,0);
			sum += alpha(i,0);
			ncontrib++;
		}
	}
	//REAL MeanAlpha = sum/ncontrib;
	alphafile << /*MeanAlpha << "\t" <<*/ "ALPHA = " << ALPHA << "\n";
	alphafile.flush();
#endif
	
	if(ALPHA > 1.){ //alpha shall be alpha <= 1
		NextW = Wn;
		NextW += DeltaW;
#ifdef PZDEBUGLINESEARCH
		alphafile << "ALPHA LIMIT APPLIED. Alpha = 1.\n";
#endif
		return 1.;
	}
	
	return ALPHA;
	
}//void

void TPZNonLinearAnalysis::IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv) {
	
	int iter = 0;
	REAL error = 1.e10;
	int numeq = fCompMesh->NEquations();
	
	TPZFMatrix<STATE> prevsol(fSolution);
	if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
	
	if(checkconv){
		TPZVec<REAL> coefs(1,1.);
		TPZFMatrix<STATE> range(numeq,1,1.);
		CheckConvergence(*this,fSolution,range,coefs);
	}
	
	while(error > tol && iter < numiter) {
		
//        fSolution.Redim(0,0);
		Assemble();
		Solve();
		if (linesearch){
			TPZFMatrix<STATE> nextSol;
			REAL LineSearchTol = 1e-3 * Norm(fSolution);
			const int niter = 10;
			this->LineSearch(prevsol, fSolution, nextSol, LineSearchTol, niter);
			fSolution = nextSol;
		}
		else{
			TPZFMatrix<STATE> sol = fSolution;
			sol += prevsol;
		}
		
		prevsol -= fSolution;
		REAL normDeltaSol = Norm(prevsol);
		prevsol = fSolution;
		this->LoadSolution(fSolution);
		this->AssembleResidual();
		double NormResLambda = Norm(fRhs);
		double norm = NormResLambda;
		out << "Iteracao n : " << (iter+1) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << NormResLambda << endl;
		
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

/** @brief Zeroes entries of val vector and deriv matrix. */
void NullForce(TPZVec<REAL> &/*point*/,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv) {
    int i,cap = val.NElements() ;
    deriv.Zero();
    for(i=0;i<cap;i++) val[i] = 0.;
}
REAL TPZNonLinearAnalysis::SolutionNorm(){
	return Norm(fSolution);
}

void TPZNonLinearAnalysis::ComputeTangent(TPZFMatrix<STATE> &tangent, TPZVec<REAL> &coefs, int icase){
	
	int neq = fCompMesh->NEquations();
	tangent.Redim(neq,neq);
	TPZFMatrix<STATE> rhs(neq,1);
	fStructMatrix->Assemble(tangent,rhs,NULL);
}

int TPZNonLinearAnalysis::NumCases() {
	return 1;
}

void TPZNonLinearAnalysis::Residual(TPZFMatrix<STATE> &residual, int icase){
	int neq = fCompMesh->NEquations();
	TPZFMatrix<STATE> tangent(neq,neq);
	residual.Redim(neq,1);
	fStructMatrix->Assemble(tangent,residual,NULL);
	residual *= -1.;
}

void TPZNonLinearAnalysis::LoadSolution(const TPZFMatrix<STATE> &state){
	fSolution = state;
	LoadSolution();
}

void TPZNonLinearAnalysis::LoadSolution(){
    TPZLinearAnalysis::LoadSolution();
}

void TPZNonLinearAnalysis::LoadState(TPZFMatrix<STATE> &state){
	this->LoadSolution(state);
}
