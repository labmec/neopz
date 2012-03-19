/**
 * \file
 * @brief Contains implementations of the TPZBlackOilAnalysis methods.
 */
//$Id: pzblackoilanalysis.cpp,v 1.7 2011-04-01 16:23:11 phil Exp $

#include "pzblackoilanalysis.h"
#include "pzblackoil2p3d.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzseqsolver.h"
#include "checkconv.h"

using namespace std;

#ifdef _AUTODIFF

TPZBlackOilAnalysis::TPZBlackOilAnalysis(TPZCompMesh *mesh, double TimeStep, std::ostream &out):TPZNonLinearAnalysis(mesh,out){
	this->fTimeStep = TimeStep;
	this->fSimulationTime = 0.;
	this->SetConvergence(0, 0.);
	this->SetNewtonConvergence(0, 0.);
	this->SetInitialSolutionAsZero();
	this->SetSaveFrequency(0,0);
}

TPZBlackOilAnalysis::~TPZBlackOilAnalysis(){
	//nothing to be done here
}

void TPZBlackOilAnalysis::SetInitialSolution(TPZFMatrix<REAL> & InitialSol){
	const int nrows = this->Mesh()->Solution().Rows();
	const int ncols = this->Mesh()->Solution().Cols();
	if ( (InitialSol.Rows() != nrows) || (InitialSol.Cols() != ncols) ){
		PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << std::endl;
	}
	else{
		this->fSolution = InitialSol;
	}
	TPZAnalysis::LoadSolution();
}

void TPZBlackOilAnalysis::SetInitialSolutionAsZero(){
	TPZFMatrix<REAL> & MeshSol = this->Mesh()->Solution();
	this->fSolution.Redim( MeshSol.Rows(), MeshSol.Cols() );
	this->fSolution.Zero();
}

void TPZBlackOilAnalysis::AssembleResidual(){
	this->SetCurrentState();
	int sz = this->Mesh()->NEquations();
	this->Rhs().Redim(sz,1);
// #warning FIX ME!!
	//  TPZStructMatrix::Assemble(this->Rhs(), *this->Mesh());
	this->fRhs += fLastState;
}//void

void TPZBlackOilAnalysis::Run(std::ostream &out, bool linesearch){
	
	ofstream File("Pocos.txt");
	File << "t(mes)\tPi(Pa)\tPp(Pa)\tQwiSC(m3/day)\tQoiSC(m3/day)\tQwpSC(m3/day)\tQopSC(m3/day)\t\tQwiFundo(m3/day)\tQoiFundo(m3/day)\tQwpFundo(m3/day)\tQopFundo(m3/day)\n";
	
	this->fSimulationTime = 0.;
	this->PostProcess(this->fDXResolution);
	
	double nextDeltaT = this->fTimeStep;
	this->SetAllMaterialsDeltaT();
	const double TotalTime = fTimeStep*fNIter;
	
	TPZFMatrix<REAL> prevsol, lastsol;
	
	for(; this->fSimulationTime < TotalTime; ){
		
		this->fTimeStep = nextDeltaT;
		this->SetAllMaterialsDeltaT();
		
		//Computing residual of last state solution
		this->SetLastState();
		this->Assemble();
		fLastState = this->fRhs;
		prevsol = fSolution;
		lastsol = fSolution;
		//Newton's method
		this->SetCurrentState();
		REAL error = this->fNewtonTol * 2. + 1.;
		int iter = 0;
		while(error > this->fNewtonTol && iter < this->fNewtonMaxIter){
			
			fSolution.Redim(0,0);
			this->Assemble();
			this->fRhs += fLastState;
			
			this->Solve();
			
			if (linesearch){
				TPZFMatrix<REAL> nextSol;
				REAL LineSearchTol = 1e-3 * Norm(fSolution);
				const int niter = 3;
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
			
			error = norm;
			iter++;
			
			if((iter % 20) == 0){
				//Computing residual of last state solution
				fSolution = prevsol;
				TPZAnalysis::LoadSolution();
				double fator = 0.1;//(sqrt(5.)-1.)/2.;
				this->TimeStep() *= fator;
				nextDeltaT = this->TimeStep();
				cout << "\nMultiplicando passo de tempo por " << fator << "\n";
				this->SetAllMaterialsDeltaT();
				this->SetLastState();
				this->Assemble();
				fLastState = this->fRhs;
				this->SetCurrentState();
				this->fNIter = fCurrentStep + (1./fator) * (fNIter - fCurrentStep) + 0.5;
			}
			
		}//Newton's iterations
		
		if(iter < 20){
			nextDeltaT = fTimeStep * 5./((double)(iter));
		}
		
		//DEBUG
		this->AssembleResidual();
		//DEBUG //
		
		prevsol = fSolution;
		this->fSimulationTime += this->TimeStep();
		
		if (this->fSaveFrequency){
			if (!(this->fCurrentStep % fSaveFrequency)){
				this->PostProcess(this->fDXResolution);
				double VazaoAguaIsc, VazaoAguaPsc, VazaoOleoIsc, VazaoOleoPsc;
				double VazaoAguaIFundo, VazaoAguaPFundo, VazaoOleoIFundo, VazaoOleoPFundo;
				Vazao(*this,-1,VazaoAguaIsc,VazaoOleoIsc,VazaoAguaIFundo,VazaoOleoIFundo);
				Vazao(*this,-2,VazaoAguaPsc,VazaoOleoPsc,VazaoAguaPFundo,VazaoOleoPFundo);
				File << this->fSimulationTime/2629743.8 << "\t" << PressaoMedia(*this,-1) << "\t" << PressaoMedia(*this,-2) << "\t"
				<< VazaoAguaIsc << "\t" 
				<< VazaoOleoIsc << "\t" 
				<< VazaoAguaPsc << "\t" 
				<< VazaoOleoPsc << "\t" << "\t"
				<< VazaoAguaIFundo << "\t" 
				<< VazaoOleoIFundo << "\t" 
				<< VazaoAguaPFundo << "\t" 
				<< VazaoOleoPFundo << "\n";
				File.flush();
			}
		}
		
		prevsol -= lastsol;
		REAL steadynorm = Norm(prevsol);
		std::cout << "*********** Steady state error at iteration " << this->fCurrentStep << ", " << this->fSimulationTime/2629743.8 << " meses " << " = " << steadynorm << "\n\n";
		if (!fForceAllSteps){
			if (steadynorm < this->fSteadyTol){
				std::cout << "Steady state solution achieved\n\n";
				this->fNIter = this->fCurrentStep;
				break;
			}
		}
		std::cout.flush();
		
	}//time step iterations
	
}//method


void TPZBlackOilAnalysis::SetLastState(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZBlackOil2P3D * blackoilmat = dynamic_cast< TPZBlackOil2P3D *>(matit->second.operator->());
		if (blackoilmat){
			blackoilmat->SetLastState();
		}
	}
}


void TPZBlackOilAnalysis::SetCurrentState(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZBlackOil2P3D * blackoilmat = dynamic_cast< TPZBlackOil2P3D *>(matit->second.operator->());
		if (blackoilmat){
			blackoilmat->SetCurrentState();
		}
	}
}

void TPZBlackOilAnalysis::SetAllMaterialsDeltaT(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZBlackOil2P3D * blackoilmat = dynamic_cast< TPZBlackOil2P3D *>(matit->second.operator->());
		if (blackoilmat){
			blackoilmat->SetTimeStep(this->TimeStep());
		}
	}
}


void TPZBlackOilAnalysis::PostProcess(int resolution, int dimension){
    REAL T = this->fSimulationTime;
    this->fTime = T;
    TPZAnalysis::PostProcess(resolution, dimension);
}//method


void TPZBlackOilAnalysis::PostProcess(TPZVec<REAL> &loc, std::ostream &out){
    REAL T = this->fSimulationTime;
    out << "\nSOLUTION #" << this->fCurrentStep << " AT TIME = " << T << std::endl;
    TPZAnalysis::PostProcess(loc, out);
    out << "\n***************************************\n" << std::endl;
}//method

void TPZBlackOilAnalysis::Assemble(){
	if(!fCompMesh || !fStructMatrix || !fSolver){
		cout << "TPZBlackOilAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh 
		<< " fStructMatrix " << (bool) fStructMatrix << " fSolver " << (bool) fSolver << " at file " 
		<< __FILE__ << " line " << __LINE__ << endl;
		return;
	}
	
	int sz = fCompMesh->NEquations();
	fRhs.Redim(sz,1);
	
	//
	bool exist = false;
	if(fSolver->Matrix()) if (fSolver->Matrix()->Rows()==sz) exist = true;
	if (exist){
		fSolver->Matrix()->Zero();
// #warning FIX ME
		//    fStructMatrix->Assemble(fSolver->Matrix(),fRhs);
	}
	else{
// #warning FIX ME
		//    TPZMatrix<REAL> *mat = fStructMatrix->CreateAssemble(fRhs);
		//    fSolver->SetMatrix(mat);
	}
	fSolver->UpdateFrom(fSolver->Matrix());
}

void TPZBlackOilAnalysis::SetConvergence(int niter, REAL eps, bool ForceAllSteps){
	this->fNIter = niter;
	this->fSteadyTol = eps;
	this->fForceAllSteps = ForceAllSteps;
}

void TPZBlackOilAnalysis::SetSaveFrequency(int SaveFrequency, int resolution){
	this->fSaveFrequency = SaveFrequency;
	this->fDXResolution = resolution;
}

void TPZBlackOilAnalysis::SetNewtonConvergence(int niter, REAL eps){
	this->fNewtonMaxIter = niter;
	this->fNewtonTol = eps;
}

REAL & TPZBlackOilAnalysis::TimeStep(){
	return this->fTimeStep;
}

void TPZBlackOilAnalysis::Solve(){
	
	const int n = this->Solver().Matrix()->Rows();
	double minP = 0, maxP = 0, minS = 0, maxS = 0, p, S;
	for(int i = 0; i < n/2; i++){
		p = this->Solver().Matrix()->operator()(2*i,2*i);
		S = this->Solver().Matrix()->operator()(2*i+1,2*i+1);
		if(p > maxP) maxP = p;
		if(p < minP) minP = p;
		if(S > maxS) maxS = S;
		if(S < minS) minS = S;
	}//for i
	
	double ScaleP = fabs(minP+maxP)/2.;
	double ScaleS = fabs(minS+maxS)/2.;
	for(int j = 0; j < n/2; j++){
		for(int i = 0; i < n; i++){
			this->Solver().Matrix()->operator()(i,2*j) *= 1./ScaleP;
			this->Solver().Matrix()->operator()(i,2*j+1) *= 1./ScaleS;
		}//i
	}//j
	
	{//DEBUG
		double minP = 0, maxP = 0, minS = 0, maxS = 0, p, S;
		for(int i = 0; i < n/2; i++){
			p = this->Solver().Matrix()->operator()(2*i,2*i);
			S = this->Solver().Matrix()->operator()(2*i+1,2*i+1);
			if(p > maxP) maxP = p;
			if(p < minP) minP = p;
			if(S > maxS) maxS = S;
			if(S < minS) minS = S;
		}//for i
//		double ScaleP = fabs(minP+maxP)/2.;
//		double ScaleS = fabs(minS+maxS)/2.;
	}
	
	TPZNonLinearAnalysis::Solve();
	
	for(int i = 0; i < n/2; i++){
		fSolution(2*i,0) *= 1./ScaleP;
		fSolution(2*i+1,0) *= 1./ScaleS;
	}//i
	
}//method

#include "TPZInterfaceEl.h"
double TPZBlackOilAnalysis::PressaoMedia(TPZBlackOilAnalysis &an, int matid){
	an.LoadSolution(an.Solution());
	TPZCompMesh * cmesh = an.Mesh();
	const int nel = cmesh->NElements();
	TPZVec<REAL> qsi(3), sol(1);
	double press = 0.;
	double AccVol = 0.;
	double locVol = 0.;
	for(int iel = 0; iel < nel; iel++){
		TPZCompEl * cel = cmesh->ElementVec()[iel];
		if(!cel) continue;
		if(cel->Material()->Id() != matid) continue;
		TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(cel);
		if(!face) continue;
		face->LeftElement()->Reference()->CenterPoint(face->LeftElement()->Reference()->NSides()-1,qsi);
		face->LeftElement()->Solution(qsi, TPZBlackOil2P3D::EOilPressure, sol);
		locVol = face->LeftElement()->Reference()->Volume();
		press += sol[0] * locVol;
		AccVol += locVol;
	}//iel
	double result = press/AccVol;
	return result;
}//method

#include "pzbndcond.h"
void TPZBlackOilAnalysis::Vazao(TPZBlackOilAnalysis &an, int matid, double & VazaoAguaSC, double  & VazaoOleoSC, double & VazaoAguaFundo, double  & VazaoOleoFundo){
	
	an.LoadSolution(an.Solution());
	TPZCompMesh * cmesh = an.Mesh();
	TPZVec<REAL> qsi(3), sol(1);
	
	TPZElementMatrix ek(cmesh, TPZElementMatrix::EK), ef(cmesh, TPZElementMatrix::EF);
	
	const int nel = cmesh->NElements();
	VazaoAguaSC = 0.;
	VazaoOleoSC = 0.;
	VazaoAguaFundo = 0.;
	VazaoOleoFundo = 0.;
	for(int iel = 0; iel < nel; iel++){
		TPZCompEl * cel = cmesh->ElementVec()[iel];
		if(!cel) continue;
		if(cel->Material()->Id() != matid) continue;
		TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(cel);
		if(!face) continue;
		face->CalcStiff(ek, ef);
		
		face->LeftElement()->Reference()->CenterPoint(face->LeftElement()->Reference()->NSides()-1,qsi);
		face->LeftElement()->Solution(qsi, TPZBlackOil2P3D::EOilPressure, sol);
		
		TPZBlackOil2P3D::BFadREAL po(sol[0],0);
		TPZBlackOil2P3D::BFadREAL Bo;
		
		TPZBndCond * bc = dynamic_cast<TPZBndCond*> (face->Material().operator->());
		TPZBlackOil2P3D * bo = dynamic_cast<TPZBlackOil2P3D *> (bc->Material().operator->());
		bo->Bo(po, Bo);
		
		VazaoOleoSC += ef.fMat(0,0)*86400.;
		VazaoAguaSC += ef.fMat(1,0)*86400.;
		
		VazaoOleoFundo += ef.fMat(0,0)*Bo.val()*86400.;
		VazaoAguaFundo += ef.fMat(1,0)*bo->Bw()*86400.;
	}//iel
}//method

#endif
