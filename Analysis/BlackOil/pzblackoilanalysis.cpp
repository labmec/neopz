/**
 * \file
 * @brief Contains implementations of the TPZBlackOilAnalysis methods.
 */

#include "pzblackoilanalysis.h"
#include "BlackOil/TPZBlackOil2P3D.h"
#include "TPZSpStructMatrix.h"
#include "pzseqsolver.h"
#include "checkconv.h"
#include "TPZElementMatrixT.h"
#include "pzcmesh.h"

using namespace std;
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

void TPZBlackOilAnalysis::SetInitialSolution(TPZFMatrix<STATE> & InitialSol){
	const int nrows = this->Mesh()->Solution().Rows();
	const int ncols = this->Mesh()->Solution().Cols();
	if ( (InitialSol.Rows() != nrows) || (InitialSol.Cols() != ncols) ){
		PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << std::endl;
	}
	else{
		this->fSolution = InitialSol;
	}
	TPZLinearAnalysis::LoadSolution();
}

void TPZBlackOilAnalysis::SetInitialSolutionAsZero(){
	TPZFMatrix<STATE> & MeshSol = this->Mesh()->Solution();
	this->fSolution.Redim( MeshSol.Rows(), MeshSol.Cols() );
	this->fSolution.Zero();
}

void TPZBlackOilAnalysis::AssembleResidual(){
	this->SetCurrentState();
	int sz = this->Mesh()->NEquations();
	TPZFMatrix<STATE> &rhs = this->Rhs();
	rhs.Redim(sz,1);
	rhs += fLastState;
}//void

void TPZBlackOilAnalysis::Run(std::ostream &out, bool linesearch){
	
	ofstream File("Pocos.txt");
	File << "t(mes)\tPi(Pa)\tPp(Pa)\tQwiSC(m3/day)\tQoiSC(m3/day)\tQwpSC(m3/day)\tQopSC(m3/day)\t\tQwiFundo(m3/day)\tQoiFundo(m3/day)\tQwpFundo(m3/day)\tQopFundo(m3/day)\n";
	
	this->fSimulationTime = 0.;
	this->PostProcess(this->fDXResolution);
	
	double nextDeltaT = this->fTimeStep;
	this->SetAllMaterialsDeltaT();
	const double TotalTime = fTimeStep*fNIter;
	
	TPZFMatrix<STATE> prevsol, lastsol;
	
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
		TPZFMatrix<STATE> &myrhs = this->Rhs();
		TPZFMatrix<STATE> &mysol = this->Solution();
		while(error > this->fNewtonTol && iter < this->fNewtonMaxIter){
			
			fSolution.Redim(0,0);
			this->Assemble();			
			myrhs += fLastState;
			
			this->Solve();
			
			if (linesearch){
				TPZFMatrix<STATE> nextSol;
				REAL LineSearchTol = 1e-3 * Norm(fSolution);
				const int niter = 3;
				this->LineSearch(prevsol, fSolution, nextSol, LineSearchTol, niter);
				fSolution = nextSol;
			}
			else{
				
				mysol += prevsol;
			}
			
			prevsol -= fSolution;
			REAL norm = Norm(prevsol);
			out << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << endl;
			
			prevsol = fSolution;
			TPZLinearAnalysis::LoadSolution();
			
			error = norm;
			iter++;
			
			if((iter % 20) == 0){
				//Computing residual of last state solution
				fSolution = prevsol;
				TPZLinearAnalysis::LoadSolution();
				double fator = 0.1;//(sqrt(5.)-1.)/2.;
				this->TimeStep() *= fator;
				nextDeltaT = this->TimeStep();
				cout << "\nMultiplicando passo de tempo por " << fator << "\n";
				this->SetAllMaterialsDeltaT();
				this->SetLastState();
				this->Assemble();
				fLastState = this->fRhs;
				this->SetCurrentState();
				this->fNIter = fCurrentStep + ((int)((1./fator) * (fNIter - fCurrentStep) + 0.5));
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
				double WaterFlowIsc, WaterFlowPsc, OilFlowIsc, OilFlowPsc;
				double WaterFlowIFundo, WaterFlowPFundo, OilFlowIFundo, OilFlowPFundo;
				Flow(*this,-1,WaterFlowIsc,OilFlowIsc,WaterFlowIFundo,OilFlowIFundo);
				Flow(*this,-2,WaterFlowPsc,OilFlowPsc,WaterFlowPFundo,OilFlowPFundo);
				File << this->fSimulationTime/2629743.8 << "\t" << AveragePressure(*this,-1) << "\t" << AveragePressure(*this,-2) << "\t"
				<< WaterFlowIsc << "\t" 
				<< OilFlowIsc << "\t" 
				<< WaterFlowPsc << "\t" 
				<< OilFlowPsc << "\t" << "\t"
				<< WaterFlowIFundo << "\t" 
				<< OilFlowIFundo << "\t" 
				<< WaterFlowPFundo << "\t" 
				<< OilFlowPFundo << "\n";
				File.flush();
			}
		}
		
		prevsol -= lastsol;
		REAL steadynorm = Norm(prevsol);
		std::cout << "*********** Steady state error at iteration " << this->fCurrentStep << ", " << this->fSimulationTime/2629743.8 << " months " << " = " << steadynorm << "\n\n";
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
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZBlackOil2P3D * blackoilmat = dynamic_cast< TPZBlackOil2P3D *>(matit->second);
		if (blackoilmat){
			blackoilmat->SetLastState();
		}
	}
}


void TPZBlackOilAnalysis::SetCurrentState(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZBlackOil2P3D * blackoilmat = dynamic_cast< TPZBlackOil2P3D *>(matit->second);
		if (blackoilmat){
			blackoilmat->SetCurrentState();
		}
	}
}

void TPZBlackOilAnalysis::SetAllMaterialsDeltaT(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZBlackOil2P3D * blackoilmat = dynamic_cast< TPZBlackOil2P3D *>(matit->second);
		if (blackoilmat){
			blackoilmat->SetTimeStep(this->TimeStep());
		}
	}
}


void TPZBlackOilAnalysis::PostProcess(int resolution, int dimension){
    REAL T = this->fSimulationTime;
    this->fTime = T;
    TPZLinearAnalysis::PostProcess(resolution, dimension);
}//method


void TPZBlackOilAnalysis::PostProcess(TPZVec<REAL> &loc, std::ostream &out){
    REAL T = this->fSimulationTime;
    out << "\nSOLUTION #" << this->fCurrentStep << " AT TIME = " << T << std::endl;
    TPZLinearAnalysis::PostProcess(loc, out);
    out << "\n***************************************\n" << std::endl;
}//method

void TPZBlackOilAnalysis::Assemble(){
	auto &solver = MatrixSolver<STATE>();
	auto solverMat = solver.Matrix();
	if(!fCompMesh || !fStructMatrix || !fSolver){
		cout << "TPZBlackOilAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh 
		<< " fStructMatrix " << (bool) fStructMatrix << " fSolver " << fSolver << " at file " 
		<< __FILE__ << " line " << __LINE__ << endl;
		return;
	}
	
	int sz = fCompMesh->NEquations();
	fRhs.Redim(sz,1);
	
	//
	bool exist = false;
	if(solverMat) if (solverMat->Rows()==sz) exist = true;
	if (exist){
		solverMat->Zero();
	}
	solver.UpdateFrom(solverMat);
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

	auto solverMat = MatrixSolver<STATE>().Matrix();
	const int n = solverMat->Rows();
	double minP = 0, maxP = 0, minS = 0, maxS = 0, p, S;
	for(int i = 0; i < n/2; i++){
		p = solverMat->operator()(2*i,2*i);
		S = solverMat->operator()(2*i+1,2*i+1);
		if(p > maxP) maxP = p;
		if(p < minP) minP = p;
		if(S > maxS) maxS = S;
		if(S < minS) minS = S;
	}//for i
	
	double ScaleP = fabs(minP+maxP)/2.;
	double ScaleS = fabs(minS+maxS)/2.;
	for(int j = 0; j < n/2; j++){
		for(int i = 0; i < n; i++){
			solverMat->operator()(i,2*j) *= 1./ScaleP;
			solverMat->operator()(i,2*j+1) *= 1./ScaleS;
		}//i
	}//j
	
	{//DEBUG
		double minP = 0, maxP = 0, minS = 0, maxS = 0, p, S;
		for(int i = 0; i < n/2; i++){
			p = solverMat->operator()(2*i,2*i);
			S = solverMat->operator()(2*i+1,2*i+1);
			if(p > maxP) maxP = p;
			if(p < minP) minP = p;
			if(S > maxS) maxS = S;
			if(S < minS) minS = S;
		}//for i
	}
	
	TPZNonLinearAnalysis::Solve();

	TPZFMatrix<STATE> &mysol = fSolution;
	for(int i = 0; i < n/2; i++){
		mysol(2*i,0) *= 1./ScaleP;
		mysol(2*i+1,0) *= 1./ScaleS;
	}//i
	
}//method

#include "TPZInterfaceEl.h"
double TPZBlackOilAnalysis::AveragePressure(TPZBlackOilAnalysis &an, int matid){
	an.LoadSolution(an.Solution());
	TPZCompMesh * cmesh = an.Mesh();
	const int nel = cmesh->NElements();
	TPZVec<REAL> qsi(3);
    TPZVec<STATE> sol(1);
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


void TPZBlackOilAnalysis::Flow(TPZBlackOilAnalysis &an, int matid, double & WaterFlowSC, double  & OilFlowSC, double & WaterFlowBottom, double  & OilFlowBottom){
	
	an.LoadSolution(an.Solution());
	TPZCompMesh * cmesh = an.Mesh();
	TPZVec<REAL> qsi(3);
    TPZVec<STATE> sol(1);
	
	TPZElementMatrixT<STATE> ek(cmesh, TPZElementMatrix::EK), ef(cmesh, TPZElementMatrix::EF);
	const int nel = cmesh->NElements();
	WaterFlowSC = 0.;
	OilFlowSC = 0.;
	WaterFlowBottom = 0.;
	OilFlowBottom = 0.;
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
		
		auto * bc =
			dynamic_cast<TPZBndCondT<STATE>*> (face->Material());
		auto * bo =
			dynamic_cast<TPZBlackOil2P3D *> (bc->Material());
		bo->Bo(po, Bo);
		
		OilFlowSC += ef.fMat(0,0)*86400.;
		WaterFlowSC += ef.fMat(1,0)*86400.;
		
		OilFlowBottom += ef.fMat(0,0)*Bo.val()*86400.;
		WaterFlowBottom += ef.fMat(1,0)*bo->Bw()*86400.;
	}//iel
}//method

