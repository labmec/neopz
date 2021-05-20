/**
 * @file
 * @brief Contains implementations of the TPZTransientAnalysis methods.
 */

#include "pztransientanalysis.h"
#include "TPZMatTransientSingleSpace.h"
#include "TPZSpStructMatrix.h"
#include "pzseqsolver.h"
#include "checkconv.h"
#include "pzcmesh.h"

using namespace std;

template<class TRANSIENTCLASS>
double TPZTransientAnalysis<TRANSIENTCLASS>::gTime = 0.;

template<class TRANSIENTCLASS>
TPZTransientAnalysis<TRANSIENTCLASS>::TPZTransientAnalysis(TPZCompMesh *mesh, bool IsLinear, std::ostream &out):/*TPZLinearAnalysis*/TPZNonLinearAnalysis(mesh,out), fSavedSolutionVec(){
	this->fTimeStep = 0.;
	this->fCurrentIter = 0;
	this->SetConvergence(0, 0.);
	this->SetNewtonConvergence(0, 0.);
	this->SetInitialSolutionAsZero();
	this->fIsLinearProblem = IsLinear;
	this->SetSaveFrequency(0,0);
	this->fSaveSolutionVecFrequency = 0; 
}

template<class TRANSIENTCLASS>
TPZTransientAnalysis<TRANSIENTCLASS>::~TPZTransientAnalysis(){
	fSavedSolutionVec.clear();
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetInitialSolution(TPZFMatrix<STATE> & InitialSol){
	const int nrows = this->Mesh()->Solution().Rows();
	const int ncols = this->Mesh()->Solution().Cols();
	if ( (InitialSol.Rows() != nrows) || (InitialSol.Cols() != ncols) ){
		PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << std::endl;
	}
	else{
		this->fSolution = InitialSol;
	}
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetInitialSolutionAsZero(){
	TPZFMatrix<STATE> & MeshSol = this->Mesh()->Solution();
	this->fSolution.Redim( MeshSol.Rows(), MeshSol.Cols() );
	this->fSolution.Zero();
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::RunTransient(std::ostream &out, bool FromBegining, bool linesearch){
	
	this->SetImplicit();
	
	if (FromBegining){
		this->fCurrentIter = 0;
		this->fSavedSolutionVec.clear();
	}
	
#ifdef _DOCHECKCONV_
	{
		TPZVec<REAL> coefs(1,1.);
		TPZFMatrix<STATE> cpSol(fSolution);
		TPZFMatrix<STATE> range(fCompMesh->NEquations(),1,1.);
		this->SetLastState();
		CheckConvergence(*this,cpSol,range,coefs);
		this->SetCurrentState();
		CheckConvergence(*this,cpSol,range,coefs);
	}
#endif
	
	this->LoadSolution(fSolution);
	
	TPZTransientAnalysis::gTime =this->GetCurrentIter() * this->TimeStep();
	//   this->PostProcess(this->fDXResolution);
	
	this->SetAllMaterialsDeltaT();
	
	if (this->fIsLinearProblem){
		this->ComputeLinearTangentMatrix();  
	}
	
	TPZFMatrix<STATE> prevsol, laststate, lastsol;
	for( ; this->fCurrentIter < this->fNIter;){
		
		this->fCurrentIter++;
		TPZTransientAnalysis::gTime = this->TimeStep() * this->fCurrentIter;
		
		//Computing residual of last state solution
		//     this->fSolution = prevsol;
		this->SetLastState();
		this->Assemble();
		laststate = this->fRhs;
		prevsol = fSolution;
		lastsol = fSolution;
		//Newton's method
		this->SetCurrentState();        
		REAL error = this->fNewtonTol * 2. + 1.;    
		int iter = 0;
		while(error > this->fNewtonTol && iter < this->fNewtonMaxIter) {
			
			fSolution.Redim(0,0);
			this->Assemble();
			TPZFMatrix<STATE> &rhs = this->fRhs;
			rhs += laststate;
			this->Solve();
			
			if (linesearch){
				TPZFMatrix<STATE> nextSol;
				REAL LineSearchTol = 1e-3 * Norm(fSolution);
				const int niter = 100;
				this->LineSearch(prevsol, fSolution, nextSol, LineSearchTol, niter);
				fSolution = nextSol;
			}
			else{
				TPZFMatrix<STATE> &sol = fSolution;
				sol += prevsol;
			}
			
			prevsol -= fSolution;
			REAL norm = Norm(prevsol);
			out << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << std::endl;
			
			prevsol = fSolution;
			TPZLinearAnalysis::LoadSolution();
			
			error = norm;
			iter++;
		}//Newton's iterations
		
		prevsol = fSolution;
		
		if (this->fSaveFrequency){
			if (!(this->fCurrentIter % fSaveFrequency)){
				this->PostProcess(this->fDXResolution);
			}
		}
		this->SaveCurrentSolutionVec();
		
		prevsol -= lastsol;
		REAL steadynorm = Norm(prevsol);
		std::cout << "*********** Steady state error at iteration " << this->fCurrentIter << " = " << steadynorm << "\n\n";
		if (!fForceAllSteps){
			if (steadynorm < this->fSteadyTol){
				std::cout << "Steady state solution achieved\n\n";
				this->fNIter = fCurrentIter;
				break;
			}
		}
		std::cout.flush();   
		
	}//time iterations
    
}//method

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetImplicit(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZMatTransientSingleSpace< TRANSIENTCLASS > * trans = dynamic_cast<TPZMatTransientSingleSpace< TRANSIENTCLASS > *>(matit->second);
		if (trans){
			trans->SetImplicit();
		}
	}
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetExplicit(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZMatTransientSingleSpace< TRANSIENTCLASS > * trans = dynamic_cast<TPZMatTransientSingleSpace< TRANSIENTCLASS > *>(matit->second);
		if (trans){
			trans->SetExplicit();
		}
	}
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetLastState(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZMatTransientSingleSpace< TRANSIENTCLASS > * trans = dynamic_cast<TPZMatTransientSingleSpace< TRANSIENTCLASS > *>(matit->second);
		if (trans){
			trans->SetLastState();
		}
	}
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetCurrentState(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZMatTransientSingleSpace< TRANSIENTCLASS > * trans = dynamic_cast<TPZMatTransientSingleSpace< TRANSIENTCLASS > *>(matit->second);
		if (trans){
			trans->SetCurrentState();
		}
	}
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetMassMatrix(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZMatTransientSingleSpace< TRANSIENTCLASS > * trans = dynamic_cast<TPZMatTransientSingleSpace< TRANSIENTCLASS > *>(matit->second);
		if (trans){
			trans->SetMassMatrix();
		}
	}
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetFluxOnly(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZMatTransientSingleSpace< TRANSIENTCLASS > * trans = dynamic_cast<TPZMatTransientSingleSpace< TRANSIENTCLASS > *>(matit->second);
		if (trans){
			trans->SetFluxOnly();
		}
	}
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetAllMaterialsDeltaT(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZMatTransientSingleSpace< TRANSIENTCLASS > * trans = dynamic_cast<TPZMatTransientSingleSpace< TRANSIENTCLASS > *>(matit->second);
		if (trans){
			trans->SetTimeStep(this->TimeStep());
		}
	}
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::PostProcess(int resolution, int dimension){
    REAL T = this->GetCurrentIter() * this->TimeStep();
    this->fTime = T;
    TPZLinearAnalysis::PostProcess(resolution, dimension);
}//method

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::PostProcess(TPZVec<REAL> &loc, std::ostream &out){
    REAL T = this->GetCurrentIter() * this->TimeStep();
    this->gTime = T;
    out << "\nSOLUTION #" << this->GetCurrentIter() << " AT TIME = " << T << std::endl;
    TPZLinearAnalysis::PostProcess(loc, out);
    out << "\n***************************************\n" << std::endl;
}//method


template<class TRANSIENTCLASS>
template<class TVar>
void TPZTransientAnalysis<TRANSIENTCLASS>::AssembleInternal()
{
	auto &mySolver = MatrixSolver<TVar>();
	auto solverMat = mySolver.Matrix();
	if(!fCompMesh || !fStructMatrix || !fSolver){
		cout << "TPZTransientAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh 
		<< " fStructMatrix " << (bool) fStructMatrix << " fSolver " << fSolver << " at file " 
		<< __FILE__ << " line " << __LINE__ << endl;
		return;
	}
	
	int sz = fCompMesh->NEquations();
	fRhs.Redim(sz,1);
	
	bool exist = false;
	if(solverMat) if (solverMat->Rows()==sz) exist = true;
	TPZAutoPointer<TPZGuiInterface> inter = new TPZGuiInterface;
	if (exist){
		if (fIsLinearProblem){
			//      TPZStructMatrix::Assemble(fRhs, *Mesh());
			fStructMatrix->Assemble(fRhs,inter);
		}
		else{
			solverMat->Zero();
			fStructMatrix->Assemble((TPZMatrix<TVar>&)solverMat,fRhs,inter);
		}
	}
	else{
		if (this->fIsLinearProblem){
			std::cout << __PRETTY_FUNCTION__ << " @ " << __LINE__ << " Error! StrMatrix must be created using" 
			<< " methodTPZTransientAnalysis::ComputeLinearTangentMatrix()"
			<< " when (this->fIsLinearProblem == true)\n";
		}
		auto *mat = fStructMatrix->CreateAssemble(fRhs,NULL);
		mySolver.SetMatrix(mat);
	}
	mySolver.UpdateFrom(solverMat);
}
template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::Assemble(){
	//TODOCOMPLEX
	AssembleInternal<STATE>();
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::ComputeLinearTangentMatrix(){
	if (!fIsLinearProblem) return;      
	this->SetCurrentState();
	const int sz = this->Mesh()->NEquations();
	fRhs.Redim(sz,1);
	//TODOCOMPLEX
	auto *mat = fStructMatrix->CreateAssemble(fRhs,NULL);
	MatrixSolver<STATE>().SetMatrix(mat);
}//method

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::ComputeMassMatrix(){
	this->SetMassMatrix();
	const int sz = this->Mesh()->NEquations();
	fRhs.Redim(sz,1);
	//TODOCOMPLEX
	auto *mat = fStructMatrix->CreateAssemble(fRhs,NULL);
	MatrixSolver<STATE>().SetMatrix(mat);
}//method

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::ComputeFluxOnly(){
	//TODOCOMPLEX
	if(!fCompMesh || !fStructMatrix || !fSolver){
		cout << "TPZTransientAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh 
		<< " fStructMatrix " << (bool) fStructMatrix << " fSolver " << fSolver << " at file " 
		<< __FILE__ << " line " << __LINE__ << endl;
		return;
	}
	
	this->SetFluxOnly();  
	int sz = fCompMesh->NEquations();
	fRhs.Redim(sz,1);
	auto solverMat = MatrixSolver<STATE>().Matrix();
	if(solverMat && solverMat->Rows()==sz){
		fStructMatrix->Assemble(fRhs,NULL);
		//    TPZStructMatrix::Assemble(fRhs, *Mesh());
	}//if
}//method

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::RunExplicit(std::ostream &out, bool FromBegining){
	
	this->SetExplicit();
	
	if (FromBegining){
		this->fCurrentIter = 0;
		this->fSavedSolutionVec.clear();
	}
	TPZTransientAnalysis::gTime = this->TimeStep() * this->fCurrentIter;
	this->PostProcess(this->fDXResolution);
	
	this->SetAllMaterialsDeltaT();
	
	TPZFMatrix<STATE> prevsol;
	for( this->fCurrentIter++ ; this->fCurrentIter < this->fNIter; this->fCurrentIter++){
		
		this->ComputeMassMatrix();
		
		TPZTransientAnalysis::gTime = this->TimeStep() * this->fCurrentIter;
		
		this->SetFluxOnly();
		
		//Computing residual of last state solution
		prevsol = fSolution;
		TPZLinearAnalysis::LoadSolution();
		this->ComputeFluxOnly();
		
		this->Solve();
		//now fSolution = deltaSol
		TPZFMatrix<STATE> &sol = fSolution;
		sol += prevsol;
		
		TPZLinearAnalysis::LoadSolution();
		if (this->fSaveFrequency){
			if (!(this->fCurrentIter % fSaveFrequency)){
				this->PostProcess(this->fDXResolution);
			}
		}
		this->SaveCurrentSolutionVec();
		
		prevsol -= fSolution;
		REAL steadynorm = Norm(prevsol);
		std::cout << "*********** Steady state error at iteration " << (this->fCurrentIter) << " = " << steadynorm << "\n\n";
		if (!fForceAllSteps){
			if (steadynorm < this->fSteadyTol){
				std::cout << "Steady state solution achieved\n\n";
				this->fNIter = fCurrentIter;
				break;
			}
		}
		std::cout.flush();   
		
	}//time iterations
	
}//method

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetSaveSolution(int SaveFrequency){
	this->fSaveSolutionVecFrequency = SaveFrequency;
}

template<class TRANSIENTCLASS>
std::list< std::pair<TPZFMatrix<STATE>, STATE> > & TPZTransientAnalysis<TRANSIENTCLASS>::GetSavedSolutions(){
	return this->fSavedSolutionVec;
}

#include <sstream>
template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SaveCurrentSolutionVec(){
	if(!this->fSaveSolutionVecFrequency) return;
	if(this->fCurrentIter % this->fSaveSolutionVecFrequency == 0){
		std::pair< TPZFMatrix<STATE>, STATE > mypair;
		mypair.first = this->Solution();
		mypair.second = TPZTransientAnalysis::gTime;
		this->fSavedSolutionVec.push_back(mypair);
		
		ofstream file("currentsol.txt");
		stringstream mess; mess << "sol( " << TPZTransientAnalysis::gTime << " ) = ";
		TPZBaseMatrix &mysol = this->Solution();
		mysol.Print(mess.str().c_str(), file);
		
	}
}