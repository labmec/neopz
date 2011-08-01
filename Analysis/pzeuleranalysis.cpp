/**
 * \file
 * @brief Contains implementations of the TPZEulerAnalysis methods.
 */
//$Id: pzeuleranalysis.cpp,v 1.50 2011-05-11 01:14:38 phil Exp $

#include "pzeuleranalysis.h"
#include "pzerror.h"
#include "TPZCompElDisc.h"
#include "pzfstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPBSpStructMatrix.h"
#include "pzbdstrmatrix.h"

#include "tpzoutofrange.h"
#include <time.h>
#include "pzlog.h"

#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("pz.converge"));
#endif

using namespace std;

TPZEulerAnalysis::TPZEulerAnalysis():
TPZAnalysis(), fFlowCompMesh(NULL),
fRhsLast(),
fLinSysEps(1e-10), fLinSysMaxIter(20),
fNewtonEps(1e-9),  fNewtonMaxIter(10),
fTimeIntEps(1e-8), fTimeIntMaxIter(100),
fEvolCFL(0), fpBlockDiag(NULL),fHasFrontalPreconditioner(0)
{
	
}

TPZEulerAnalysis::TPZEulerAnalysis(TPZFlowCompMesh *mesh, std::ostream &out):
TPZAnalysis(mesh, out), fFlowCompMesh(mesh),
fRhsLast(),
fLinSysEps(1e-10), fLinSysMaxIter(20),
fNewtonEps(1e-9),  fNewtonMaxIter(10),
fTimeIntEps(1e-8), fTimeIntMaxIter(100),
fEvolCFL(0), fpBlockDiag(NULL),fHasFrontalPreconditioner(0)
{
	
}

TPZEulerAnalysis::~TPZEulerAnalysis()
{
}

void TPZEulerAnalysis::SetAdvancedState()
{
	SetContributionTime(Advanced_CT);
	fCompMesh->LoadSolution(fSolution);
}

void TPZEulerAnalysis::SetLastState()
{
	SetContributionTime(Last_CT);
	fCompMesh->LoadSolution(fSolution);
}

void TPZEulerAnalysis::SetContributionTime(TPZContributeTime time)
{
	fFlowCompMesh->SetContributionTime(time);
}

void TPZEulerAnalysis::UpdateSolAndRhs(TPZFMatrix & deltaSol, REAL & epsilon)
{
    REAL initEpsilon = epsilon;
    int outofrange = 0;
    try
    {
        fSolution += deltaSol;
        fCompMesh->LoadSolution(fSolution);
        AssembleRhs();
        epsilon = Norm(fRhs);
    }
	catch(TPZOutofRange obj)
	{
		outofrange = 1;
		fSolution -= deltaSol;
		epsilon = initEpsilon;
		fCompMesh->LoadSolution(fSolution);
	}
	
    if(epsilon > initEpsilon)
    {
		fSolution -= deltaSol;
		fCompMesh->LoadSolution(fSolution);
		/*int resultlin = */
		LineSearch(initEpsilon ,fSolution, deltaSol);
		fSolution += deltaSol;
		fCompMesh->LoadSolution(fSolution);
		epsilon = Norm(fRhs);
    }
	
    if(outofrange)
    {
		/*int resultlin = */LineSearch(initEpsilon ,fSolution, deltaSol);
		fSolution += deltaSol;
		fCompMesh->LoadSolution(fSolution);
		epsilon = Norm(fRhs);
		fFlowCompMesh->ScaleCFL(.5);
    }
}

void TPZEulerAnalysis::UpdateHistory()
{
	
#ifdef RESTART_ZEROED
	//Zeroeing the newest iterative solution
	fSolution.Zero();
#else
	// The last state should be copied to the advanced
	// state vector.
	// These are in fact the same storage, so that the
	// memory desn't need to be copied.
#endif
}

void TPZEulerAnalysis::BufferLastStateAssemble()
{
	fRhsLast.Zero();
	fRhsLast.Redim(fCompMesh->NEquations(),1);
	SetLastState();
	fStructMatrix->Assemble(fRhsLast,NULL);
	//   TPZStructMatrix::Assemble(fRhsLast, *fFlowCompMesh);
	SetAdvancedState();
	UpdateHistory();
}

REAL TPZEulerAnalysis::EvaluateFluxEpsilon()
{
	// enabling the flux evaluation type
	fFlowCompMesh->SetResidualType(Flux_RT);
	TPZFMatrix Flux;
	Flux.Redim(fCompMesh->NEquations(),1);
	Flux.Zero();
	// contribution of last state
	SetLastState();
	fStructMatrix->Assemble(Flux,NULL);
	//   TPZStructMatrix::Assemble(Flux, *fFlowCompMesh);
	
	// constribution of adv state
	SetAdvancedState();
	fStructMatrix->Assemble(Flux,NULL);
	//   TPZStructMatrix::Assemble(Flux, *fFlowCompMesh);
	
	// reseting the residual type to residual
	fFlowCompMesh->SetResidualType(Residual_RT);
	
	return Norm(Flux);
}

void TPZEulerAnalysis::Assemble()
{
	if(!fCompMesh)
	{
		PZError << "TPZEulerAnalysis::Assemble Error: No Computational Mesh\n";
		//      return;
		DebugStop();
	}
	
	if(!fStructMatrix)
	{
		PZError << "TPZEulerAnalysis::Assemble Error: No Structural Matrix\n";
		DebugStop();
		//      return;
	}
	
	if(!fSolver)
	{
		PZError << "TPZEulerAnalysis::Assemble Error: No Solver\n";
		DebugStop();
		//      return;
	}
	
	// contributing referring to the last state
	fRhs = fRhsLast;
	
	TPZAutoPointer<TPZMatrix>  pTangentMatrix = fSolver->Matrix();
	
	if(!pTangentMatrix || dynamic_cast<TPZParFrontStructMatrix <TPZFrontNonSym> *>(fStructMatrix.operator->()))
	{
		pTangentMatrix = fStructMatrix->CreateAssemble(fRhs,NULL);
		fSolver->SetMatrix(pTangentMatrix);
	}
	else
	{
		if(!pTangentMatrix)
		{
			PZError << "TPZEulerAnalysis::Assemble Error: No Structural Matrix\n";
			DebugStop();
			//         return;
			
		}
		
		pTangentMatrix->Zero();
		
		// Contributing referring to the advanced state
		// (n+1 index)
		fStructMatrix->Assemble(pTangentMatrix, fRhs,NULL);
	}
	
	if(fpBlockDiag)
	{
		fpBlockDiag->Zero();
		//fpBlockDiag->SetIsDecomposed(0); // Zero already makes it
		fpBlockDiag->BuildFromMatrix(pTangentMatrix);
	}
	
	/*
	 ofstream Mout("Matriz.out");
	 ofstream Vout("Vetor.out");
	 
	 pTangentMatrix->Print("Matrix", Mout);//EMathematicaInput);
	 
	 fRhs.Print("Rhs", Vout);
	 
	 Mout.close();
	 Vout.close();*/
}


void TPZEulerAnalysis::AssembleRhs()
{
	if(!fCompMesh) return;
	
	// Contributing referring to the last state (n index)
	fRhs = fRhsLast;
	
	// Contributing referring to the advanced state
	// (n+1 index)
	fStructMatrix->Assemble(fRhs,NULL);
	//   TPZStructMatrix::Assemble(fRhs, *fFlowCompMesh);
	
}

//ofstream eulerout("Matrizes.out");

int TPZEulerAnalysis::Solve(REAL & res, TPZFMatrix * residual, TPZFMatrix & delSol) {
	int numeq = fCompMesh->NEquations();
	if(fRhs.Rows() != numeq ) return 0;
	
	TPZFMatrix rhs(fRhs);
	
	fSolver->Solve(rhs, delSol);
	
	if(residual)
	{   // verifying the inversion of the linear system
		residual->Redim(numeq,1);
		/*
		 fSolver->Matrix()->Residual(delSol, fRhs, *residual);
		 res = Norm(*residual);
		 
		 if(res > fLinSysEps)
		 {
		 cout << "Linear system invertion did not achieve expected tolerance:" << res << endl;
		 cout.flush();
		 }
		 */
	}
	
	/*
	 fSolver->Matrix()->Print("Matrix", eulerout, EMathematicaInput);
	 delu.Print("delu", eulerout, EMathematicaInput);
	 fRhs.Print("Rhs", eulerout, EMathematicaInput);
	 eulerout.flush();
	 */
    return 1;
}

int TPZEulerAnalysis::RunNewton(REAL & epsilon, int & numIter)
{
	TPZFMatrix delSol(fRhs.Rows(),1);
	
	int i = 0;
	REAL res;// residual of linear invertion.
	epsilon = fNewtonEps * REAL(2.);// ensuring the loop will be
	// performed at least once.
	
	TPZFMatrix residual;
	//    REAL res1,res2 = 0.;
	if(fHasFrontalPreconditioner)
	{
		TPZFrontStructMatrix <TPZFrontNonSym> StrMatrix(Mesh());
		StrMatrix.SetQuiet(1);
		TPZMatrix *front = StrMatrix.CreateAssemble(fRhs,NULL);
		TPZStepSolver FrontSolver;
		FrontSolver.SetDirect(ELU);
		FrontSolver.SetMatrix(front);
		TPZStepSolver *step = dynamic_cast<TPZStepSolver *>(fSolver);
		step->SetPreconditioner(FrontSolver);
	}
	
	while(i < fNewtonMaxIter && epsilon > fNewtonEps)
	{
		// Linearizes the system with the newest iterative solution
		Assemble();
		if(i==0)
		{
			epsilon = Norm(fRhs);
			cout << "\tEntry NonLinEpsilon:" << epsilon << endl;
		}
		
		// Testing whether the rhs computed by assemble is the same as the rhs computed by AssembleRhs
		// Both are different, because expression templates change the order of the computations
		/*
		 res1 = Norm(fRhs);
		 if(res2 != 0. && res1 != res2)
		 {
		 this->CompareRhs();
		 }
		 */
		
		//Solves the linearized system
		fTotalNewton++;
		if(Solve(res, &residual, delSol) == 0) return 0;
		
		// Updates the solution, attempts to update Rhs (vector only) and
		// returns the nonlinear residual (epsilon).
		// According to the initial and final values of epsilon, the
		// method may perform a line search.
		UpdateSolAndRhs(delSol, epsilon);
		
		cout << "\tNonLinEpsilon:" << epsilon << endl;
		i++;
	}
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Numero de iteracoes de Newton " << i;
		LOGPZ_DEBUG(logger,sout.str().c_str());
	}
#endif
	numIter = i;
	
	if(epsilon > fNewtonEps)return 0;
	return 1;
}

TPZDXGraphMesh * TPZEulerAnalysis::PrepareDXMesh(const std::string &dxout, int dxRes)
{
	TPZVec<std::string> scalar(4),vector(0);
	scalar[0] = "density";
	scalar[1] = "pressure";
	scalar[2] = "normvelocity";
	scalar[3] = "Mach";
	
	
	TPZAutoPointer<TPZMaterial>  mat = fFlowCompMesh->GetFlowMaterial();
	int dim = mat->Dimension();
	//ResetReference(Mesh());//retira refer�ncias para criar graph consistente
	TPZDXGraphMesh * graph = new TPZDXGraphMesh (Mesh(),dim,mat,scalar,vector);
	//SetReference(Mesh());//recupera as refer�ncias retiradas
	//ofstream *dxout = new ofstream("ConsLaw.dx");
	//cout << "\nDX output file : ConsLaw.dx\n";
	graph->SetFileName(dxout);
	//int resolution = 2;
	graph->SetResolution(dxRes);
	graph->DrawMesh(dim);
	
	return graph;
}

void TPZEulerAnalysis::Run(std::ostream &out, const std::string & dxout, int dxRes)
{
	// this analysis loop encloses several calls
	// to Newton's linearizations, updating the
	// time step.
	
	this->fTotalNewton = 0;
	out << "\nBeginning time integration";
	
	time_t startTime_t = time(NULL);
	
	TPZDXGraphMesh * graph = PrepareDXMesh(dxout, dxRes);
	int numIterDX = 1;
	
	REAL epsilon_Newton/*, epsilon_Global*/;
	int numIter_Newton/*, numIter_Global*/;
	REAL epsilon, lastEpsilon=0;
	
	int i = 0;
	epsilon = 2. * fTimeIntEps;
	
#ifdef RESTART_ZEROED
	fSolution->Zero();
#else
	// the history equals the solution
	// nothing must be done.
#endif
	
	// evaluates the time step based on the solution
	// (last Sol, since Curr sol equals it)
	REAL nextTimeStep = ComputeTimeStep();
	REAL AccumTime = 0.;
	
	// Buffers the contribution of the last state
	BufferLastStateAssemble();
	
	out << "iter\ttime\teps(dU/dt)\tNewtonEps=\tnNewtonIter\n";
	
	while(i < fTimeIntMaxIter && epsilon > fTimeIntEps)
	{
		if(i%numIterDX==0)
		{
			graph->DrawSolution(i,AccumTime);
			graph->Out().flush();
			// increases the accumulated time and
			AccumTime += nextTimeStep;
		}
		
		// Solves the nonlinear system, updates the solution,
		// history and last state assemble buffer.
		/*int newtonreturn = */
		RunNewton(epsilon_Newton, numIter_Newton);
		
		// resetting the forcing function for iterations greater than the 1st
		fFlowCompMesh->SetFlowforcingFunction(NULL);
		
		// buffering the contribution to the RHS
		BufferLastStateAssemble();
		
		// zeroes the fSolution vector if necessary
		UpdateHistory();
		
		// evaluates the norm of fluxes across interfaces
		epsilon = EvaluateFluxEpsilon();
		
		// computing the time step
		nextTimeStep = ComputeTimeStep();
		
		CFLControl(lastEpsilon, epsilon, epsilon_Newton,  nextTimeStep);
		
#ifdef LOG4CXX
		if(logger->isDebugEnabled())
		{
			std::stringstream sout;
			sout << "iteration " << i << "Du/Dt " << epsilon << " Total Newton " << fTotalNewton;
			LOGPZ_DEBUG(logger,sout.str().c_str());
		}    
#endif
		// output
		{
			out << "\n" << i
			<< "\t" << AccumTime
			<< "\t" << epsilon
			<< "\t" << epsilon_Newton
			<< "\t" << numIter_Newton;
			
			cout << "iter:" << i
			<< "\ttime:" << AccumTime
			<< "\t eps(dU/dt)=" << epsilon
			<< "\t |NewtonEps=" << epsilon_Newton
			<< "\t nIter=" << numIter_Newton << endl;
		}
		i++;
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Total number of newton iterations " << this->fTotalNewton;
		LOGPZ_DEBUG(logger,sout.str().c_str());
		
	}
#endif
	
	fSolver->ResetMatrix(); // deletes the memory allocated
	// for the storage of the tangent matrix.
	
	graph->DrawSolution(i,AccumTime);
	graph->Out().flush();
	
	delete graph;
	
	time_t endTime_t = time(NULL);
	
	out  << "\nElapsed time in seconds: " << (endTime_t - startTime_t) << endl;
	
	cout << "\nElapsed time in seconds: " << (endTime_t - startTime_t) << endl;
	
	
}

void TPZEulerAnalysis::CFLControl(REAL & lastEpsilon, REAL & epsilon, REAL & epsilon_Newton, REAL & timeStep)
{
	// CFL control based on Flux reduction
	if((lastEpsilon>0. && epsilon>0.) && fEvolCFL >= 1)
	{
		fFlowCompMesh->ScaleCFL(lastEpsilon/epsilon);
		timeStep *= lastEpsilon/epsilon;
	}
	
	// CFL control based on Nonlinear invertion
	if(epsilon_Newton < fNewtonEps && fEvolCFL >= 2)
	{
		fFlowCompMesh->ScaleCFL(2.);
		timeStep *= 2.;
	}
	
	if(epsilon_Newton > fNewtonEps && fEvolCFL >= 3)
	{
		fFlowCompMesh->ScaleCFL(.5);
		timeStep *= .5;
	}
	
	lastEpsilon = epsilon;
}


REAL TPZEulerAnalysis::ComputeTimeStep()
{
	return fFlowCompMesh->ComputeTimeStep();
}

void TPZEulerAnalysis::SetLinSysCriteria(REAL epsilon, int maxIter)
{
	fLinSysEps     = epsilon;
	fLinSysMaxIter = maxIter;
}

void TPZEulerAnalysis::SetNewtonCriteria(REAL epsilon, int maxIter)
{
	fNewtonEps     = epsilon;
	fNewtonMaxIter = maxIter;
}

void TPZEulerAnalysis::SetTimeIntCriteria(REAL epsilon, int maxIter)
{
	fTimeIntEps     = epsilon;
	fTimeIntMaxIter = maxIter;
}

void TPZEulerAnalysis::SetEvolCFL(int EvolCFL)
{
	fEvolCFL = EvolCFL;
}

void TPZEulerAnalysis::SetBlockDiagonalPrecond(TPZBlockDiagonal * blockDiag)
{
	fpBlockDiag = blockDiag;
}

void TPZEulerAnalysis::WriteCMesh( const char * str)
{
	TPZFileStream fstr;
	fstr.OpenWrite(str);
	fGeoMesh->Write(fstr,1);
	fFlowCompMesh->Write(fstr,1);
}


int TPZEulerAnalysis::LineSearch(REAL &residual, TPZFMatrix &sol0, TPZFMatrix &direction)
{
	REAL smallestincr = 1.e-3;
	REAL dist = 0.;
	REAL incr = 1.0;
	TPZFMatrix solution;
	fFlowCompMesh->LoadSolution(sol0);
	AssembleRhs();
	REAL preverr = Norm(fRhs);
	REAL error;
	dist += incr;
	while (fabs(incr) > smallestincr) {
		solution = sol0;
		solution += (direction *dist);
		fFlowCompMesh->LoadSolution(solution);
		try
		{
			AssembleRhs();
			error = Norm(fRhs);
		}
		catch(TPZOutofRange obj)
		{
			error = 2*preverr+2.;
		}
		if (error < preverr && fabs(dist - 1.0) < 1.e-9) {
			incr = 0.;
		}
		else if (error < preverr && dist < 1.00 && dist > 0.) {
			dist += incr;
			preverr = error;
		}
		else 
		{
			dist -= incr;
			incr *= 0.1;
			dist += incr;
		}
	}
	dist -= incr;
	if (dist <= 0.) {
		cout << "TPBNonLinearAnalysis improper search direction\n";
		direction *= smallestincr;
	}
	else {
		direction *= dist;
		if (fabs(dist - 1.) > smallestincr) {
			// cout << "MaxDescent scale = " << dist << endl;
			// cout.flush();
		}
	}
	fFlowCompMesh->LoadSolution(sol0);
	if(dist > 0.) return 1;
	else return 0;
}



/*!
 \fn TPZEulerAnalysis::CompareRhs()
 */
void TPZEulerAnalysis::CompareRhs()
{
	int iel;
	//int numel = 0;
	int nelem = fCompMesh->NElements();
	TPZElementMatrix ek1(fCompMesh, TPZElementMatrix::EK) ,ef1(fCompMesh, TPZElementMatrix::EF),ef2(fCompMesh, TPZElementMatrix::EF);
	
	TPZAdmChunkVector<TPZCompEl *> &elementvec = fCompMesh->ElementVec();
	TPZFNMatrix<64> diff(8,8);
	//  REAL diffnorm;
	
	for(iel=0; iel < nelem; iel++) {
		TPZCompEl *el = elementvec[iel];
		if(!el) continue;
		el->CalcStiff(ek1,ef1);
		el->CalcResidual(ef2);
		diff =ef1.fMat;
		diff -= ef2.fMat;
		//   diffnorm = Norm(diff);
	}
}


/*!
 \fn TPZEulerAnalysis::SetGMResFront(REAL tol, int numiter, int numvectors)
 */
void TPZEulerAnalysis::SetGMResFront(REAL tol, int numiter, int numvectors)
{
	TPZFrontStructMatrix <TPZFrontNonSym> strfront(Mesh());
	strfront.SetQuiet(1);
	TPZMatrix *front = strfront.CreateAssemble(fRhs,NULL);
	
	TPZStepSolver FrontSolver;
	FrontSolver.SetDirect(ELU);
	FrontSolver.SetMatrix(front);
	
	
	TPZSpStructMatrix StrMatrix(Mesh());
	//TPZFStructMatrix StrMatrix(cmesh);
	SetStructuralMatrix(StrMatrix);
	
	TPZMatrix * mat = StrMatrix.Create();
	TPZStepSolver Solver;
	Solver.SetGMRES(numiter,
					numvectors,
					FrontSolver,
					tol,
					0);
	Solver.SetMatrix(mat);
	SetSolver(Solver);
	fHasFrontalPreconditioner = 0;
	
}


/*!
 \fn TPZEulerAnalysis::SetFrontalSolver()
 */
void TPZEulerAnalysis::SetFrontalSolver()
{
	TPZFrontStructMatrix <TPZFrontNonSym> *StrMatrix = new TPZFrontStructMatrix <TPZFrontNonSym>(Mesh());
	StrMatrix->SetQuiet(1);
	//  TPZMatrix *front = StrMatrix.CreateAssemble(fRhs);
	SetStructuralMatrix(StrMatrix);
	TPZStepSolver FrontSolver;
	FrontSolver.SetDirect(ELU);
	//  FrontSolver.SetMatrix(front);
	SetSolver(FrontSolver);
	fHasFrontalPreconditioner = 1;
}


/*!
 \fn TPZEulerAnalysis::SetGMResBlock(REAL tol, int numiter, int numvec)
 */
void TPZEulerAnalysis::SetGMResBlock(REAL tol, int numiter, int numvec)
{
	TPZSpStructMatrix StrMatrix(Mesh());
	//TPZFStructMatrix StrMatrix(cmesh);
	SetStructuralMatrix(StrMatrix);
	
	TPZMatrix * mat = StrMatrix.Create();
	TPZBlockDiagonalStructMatrix strBlockDiag(Mesh());
	TPZStepSolver Pre;
	TPZBlockDiagonal * block = new TPZBlockDiagonal();//blockDiag.Create();
	strBlockDiag.AssembleBlockDiagonal(*block); // just to initialize structure
	Pre.SetMatrix(block);
	Pre.SetDirect(ELU);
	TPZStepSolver Solver;
	Solver.SetGMRES(numiter,
					numvec,
					Pre,
					tol,
					0);
	Solver.SetMatrix(mat);
	SetSolver(Solver);
	SetBlockDiagonalPrecond(block);
	fHasFrontalPreconditioner = 0;
	
}
