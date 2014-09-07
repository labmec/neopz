//$Id: pzelastoplasticanalysis.cpp,v 1.27 2010-11-23 18:58:05 diogo Exp $
#include "pzelastoplasticanalysis.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "checkconv.h"
#include "pzstrmatrix.h"
#include "pzelastoplastic.h"
#include "tpzautopointer.h"
#include "pzcompelwithmem.h"
#include "pzelastoplasticmem.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzbdstrmatrix.h"
#include "pzstepsolver.h"
#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelastoplastic2D.h"

#include <map>
#include <set>
#include <stdio.h>
#include <fstream>

#include "pzsolve.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr EPAnalysisLogger(Logger::getLogger("pz.analysis.elastoplastic"));
static LoggerPtr loggertest(Logger::getLogger("testing"));
#endif

using namespace std;


TPZElastoPlasticAnalysis::TPZElastoPlasticAnalysis() : TPZNonLinearAnalysis(), fPrecond(NULL) {
	//Mesh()->Solution().Zero(); already performed in the nonlinearanalysis base class
	//fSolution.Zero();
}

TPZElastoPlasticAnalysis::TPZElastoPlasticAnalysis(TPZCompMesh *mesh,std::ostream &out) : TPZNonLinearAnalysis(mesh,out), fPrecond(NULL) {

	int numeq = fCompMesh->NEquations();
	fCumSol.Redim(numeq,1);
	fCumSol.Zero();
	fSolution.Redim(numeq,1);
	fSolution.Zero();
	
	TPZAnalysis::LoadSolution();
}

TPZElastoPlasticAnalysis::~TPZElastoPlasticAnalysis()
{
	if(fPrecond)delete fPrecond;
	
#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << "<<< TPZElastoPlasticAnalysis::~TPZElastoPlasticAnalysis() *** Killing Object\n";
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif
}

REAL TPZElastoPlasticAnalysis::LineSearch(const TPZFMatrix<REAL> &Wn, const TPZFMatrix<REAL> &DeltaW, TPZFMatrix<REAL> &NextW, REAL RhsNormPrev, REAL &RhsNormResult, int niter, bool & converging){

    TPZFMatrix<REAL> Interval = DeltaW;

#ifdef DEBUG
    {
        LoadSolution(Wn);
        AssembleResidual();
        AdjustResidual(fRhs);
        STATE normprev = Norm(fRhs);
        if (fabs(normprev - RhsNormPrev) > 1.e-6) {
            std::stringstream sout;
            sout << "Norm of Wn " << Norm(Wn) << std::endl;
            sout << "Input previous norm " << RhsNormPrev << " Computed Norm " << normprev;
            LOGPZ_ERROR(EPAnalysisLogger, sout.str())
        }
    }
#endif
    REAL scalefactor = 1.;
    int iter = 0;
    do {
        Interval *= scalefactor;
        NextW = Wn;
        NextW += Interval;
        LoadSolution(NextW);
        AssembleResidual();
        AdjustResidual(fRhs);
        RhsNormResult = Norm(fRhs);
#ifndef PLASTICITY_CLEAN_OUT
        std::cout << "Scale factor " << scalefactor << " resnorm " << RhsNormResult << std::endl;
#endif
        scalefactor *= 0.5;
        iter++;
    } while (RhsNormResult > RhsNormPrev && iter < 30);
    if(fabs(RhsNormResult - RhsNormPrev)<1.e-6 )
    {
        converging=false;
    }
    else
    {
        converging=true;
    }
    scalefactor *= 2.;
	return scalefactor;
	
}//void

/// Iterative process using the linear elastic material as tangent matrix
void TPZElastoPlasticAnalysis::IterativeProcess(std::ostream &out, TPZAutoPointer<TPZMatrix<STATE> > linearmatrix, REAL tol, int numiter, bool linesearch)
{
	int iter = 0;
	REAL error = 1.e10;
	int numeq = fCompMesh->NEquations();
	Mesh()->Solution().Zero();
	fSolution.Zero();
	
	TPZFMatrix<REAL> prevsol(fSolution);
	if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
    
    if (linearmatrix) {
        std::cout << __FILE__ << ":" << __LINE__ << "Decomposed " << linearmatrix->IsDecomposed() << std::endl;
    }
    
#ifdef LOG4CXX
    if(EPAnalysisLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Solution norm of fSolution " << Norm(fSolution);
        LOGPZ_DEBUG(EPAnalysisLogger, sout.str())
    }
#endif
	
    TPZAnalysis::AssembleResidual();
    AdjustResidual(fRhs);
    REAL RhsNormPrev = Norm(fRhs);
    
//    std::cout << __LINE__ << " Norm prevsol " << Norm(prevsol) << std::endl;
	bool linesearchconv=true;
	while(error > tol && iter < numiter) {
		
		//fSolution.Redim(0,0);
        REAL RhsNormResult = 0.;
#ifdef DEBUG
        if (linearmatrix) {
            std::cout << __FILE__ << ":" << __LINE__ << "Decomposed " << linearmatrix->IsDecomposed() << std::endl;
        }
#endif
		LocalSolve();
        TPZFMatrix<STATE> solkeep(fSolution);
//        std::cout << __LINE__ << " Norm prevsol " << Norm(prevsol) << std::endl;
		if (linesearch){
            {
//                std::cout << __LINE__ << " Norm prevsol " << Norm(prevsol) << std::endl;
                TPZFMatrix<STATE> nextsol(prevsol);
                nextsol += solkeep;
                LoadSolution(nextsol);
                AssembleResidual();
                AdjustResidual(fRhs);
                RhsNormResult = Norm(fRhs);
//                std::cout << __LINE__ << " Norm prevsol " << Norm(prevsol) << std::endl;
//                std::cout << "RhsNormResult " << RhsNormResult << std::endl;
            }
            if (RhsNormResult > tol && RhsNormResult > RhsNormPrev) {
                fSolution = prevsol;
                TPZFMatrix<REAL> nextSol;
//                std::cout << __LINE__ << " Norm prevsol " << Norm(prevsol) << std::endl;
#ifdef LOG4CXX
                if (EPAnalysisLogger->isDebugEnabled()) {
                    std::stringstream sout;
                    std::cout << __LINE__ << " Norm prevsol " << Norm(prevsol) << std::endl;
                    LOGPZ_DEBUG(EPAnalysisLogger, sout.str())
                }
#endif
                const int niter = 10;
                this->LineSearch(prevsol, solkeep, nextSol, RhsNormPrev, RhsNormResult, niter,linesearchconv);
                fSolution = nextSol;
            }
            else
            {
//                fSolution = solkeep;
//                fSolution += prevsol;
            }
		}
		else{
			fSolution += prevsol;
            LoadSolution(fSolution);
            AssembleResidual();
            AdjustResidual(fRhs);
            RhsNormResult = Norm(fRhs);
		}
		
		prevsol -= fSolution;
//        std::cout << __LINE__ << " Norm prevsol " << Norm(prevsol) << std::endl;
		REAL normDeltaSol = Norm(prevsol);
		prevsol = fSolution;
//        std::cout << __LINE__ << " Norm prevsol " << Norm(prevsol) << std::endl;
		REAL norm = RhsNormResult;
        
#ifdef DEBUG
        {
            LoadSolution(fSolution);
            AssembleResidual();
            AdjustResidual(fRhs);
            REAL rhsnorm = Norm(fRhs);
            std::cout << "Norm rhs reported " << RhsNormResult << " now computed " << rhsnorm << std::endl;
        }
#endif
        
        RhsNormPrev = RhsNormResult;
		//       out << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << endl;
        std::cout << "Iteracao n : " << (iter+1) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << RhsNormResult << endl;
        //        std::cout << "Iteracao n : " << (iter+1) << " : fRhs : " << fRhs << endl;
		
		if(norm < tol) {
            std::cout << "\nTolerancia atingida na iteracao : " << (iter+1) << endl;
            std::cout << "\n\nNorma da solucao |Delta(Un)|  : " << norm << endl << endl;
			
		} else
			if( (norm - error) > 1.e-9 ) {
                std::cout << "\nDivergent Method\n";
			}
		error = norm;
		iter++;
		out.flush();
	}
    
}


void TPZElastoPlasticAnalysis::IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv,bool &ConvOrDiverg) {
	
	int iter = 0;
	REAL error = 1.e10;
	int numeq = fCompMesh->NEquations();
	//Mesh()->Solution().Zero();
	//fSolution->Zero();
    

	
	TPZFMatrix<REAL> prevsol(fSolution);
	if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
    
#ifdef LOG4CXX_keep
    {
        std::stringstream sout;
        fSolution.Print("Solution for checkconv",sout);
        LOGPZ_DEBUG(EPAnalysisLogger, sout.str())
    }
#endif
	
	if(checkconv){
		TPZVec<REAL> coefs(1,1.);
		TPZFMatrix<REAL> range(numeq,1,1.e-5);
		CheckConvergence(*this,fSolution,range,coefs);
	}
    
    REAL RhsNormPrev = LocalAssemble(0);
	bool linesearchconv=true;
    
	while(error > tol && iter < numiter) {
        
        if(iter!=0)
        {
            LocalAssemble(0);
        }
		
		fSolution.Redim(0,0);
        REAL RhsNormResult = 0.;
		LocalSolve();
		if (linesearch){
			TPZFMatrix<REAL> nextSol;
			const int niter = 10;
			this->LineSearch(prevsol, fSolution, nextSol, RhsNormPrev, RhsNormResult, niter,linesearchconv);
			fSolution = nextSol;
		}
		else{
			fSolution += prevsol;
            LoadSolution(fSolution);
            AssembleResidual();
            AdjustResidual(fRhs);
            RhsNormResult = Norm(fRhs);
		}
		
		prevsol -= fSolution;
		REAL normDeltaSol = Norm(prevsol);
		prevsol = fSolution;
		REAL norm = RhsNormResult;
        RhsNormPrev = RhsNormResult;
		//       out << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << endl;
        std::cout << "Iteracao n : " << (iter+1) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << RhsNormResult << endl;
//        std::cout << "Iteracao n : " << (iter+1) << " : fRhs : " << fRhs << endl;
		
		if(norm < tol) {
			 std::cout << "\nTolerancia atingida na iteracao : " << (iter+1) << endl;
			 std::cout << "\n\nNorma da solucao |Delta(Un)|  : " << norm << endl << endl;
            ConvOrDiverg=true;
			
		} else
			if( (norm - error) > 1.e-9  || linesearchconv ==false) {
                std::cout << "\nDivergent Method -- Exiting Consistent Tangent Iterative Process \n";
                std::cout << "\n Trying linearMatrix IterativeProcess \n\n";
                ConvOrDiverg=false;
                return;
			}
		error = norm;
		iter++;
		out.flush();
	}
    
}

void TPZElastoPlasticAnalysis::IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv) {
	
	int iter = 0;
	REAL error = 1.e10;
	int numeq = fCompMesh->NEquations();
	//Mesh()->Solution().Zero();
	//fSolution->Zero();
    
    
	
	TPZFMatrix<REAL> prevsol(fSolution);
	if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
    
#ifdef LOG4CXX_keep
    {
        std::stringstream sout;
        fSolution.Print("Solution for checkconv",sout);
        LOGPZ_DEBUG(EPAnalysisLogger, sout.str())
    }
#endif
	
	if(checkconv){
		TPZVec<REAL> coefs(1,1.);
		TPZFMatrix<REAL> range(numeq,1,1.e-5);
		CheckConvergence(*this,fSolution,range,coefs);
	}
    
    REAL RhsNormPrev = LocalAssemble(0);
	bool linesearchconv=true;
	while(error > tol && iter < numiter) {
		
		fSolution.Redim(0,0);
        REAL RhsNormResult = 0.;
		LocalSolve();
		if (linesearch){
			TPZFMatrix<REAL> nextSol;
			const int niter = 10;
			this->LineSearch(prevsol, fSolution, nextSol, RhsNormPrev, RhsNormResult, niter,linesearchconv);
			fSolution = nextSol;
		}
		else{
			fSolution += prevsol;
            LoadSolution(fSolution);
            AssembleResidual();
            AdjustResidual(fRhs);
            RhsNormResult = Norm(fRhs);
		}
		
		prevsol -= fSolution;
		REAL normDeltaSol = Norm(prevsol);
		prevsol = fSolution;
		REAL norm = RhsNormResult;
        RhsNormPrev = RhsNormResult;
		//       out << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << endl;
        std::cout << "Iteracao n : " << (iter+1) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << RhsNormResult << endl;
        //        std::cout << "Iteracao n : " << (iter+1) << " : fRhs : " << fRhs << endl;
		
		if(norm < tol) {
            std::cout << "\nTolerancia atingida na iteracao : " << (iter+1) << endl;
            std::cout << "\n\nNorma da solucao |Delta(Un)|  : " << norm << endl << endl;
			
		} else
			if( (norm - error) > 1.e-9 ) {
                std::cout << "\nDivergent Method \n";

			}
		error = norm;
		iter++;
		out.flush();
	}
    
}



REAL TPZElastoPlasticAnalysis::LocalAssemble(int precond)
{
	int size = fSolution.Rows();
	
//	#ifdef LOG4CXX
//	{
//	   std::stringstream sout;
//	   sout << ">>> TPZElastoPlasticAnalysis::LocalAssemble() *** "
//	        << "About to ASSEMBLE the ek and ef stiffness matrix and load vector (rank = " << size << ")";
//	   LOGPZ_DEBUG(EPAnalysisLogger,sout.str().c_str());
//	}
//	#endif
	
    TPZAnalysis::Assemble();
	
	REAL norm = 0;
	
    TPZMatrix<REAL> * pMatrix = TPZAnalysis::fSolver->Matrix().operator->();
	
    AdjustTangentMatrix(*pMatrix);
    AdjustResidual(fRhs);
    
    if(precond && pMatrix)
	{
		TPZFMatrix<REAL> localRhs = fRhs;
		int i;
		for (i = 0; i < size; i++)localRhs(i,0) /= pMatrix->operator()(i,i);
		norm = Norm(localRhs);
	}
	
	else
	{
		norm = Norm(fRhs);
	}
    
#ifdef LOG4CXX
    if (loggertest->isDebugEnabled()) {
        std::stringstream sout;
        pMatrix->Print("Global Matrix",sout);
        fRhs.Print("Rhs", sout);
        LOGPZ_DEBUG(loggertest, sout.str())
    }
#endif

    return norm;
}

/// Apply zero on the lines and columns of the Dirichlet boundary conditions
void TPZElastoPlasticAnalysis::AdjustTangentMatrix(TPZMatrix<STATE> &matrix)
{
    std::set<long>::iterator it;
    long size = matrix.Rows();
    for (it = fEquationstoZero.begin(); it != fEquationstoZero.end(); it++) {
        int eq = *it;
        for (int i=0; i<size; i++) {
            matrix.Put(eq, i, 0.);
            matrix.Put(i, eq, 0.);
        }
        matrix.Put(eq, eq, 1.);
    }
    
}

/// Apply zero to the equations of the Dirichlet boundary conditions
void TPZElastoPlasticAnalysis::AdjustResidual(TPZFMatrix<STATE> &rhs)
{
    std::set<long>::iterator it;
    long size = rhs.Rows();
    long cols = rhs.Cols();
    for (it = fEquationstoZero.begin(); it != fEquationstoZero.end(); it++) {
        long eq = *it;
        for (long i=0; i<size; i++) {
            for (long j=0; j<cols; j++) 
            {
                rhs.Put(eq, j, 0.);
            }
        }
    }    
}



REAL TPZElastoPlasticAnalysis::LocalSolve()
{
	#ifdef LOG4CXX
	{
	   std::stringstream sout;
	   sout << ">>> TPZElastoPlasticAnalysis::Solve() *** "
	        << "About to SOLVE the linear system";
	   LOGPZ_DEBUG(EPAnalysisLogger,sout.str().c_str());
	}
	#endif
	
    std::set<long>::iterator it;
    for (it = fEquationstoZero.begin(); it != fEquationstoZero.end(); it++) {
        fRhs(*it,0) = 0.;
    }
//	cout << "\nLocalSolve: fSolution=\n";
//	for(i = 0; i < 4; i++)cout << "\t" << fSolution(i);

#ifdef LOG4CXX_K
    {
        std::stringstream sout;
        fRhs.Print("Right hand side",sout);
        LOGPZ_DEBUG(EPAnalysisLogger, sout.str())
    }
#endif
    TPZAnalysis::Solve();
#ifdef LOG4CXX_K
    {
        std::stringstream sout;
        fSolution.Print("Solution",sout);
        LOGPZ_DEBUG(EPAnalysisLogger, sout.str())
    }
#endif

//    cout << "\n DisplacementSIZE = "<< fSolution << endl; 
  //  cout << "\n Displacement = "<< fSolution << endl; 
    
    REAL norm = Norm(fSolution);
	#ifdef LOG4CXX
	{
	   std::stringstream sout;
	   sout << "<<< TPZElastoPlasticAnalysis::Solve() *** "
	        << " with Norm(DeltaU) = " << norm;
	//   sout << Rhs();
	   LOGPZ_DEBUG(EPAnalysisLogger,sout.str().c_str());
	}
	#endif

    return norm;
}

void TPZElastoPlasticAnalysis::SetUpdateMem(int update)
{
	if(!fCompMesh)return;
	
	std::map<int, TPZMaterial *> & refMatVec = fCompMesh->MaterialVec();

    std::map<int, TPZMaterial * >::iterator mit;
	
	TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem; // defined in file pzelastoplastic.h
	TPZMatWithMem<TPZPoroElastoPlasticMem> * pMatWithMem2; // define in file pzporous.h

    for(mit=refMatVec.begin(); mit!= refMatVec.end(); mit++)
    {
        pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>( mit->second );
		if(pMatWithMem != NULL)
        {
           pMatWithMem->SetUpdateMem(update); 
        }
        pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZPoroElastoPlasticMem> *>( mit->second);
		if(pMatWithMem2 != NULL)
        {
            pMatWithMem2->SetUpdateMem(update);
        }
    }
	
}

#include "pzelasmat.h"

REAL TPZElastoPlasticAnalysis::AcceptSolution(const int ResetOutputDisplacements)
{	
    
    TPZMaterial *mat = fCompMesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }
    TPZElasticityMaterial *elasmat = dynamic_cast<TPZElasticityMaterial *>(mat);
    if(elasmat)
    {
        // the material is linear
        return 0.;
    }
    
    
	if(ResetOutputDisplacements)
	{
		fCumSol.Zero();
	}else{
		fCumSol += fSolution;
	}

	#ifdef LOG4CXX
	{
	   std::stringstream sout;
	   sout << ">>> TTPZElastoPlasticAnalysis::AcceptSolution *** "
	        << " with Norm(fCumSol) = " << Norm(fCumSol);
	   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
	}
	#endif
	
	this->SetUpdateMem(true);
	
	fRhs.Zero();
	
    AssembleResidual();
    AdjustResidual(fRhs);
	REAL norm = Norm(fRhs);
	
	this->SetUpdateMem(false);
	
	fSolution.Zero();
	
	TPZAnalysis::LoadSolution();
    
	
	return norm;
}

void TPZElastoPlasticAnalysis::CheckConv(std::ostream &out, REAL range) {

#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << ">>> TPZElastoPlasticAnalysis::CheckConv() ***"
        << "\nEntering method with parameters:"
	    << "\n range = " << range;
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif
	
   int numeq = fCompMesh->NEquations();

   TPZFMatrix<REAL> rangeMatrix(numeq, 1, range);
   
   TPZVec<REAL> coefs(1,1.);

   CheckConvergence(*this,fSolution,rangeMatrix,coefs);
   
}

void TPZElastoPlasticAnalysis::ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &coefs, int icase){

	int neq = fCompMesh->NEquations();
	tangent.Redim(neq,neq);
	TPZFMatrix<REAL> rhs(neq,1);
	TPZFStructMatrix substitute(Mesh());
	TPZAutoPointer<TPZGuiInterface> guiInterface(0);
	substitute.Assemble(tangent,rhs,guiInterface);
//	TPZStructMatrix::Assemble(tangent, rhs, *Mesh());
}

int TPZElastoPlasticAnalysis::NumCases(){
	return 1;
}

void TPZElastoPlasticAnalysis::Residual(TPZFMatrix<REAL> &residual, int icase){
	int neq = fCompMesh->NEquations();
//	TPZFMatrix<REAL> tangent(neq,neq);
	residual.Redim(neq,1);
	TPZFStructMatrix substitute(Mesh());
	TPZAutoPointer<TPZGuiInterface> guiInterface(0);
	substitute.Assemble(residual,guiInterface);
//	TPZStructMatrix::Assemble(/*tangent,*/ residual, *Mesh());
	residual *= -1;
}

void TPZElastoPlasticAnalysis::SetPrecond(TPZMatrixSolver<REAL> &precond){
  if(fPrecond) delete fPrecond;
    fPrecond = (TPZMatrixSolver<REAL> *) precond.Clone();
}

void TPZElastoPlasticAnalysis::UpdatePrecond()
{
   if(fPrecond)
   {
		TPZMatrix<REAL> * pMatrix = TPZAnalysis::fSolver->Matrix().operator->();
		TPZMatrix<REAL> * pPrecondMat = fPrecond->Matrix().operator->();
		pPrecondMat->Zero();
		TPZBlockDiagonal<REAL> *pBlock = dynamic_cast<TPZBlockDiagonal<REAL> *>(pPrecondMat);
		pBlock->BuildFromMatrix(*pMatrix);
   }
}

void TPZElastoPlasticAnalysis::SetBiCGStab(int numiter, REAL tol)
{
#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << ">>> TPZElastoPlasticAnalysis::SetBiCGStab() *** numiter = " << numiter << " and tol=" << tol;
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

	TPZSpStructMatrix StrMatrix(Mesh());
    this->SetStructuralMatrix(StrMatrix);
	TPZMatrix<REAL> * mat = StrMatrix.Create();

    TPZBlockDiagonalStructMatrix strBlockDiag(Mesh());
    TPZStepSolver<REAL> Pre;
    TPZBlockDiagonal<REAL> * block = new TPZBlockDiagonal<REAL>();
	
#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << "*** TPZElastoPlasticAnalysis::SetBiCGStab() *** Assembling Block Diagonal Preconditioning matrix\n";
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

    strBlockDiag.AssembleBlockDiagonal(*block); // just to initialize structure
	Pre.SetMatrix(block);
    Pre.SetDirect(ELU);
    TPZStepSolver<REAL> Solver;
 	Solver.SetBiCGStab(numiter, Pre, tol, 0);
    Solver.SetMatrix(mat);
    this->SetSolver(Solver);
	this->SetPrecond(Pre);

#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << "<<< TPZElastoPlasticAnalysis::SetBiCGStab() *** Exiting\n";
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

}


void TPZElastoPlasticAnalysis::SetBiCGStab_Jacobi(int numiter, REAL tol)
{
#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << ">>> TPZElastoPlasticAnalysis::SetBiCGStab_Jacobi() *** numiter = " << numiter << " and tol=" << tol;
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

	TPZSpStructMatrix StrMatrix(Mesh());
//	TPZFStructMatrix StrMatrix(Mesh());
    this->SetStructuralMatrix(StrMatrix);
	TPZMatrix<REAL> * mat = StrMatrix.Create();

    TPZBlockDiagonalStructMatrix strBlockDiag(Mesh());
    TPZStepSolver<REAL> Pre;
    TPZBlockDiagonal<REAL> * block = new TPZBlockDiagonal<REAL>();
	
#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << "*** TPZElastoPlasticAnalysis::SetBiCGStab_Jacobi() *** Assembling Block Diagonal Preconditioning matrix\n";
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

    strBlockDiag.AssembleBlockDiagonal(*block); // just to initialize structure
	Pre.SetMatrix(block);
    //    Pre.SetDirect(ELU);
    //Pre.SetDirect(ELDLt);
	Pre.SetJacobi(numiter, tol, 0);
    TPZStepSolver<REAL> Solver;
 	Solver.SetBiCGStab(numiter, Pre, tol, 0);
    Solver.SetMatrix(mat);
    this->SetSolver(Solver);
	this->SetPrecond(Pre);

#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << "<<< TPZElastoPlasticAnalysis::SetBiCGStab_Jacobi() *** Exiting\n";
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif
}

void TPZElastoPlasticAnalysis::SetLU()
{
#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << ">>> TPZElastoPlasticAnalysis::SetLU() ***\n";
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif
	
    TPZFStructMatrix StrMatrix(Mesh());
    this->SetStructuralMatrix(StrMatrix);

    TPZMatrix<REAL> * mat = StrMatrix.Create();

    TPZStepSolver<REAL> Solver;
    //Solver.SetDirect(ELU);// ECholesky -> simétrica e positiva definida
	Solver.SetDirect(ELU);
    Solver.SetMatrix(mat);

    this->SetSolver(Solver);
}

void TPZElastoPlasticAnalysis::TransferSolution(TPZPostProcAnalysis & ppanalysis)
{
	TPZFMatrix<REAL> bkpSolution = fSolution;
   
	
	fSolution = fCumSol;
//	 fSolution.Print();
	TPZAnalysis::LoadSolution();//Carrega a solucao convergida no analysis
	//passa o cum sol para o post
	ppanalysis.TransferSolution();//Transfere solucao convergida para o pos processamento
	
    
	fSolution = bkpSolution;
	
	TPZAnalysis::LoadSolution();	
}

void TPZElastoPlasticAnalysis::ManageIterativeProcess(std::ostream &out,REAL tol,int numiter,
									int BCId, int nsteps, REAL PGRatio,
									TPZFMatrix<REAL> & val1Begin, TPZFMatrix<REAL> & val1End,
									TPZFMatrix<REAL> & val2Begin, TPZFMatrix<REAL> & val2End,
									TPZPostProcAnalysis * ppAnalysis, int res)
{
										
	if(!fCompMesh)return;
										
#ifdef LOG4CXX
{
	
   std::stringstream sout;
   sout << "<<< TPZElastoPlasticAnalysis::ManageIterativeProcess() ***";
   sout << "\nWith parameters:\n";
   sout << "\ntol = " << tol;
   sout << "\nnumiter = " << numiter;
   sout << "\nBCId = " << BCId;
   sout << "\nnsteps = " << nsteps;
   sout << "\nPGRatio = " << PGRatio;
   sout << "\nval1Begin = " << val1Begin;
   sout << "\nval1End = " << val1End;
   sout << "\nval2Begin = " << val2Begin;
   sout << "\nval2End = " << val2End;
   if(ppAnalysis)
	{
		sout << "\nppanalysis set";
	}else
	{
		sout << "\nppanalysis NOT set";
	}
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif
										
	// computing the initial value for the PG progression such that its sum equals one;
	REAL a0;
	if(fabs(PGRatio - 1.) < 1.e-3)
	{
	    a0 = 1. / REAL(nsteps);
	}else{
		a0 = (PGRatio - 1) / (pow(PGRatio,nsteps) - 1.);
	}
	TPZFNMatrix<36> val1(6,6,0.), deltaVal1(6,6,0.);
	TPZFNMatrix< 6> val2(6,1,0.), deltaVal2(6,1,0.);
	
	deltaVal1 = val1End;
	deltaVal1.ZAXPY(-1., val1Begin);
	deltaVal2 = val2End;
	deltaVal2.ZAXPY(-1., val2Begin);
	
	// ZAXPY operation: *this += alpha * p			

	TPZMaterial * mat = fCompMesh->FindMaterial(BCId);
	TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat);
	if(!pBC)return;
 
    int i;
	for(i = 0; i < nsteps; i++)
	{
		REAL stepLen;
		if(fabs(PGRatio - 1.) < 1.e-3)
		{
			stepLen = REAL(i+1) / REAL(nsteps);
		}else{
		    stepLen = a0 * (pow(PGRatio,i+1) - 1) / (PGRatio - 1.);
		}
		
		val1 = val1Begin;
		val1.ZAXPY(stepLen, deltaVal1);
		val2 = val2Begin;
		val2.ZAXPY(stepLen, deltaVal2);
		
		pBC->Val1() = val1;
		pBC->Val2() = val2;
		
		#ifdef LOG4CXX
		{
		   std::stringstream sout;
		   sout << "*** TPZElastoPlasticAnalysis::ManageIterativeProcess() *** load step " << i;
		   sout << " stepLen = " << stepLen;
		   sout << "\nBC.val1() = " << val1;
		   sout << "\nBC.val2() = " << val2;
		   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
		}
		#endif
		
        bool linesearch = false;
        bool checkconv = false;
            bool convordiv;
		IterativeProcess(out, tol, numiter, linesearch, checkconv,convordiv);
        
		
		#ifdef LOG4CXX
		{
		   std::stringstream sout;
		   sout << "*** TPZElastoPlasticAnalysis::ManageIterativeProcess() *** load step " << i << " ended";
		   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
		}
		#endif
		
		AcceptSolution();
		
		if(ppAnalysis)
		{
			#ifdef LOG4CXX
			{
			   std::stringstream sout;
			   sout << "*** TPZElastoPlasticAnalysis::ManageIterativeProcess() *** PostProcessing ";
			   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
			}
			#endif
			TransferSolution(*ppAnalysis);
			ppAnalysis->PostProcess(res);
		}
	}
		
	#ifdef LOG4CXX
	{
	   std::stringstream sout;
	   sout << "<<< TPZElastoPlasticAnalysis::ManageIterativeProcess() *** Exiting";
	   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
	}
	#endif
}

// CompEl create Functions setup

#include "pzintel.h"

#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "tpzpoint.h"

#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "tpzline.h"

#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "tpztriangle.h"

#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"

#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "tpzprism.h"

#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "tpztetrahedron.h"

#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "tpzpyramid.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "tpzcube.h"

#include "pzelctemp.h"


void TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(TPZCompMesh *cmesh)
{
/*	pzgeom::TPZGeoPoint::fp = TPZElastoPlasticAnalysis::CreatePointElWithMem;
	 pzgeom::TPZGeoQuad::fp = TPZElastoPlasticAnalysis::CreateQuadElWithMem;
	pzgeom::TPZGeoTriangle::fp = TPZElastoPlasticAnalysis::CreateTriangElWithMem;
	pzgeom::TPZGeoPrism::fp = TPZElastoPlasticAnalysis::CreatePrismElWithMem;
	pzgeom::TPZGeoTetrahedra::fp = TPZElastoPlasticAnalysis::CreateTetraElWithMem;
	pzgeom::TPZGeoPyramid::fp = TPZElastoPlasticAnalysis::CreatePyramElWithMem;
	pzgeom::TPZGeoCube::fp = TPZElastoPlasticAnalysis::CreateCubeElWithMem;
*/
    TPZManVector<TCreateFunction,10> functions(8);
    functions[EPoint] = &TPZElastoPlasticAnalysis::CreatePointElWithMem;
	functions[EOned] = TPZElastoPlasticAnalysis::CreateLinearElWithMem;
	functions[EQuadrilateral] = TPZElastoPlasticAnalysis::CreateQuadElWithMem;
	functions[ETriangle] = TPZElastoPlasticAnalysis::CreateTriangElWithMem;
	functions[EPrisma] = TPZElastoPlasticAnalysis::CreatePrismElWithMem;
	functions[ETetraedro] = TPZElastoPlasticAnalysis::CreateTetraElWithMem;
	functions[EPiramide] = TPZElastoPlasticAnalysis::CreatePyramElWithMem;
	functions[ECube] = TPZElastoPlasticAnalysis::CreateCubeElWithMem;
    cmesh->ApproxSpace().SetCreateFunctions(functions);

}

TPZCompEl * TPZElastoPlasticAnalysis::CreateCubeElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, long &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeCube > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateLinearElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, long &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeLinear > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreatePointElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, long &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapePoint > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreatePrismElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, long &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapePrism > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreatePyramElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, long &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapePiram > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateQuadElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, long &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeQuad > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateTetraElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, long &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeTetra > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateTriangElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, long &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeTriang > >(mesh,gel,index);
}

void TPZElastoPlasticAnalysis::IdentifyEquationsToZero()
{
    fEquationstoZero.clear();
    long nel = fCompMesh->NElements();
    for (long iel=0; iel<nel; iel++) {
        TPZCompEl *cel = fCompMesh->ElementVec()[iel];
        if (!cel) {
            continue;
        }
        TPZMaterial *mat = cel->Material();
        if (!mat) {
            continue;
        }
        int matid = mat->Id();
        if (fMaterialIds.find(matid) == fMaterialIds.end()) {
            continue;
        }
        int direction = fMaterialIds[matid];
        long nc = cel->NConnects();
        for (long ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            long seqnum = c.SequenceNumber();
            long pos = fCompMesh->Block().Position(seqnum);
            int blsize = fCompMesh->Block().Size(seqnum);
            for (long i=pos+direction; i<pos+blsize; i+=2) {
                fEquationstoZero.insert(i);
            }
        }
    }
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Equations to zero ";
        std::set<long>::iterator it;
        for (it=fEquationstoZero.begin(); it!= fEquationstoZero.end(); it++) {
            sout << *it << " ";
        }
        LOGPZ_DEBUG(EPAnalysisLogger, sout.str())
    }
#endif
}