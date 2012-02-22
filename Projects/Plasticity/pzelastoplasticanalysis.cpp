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
static LoggerPtr EPAnalysisLogger(Logger::getLogger("analysis.elastoplastic"));
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

void TPZElastoPlasticAnalysis::IterativeProcess(std::ostream &out,REAL tol,int numiter)
{

#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << ">>> TPZElastoPlasticAnalysis::IterativeProcess() ***"
        << "\nEntering method with parameters:"
	    << "\n tol = " << tol
		<< "\n numiter = " << numiter;
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif
	
   int iter = 0;
   REAL normRhs, normDelU = 0;
   int numeq = fCompMesh->NEquations();

   fSolution.Zero();

   TPZAnalysis::LoadSolution();

   TPZFMatrix prevSol(fSolution);
   if(prevSol.Rows() != numeq) prevSol.Redim(numeq,1);

   fRhs.Zero();

   normRhs = this->LocalAssemble(0/*precondNormRhs*/);

   UpdatePrecond();

#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << "*** TPZElastoPlasticAnalysis::IterativeProcess() *** "
        << " Start of newton method with Norm(Rhs) = " << normRhs;
   //sout << Rhs();
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

   do{
	  iter++;

	  normDelU = this->LocalSolve();

	  //* Performs an ZAXPY operation being *this += alpha * p
   	  //* @param alpha Being alpha on above opereation
      //* @param p Being p on above operation
      // ZAXPY(const REAL alpha, const TPZFMatrix &p);

/* 
	  REAL alpha = 1 / log(normRhs);
	  if(alpha > 1.)alpha = 1.;
	   
	  cout << "\n Alpha = " << alpha;
	   
	  prevsol.ZAXPY(alpha, fSolution);
	  fSolution = prevsol;
	   
*/
	  fSolution += prevSol; 
	  prevSol = fSolution;
	   
	  TPZAnalysis::LoadSolution();
	   
	  fRhs.Zero();
	   
	  normRhs = this->LocalAssemble(0/*precondNormRhs*/);
	   
	   cout << "ITERACAO DE NEWTON = "<< iter <<endl;
	   cout << "NORMA VETOR = "<< normRhs <<endl;
	   cout << "NORMA normDelU = "<< normDelU <<endl;
	   ofstream arg1("convVec.txt");
	   ofstream arg2("convU.txt");
	   arg1 << iter <<"\t"<<normRhs;
	   arg2 << iter <<"\t"<<normDelU;
	   
	  UpdatePrecond();

#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << "*** TPZElastoPlasticAnalysis::IterativeProcess() *** "
        << " End of " << iter << "-th newton step with Norm(Rhs) = " << normRhs
		<< " and Norm(DeltaU) = " << normDelU;
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

   }while(normRhs > tol && iter < numiter);

#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << "<<< TPZElastoPlasticAnalysis::IterativeProcess() *** ";
   if(normRhs > tol)
	{
		sout << " #### Truncating Method with normRhs = " << normRhs << " after " << iter << " Newton Steps";
	}else{
		sout << " Exiting Converged Method after " << iter << " Newton Steps with normRhs = " << normRhs;
	}
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

   return;

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
	
    TPZMatrix * pMatrix = TPZAnalysis::fSolver->Matrix().operator->();
	
    if(precond && pMatrix)
	{
		TPZFMatrix localRhs = fRhs;
		int i;
		for (i = 0; i < size; i++)localRhs(i,0) /= pMatrix->operator()(i,i);
		norm = Norm(localRhs);
	}
	
	else
	{
		norm = Norm(fRhs);
	}
//	ofstream fileout("rigidez.nb");
//	fSolver->Matrix()->Print("Rigidez = ", fileout, EMathematicaInput);
//	#ifdef LOG4CXX
//	{
//	   std::stringstream sout;
//	   sout << "<<< TPZElastoPlasticAnalysis::LocalAssemble() *** "
//	        << " with Norm(Rhs) = " << norm;
//	   sout << "\n Rhs =\n" << TPZAnalysis::Rhs()
//			<< "\n Matrix =\n";
//	   if( TPZAnalysis::fSolver )
//		{
//			TPZAnalysis::fSolver->Matrix()->Print(0, sout);
//		}else
//		{
//			sout << " none";
//		}
//	   LOGPZ_DEBUG(EPAnalysisLogger,sout.str().c_str());
//	}
//	#endif

    return norm;
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
	
	int i;
//	cout << "\nLocalSolve: fSolution=\n";
//	for(i = 0; i < 4; i++)cout << "\t" << fSolution(i);

    TPZAnalysis::Solve();

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
	
	std::map<int, TPZAutoPointer<TPZMaterial > > & refMatVec = fCompMesh->MaterialVec();

    std::map<int, TPZAutoPointer<TPZMaterial > >::iterator mit;
	
	TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem; // defined in file pzelastoplastic.h
	TPZMatWithMem<TPZPoroElastoPlasticMem> * pMatWithMem2; // define in file pzporous.h

    for(mit=refMatVec.begin(); mit!= refMatVec.end(); mit++){
        pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>( mit->second.operator->() );
		if(pMatWithMem != NULL) pMatWithMem->SetUpdateMem(update);
        pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZPoroElastoPlasticMem> *>( mit->second.operator->() );
		if(pMatWithMem2 != NULL) pMatWithMem2->SetUpdateMem(update);
    }
	
}

REAL TPZElastoPlasticAnalysis::AcceptSolution(const int ResetOutputDisplacements)
{	
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
	
	REAL norm = this->LocalAssemble(0);
	
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

   TPZFMatrix rangeMatrix(numeq, 1, range);
   
   TPZVec<REAL> coefs(1,1.);

   CheckConvergence(*this,fSolution,rangeMatrix,coefs);
   
}

void TPZElastoPlasticAnalysis::ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &coefs, int icase){

	int neq = fCompMesh->NEquations();
	tangent.Redim(neq,neq);
	TPZFMatrix rhs(neq,1);
	TPZFStructMatrix substitute(Mesh());
	TPZAutoPointer<TPZGuiInterface> guiInterface(0);
	substitute.Assemble(tangent,rhs,guiInterface);
//	TPZStructMatrix::Assemble(tangent, rhs, *Mesh());
}

int TPZElastoPlasticAnalysis::NumCases(){
	return 1;
}

void TPZElastoPlasticAnalysis::Residual(TPZFMatrix &residual, int icase){
	int neq = fCompMesh->NEquations();
	TPZFMatrix tangent(neq,neq);
	residual.Redim(neq,1);
	TPZFStructMatrix substitute(Mesh());
	TPZAutoPointer<TPZGuiInterface> guiInterface(0);
	substitute.Assemble(residual,guiInterface);
//	TPZStructMatrix::Assemble(/*tangent,*/ residual, *Mesh());
	residual *= -1;
}

void TPZElastoPlasticAnalysis::SetPrecond(TPZMatrixSolver &precond){
  if(fPrecond) delete fPrecond;
    fPrecond = (TPZMatrixSolver *) precond.Clone();
}

void TPZElastoPlasticAnalysis::UpdatePrecond()
{
   if(fPrecond)
   {
		TPZMatrix * pMatrix = TPZAnalysis::fSolver->Matrix().operator->();
		TPZMatrix * pPrecondMat = fPrecond->Matrix().operator->();
		pPrecondMat->Zero();
		TPZBlockDiagonal *pBlock = dynamic_cast<TPZBlockDiagonal *>(pPrecondMat);
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
	TPZMatrix * mat = StrMatrix.Create();

    TPZBlockDiagonalStructMatrix strBlockDiag(Mesh());
    TPZStepSolver Pre;
    TPZBlockDiagonal * block = new TPZBlockDiagonal();
	
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
    TPZStepSolver Solver;
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
	TPZMatrix * mat = StrMatrix.Create();

    TPZBlockDiagonalStructMatrix strBlockDiag(Mesh());
    TPZStepSolver Pre;
    TPZBlockDiagonal * block = new TPZBlockDiagonal();
	
#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << "*** TPZElastoPlasticAnalysis::SetBiCGStab_Jacobi() *** Assembling Block Diagonal Preconditioning matrix\n";
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

    strBlockDiag.AssembleBlockDiagonal(*block); // just to initialize structure
	Pre.SetMatrix(block);
    //Pre.SetDirect(ELU);
	Pre.SetJacobi(numiter, tol, 0);
    TPZStepSolver Solver;
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

    TPZMatrix * mat = StrMatrix.Create();

    TPZStepSolver Solver;
    //Solver.SetDirect(ELU);// ECholesky -> simétrica e positiva definida
	Solver.SetDirect(ELU);
    Solver.SetMatrix(mat);

    this->SetSolver(Solver);
}

void TPZElastoPlasticAnalysis::TransferSolution(TPZPostProcAnalysis & ppanalysis)
{
	TPZFMatrix bkpSolution = fSolution;
	
	fSolution = fCumSol;
	
	TPZAnalysis::LoadSolution();
	
	ppanalysis.TransferSolution();
	
	fSolution = bkpSolution;
	
	TPZAnalysis::LoadSolution();	
}

void TPZElastoPlasticAnalysis::ManageIterativeProcess(std::ostream &out,REAL tol,int numiter,
									int BCId, int nsteps, REAL PGRatio,
									TPZFMatrix & val1Begin, TPZFMatrix & val1End,
									TPZFMatrix & val2Begin, TPZFMatrix & val2End,
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

	TPZAutoPointer<TPZMaterial> mat = fCompMesh->FindMaterial(BCId);
	TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat.operator->());
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
		
		IterativeProcess(out, tol, numiter);
		
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

TPZCompEl * TPZElastoPlasticAnalysis::CreateCubeElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeCube > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateLinearElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeLinear > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreatePointElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapePoint > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreatePrismElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapePrism > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreatePyramElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapePiram > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateQuadElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeQuad > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateTetraElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeTetra > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateTriangElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int &index)
{
	return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeTriang > >(mesh,gel,index);
}
