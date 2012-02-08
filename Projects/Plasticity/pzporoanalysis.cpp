//$Id: pzporoanalysis.cpp,v 1.7 2009-10-05 03:49:58 erick Exp $
#include "pzporoanalysis.h"

#include "pzlog.h"
#include "pzmattemporal.h"
#include "pzsolve.h" 

#ifdef LOG4CXX
static LoggerPtr PoroAnalysisLogger(Logger::getLogger("analysis.poro"));
#endif

using namespace std;


TPZPoroElastoPlasticAnalysis::TPZPoroElastoPlasticAnalysis() : 
     TPZElastoPlasticAnalysis(), fPorousMaterialIds(), fRhsLast() 
{
	FindPorousMaterials();
}

TPZPoroElastoPlasticAnalysis::TPZPoroElastoPlasticAnalysis(TPZCompMesh *mesh,std::ostream &out)
  : TPZElastoPlasticAnalysis(mesh, out), fPorousMaterialIds(), fRhsLast()
{	
	int numeq = fCompMesh->NEquations();
	fRhsLast.Redim(numeq,1);
	
	FindPorousMaterials();
}

TPZPoroElastoPlasticAnalysis::~TPZPoroElastoPlasticAnalysis()
{
#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << "<<< TPZPoroElastoPlasticAnalysis::~TPZPoroElastoPlasticAnalysis() *** Killing Object\n";
   LOGPZ_INFO(PoroAnalysisLogger,sout.str().c_str());
}
#endif
}

/*
REAL TPZPoroElastoPlasticAnalysis::LocalAssemble(int precond)
{		
	TPZElastoPlasticAnalysis::LocalAssemble(precond);

	REAL norm = Norm(fRhs);
	
	#ifdef LOG4CXX
	{
	   std::stringstream sout;
	   sout << "<<< TPZPoroElastoPlasticAnalysis::LocalAssemble() *** "
	        << " with Norm(Rhs) = " << norm;
	   LOGPZ_INFO(PoroAnalysisLogger,sout.str().c_str());
	}
	#endif
	

int i;
	for(i = 0; i < 4; i++)cout << "\t" << fSolution(i,0);
	cout << endl;
	for(i = 0; i < 4; i++)cout << "\t" << fRhs(i,0);
	cout << endl;
	
	return norm;
}
*/
void TPZPoroElastoPlasticAnalysis::SetAdvancedState()
{
   SetContributionTime(Advanced_CT);
   fCompMesh->LoadSolution(fSolution);
}

void TPZPoroElastoPlasticAnalysis::SetLastState()
{
   SetContributionTime(Last_CT);
   fCompMesh->LoadSolution(fSolution);
}

void TPZPoroElastoPlasticAnalysis::SetContributionTime(TPZContributeTime time)
{
	int i, n = fPorousMaterialIds.NElements();
	
	for(i = 0; i < n; i++)
	{
		TPZAutoPointer<TPZMaterial> pMat = fCompMesh->FindMaterial(fPorousMaterialIds[i]);
		TPZMatTemporal * pMatTemp = dynamic_cast<TPZMatTemporal *>(pMat.operator->());
		if(pMatTemp)pMatTemp->SetContributionTime(time);
	}
}

void TPZPoroElastoPlasticAnalysis::SetDeltaT(const REAL deltaT)
{
	int i, n = fPorousMaterialIds.NElements();
	
	for(i = 0; i < n; i++)
	{
		TPZAutoPointer<TPZMaterial> pMat = fCompMesh->FindMaterial(fPorousMaterialIds[i]);
		TPZMatTemporal * pMatTemp = dynamic_cast<TPZMatTemporal *>(pMat.operator->());
		if(pMatTemp)pMatTemp->SetDeltaT(deltaT);
	}
}
		
int::TPZPoroElastoPlasticAnalysis::FindPorousMaterials()
{
	fPorousMaterialIds.Resize(0);
	int i, n = fCompMesh->NMaterials();
	
	for(i = 0; i < n; i++)
	{
		TPZAutoPointer<TPZMaterial> pMat = fCompMesh->MaterialVec()[i];
		TPZMatTemporal * pMatTemp = dynamic_cast<TPZMatTemporal *>(pMat.operator->());
		if(pMatTemp)fPorousMaterialIds.Push(pMat->Id() );
	}
	
	return fPorousMaterialIds.NElements();
}

void TPZPoroElastoPlasticAnalysis::Run(std::ostream &out,REAL tol ,int numiter,
									TPZPostProcAnalysis * ppAnalysis, int res)
{
#ifdef LOG4CXX
{
   std::stringstream sout;
   sout << ">>> TPZPoroElastoPlasticAnalysis::Run() ***";
   LOGPZ_INFO(PoroAnalysisLogger,sout.str().c_str());
}
#endif

	IterativeProcess(out, tol, numiter);
		
	#ifdef LOG4CXX
	{
	   std::stringstream sout;
	   sout << "*** TPZPoroElastoPlasticAnalysis::Run() *** IterativeProcess() completed";
	   LOGPZ_INFO(PoroAnalysisLogger,sout.str().c_str());
	}
	#endif

	AcceptSolution();
		
	if(ppAnalysis)
	{
		#ifdef LOG4CXX
		{
		   std::stringstream sout;
		   sout << "*** TPZPoroElastoPlasticAnalysis::Run() *** PostProcessing ";
		   LOGPZ_INFO(PoroAnalysisLogger,sout.str().c_str());
		}
		#endif
		TransferSolution(*ppAnalysis);
		ppAnalysis->PostProcess(res);
	}
		
	#ifdef LOG4CXX
	{
	   std::stringstream sout;
	   sout << "<<< TPZPoroElastoPlasticAnalysis::ManageIterativeProcess() *** Exiting";
	   LOGPZ_INFO(PoroAnalysisLogger,sout.str().c_str());
	}
	#endif
}


REAL TPZPoroElastoPlasticAnalysis::AcceptSolution(const int ResetOutputDisplacements)
{	
	
	int i, n = fPorousMaterialIds.NElements();
	int nstate = 0;
	
	for(i = 0; i < n; i++)
	{
		TPZAutoPointer<TPZMaterial> pMat = fCompMesh->FindMaterial(fPorousMaterialIds[i]);
		if(pMat->NStateVariables() > nstate)
			nstate = pMat->NStateVariables();
	}	
	
	n = fSolution.Rows();

	if(ResetOutputDisplacements)
	{
		fCumSol.Zero();
	}else{
		fCumSol += fSolution; 
	}
	
	#ifdef LOG4CXX
	{
	   std::stringstream sout;
	   sout << ">>> TPZPoroElastoPlasticAnalysis::AcceptSolution *** "
	        << " with Norm(fCumSol) = " << Norm(fCumSol);
	   LOGPZ_INFO(PoroAnalysisLogger,sout.str().c_str());
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
