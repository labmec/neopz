/***************************************************************************
 *   Copyright (C) 2006 by Philippe Devloo   *
 *   phil@fec.unicamp.br   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>

#include "tpzdohrsubstruct.h"
#include "tpzdohrmatrix.h"
#include "tpzdohrprecond.h"
#include "pzdohrstructmatrix.h"
#include "pzstepsolver.h"
#include "pzcompel.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "TPZParFrontStructMatrix.h"

#include "pzelast3d.h"
#include "pzbndcond.h"

#include "tpzdohrassembly.h"

#include "pzlog.h"
#include "tpzgensubstruct.h"
#include "pzlog.h"

#include "TPZfTime.h"
#include "TPZTimeTemp.h"
//#include "TPZDataBase.h"
#include "tpzpairstructmatrix.h"

#include "TPZVTKGeoMesh.h"
#include <fstream>
#include <string>

#ifdef LOG4CXX
static LoggerPtr loggerconverge(Logger::getLogger("pz.converge"));
static LoggerPtr logger(Logger::getLogger("main"));
#endif
#include "tpzdohrmatrix.h"

void InsertElasticity(TPZAutoPointer<TPZCompMesh> mesh);
TPZGeoMesh *MalhaPredio();
REAL Height(TPZGeoMesh *gmesh);
int SubStructure(TPZAutoPointer<TPZCompMesh> cmesh, REAL height);

using namespace std;

#include "clock_timer.h"
#include "timing_analysis.h"
#include <sys/resource.h> // getrusage                                                                                                                                                                    
#define TIMING_ANALYSIS

#ifdef LOG4CXX
static LoggerPtr perflog(Logger::getLogger("perf"));
#else
#define LOG4CXX_TRACE(logger, expression)
#define LOG4CXX_DEBUG(logger, expression)
#define LOG4CXX_INFO(logger, expression)
#define LOG4CXX_WARN(logger, expression)
#define LOG4CXX_ERROR(logger, expression)
#define LOG4CXX_FATAL(logger, expression)
#endif

void printrusage2(ostream& out, const char* msg1, const struct rusage& ru1, 
		  const char* msg2, const struct rusage& ru2)
{
  double elapsed1, elapsed2;
  out << "==============================================================" << std::endl;
  out << setw(12) << msg1 << " | " << setw(12) << msg2 << " | " << std::endl;

  elapsed1 = (ru1.ru_utime.tv_sec * 1000.0) + (ru1.ru_utime.tv_usec / 1000.0);
  elapsed2 = (ru2.ru_utime.tv_sec * 1000.0) + (ru2.ru_utime.tv_usec / 1000.0);
  out << setw(12) << setprecision(2) << fixed << elapsed1 << " | "
      << setw(12) << setprecision(2) << fixed << elapsed2 
      << " | ms. user time used" << std::endl;

  elapsed1 = (ru1.ru_stime.tv_sec * 1000.0) + (ru1.ru_stime.tv_usec / 1000.0);
  elapsed2 = (ru2.ru_stime.tv_sec * 1000.0) + (ru2.ru_stime.tv_usec / 1000.0);
  out << setw(12) << setprecision(2) << fixed << elapsed1  << " | "
      << setw(12) << setprecision(2) << fixed << elapsed2 
      << " | ms. system time used" << std::endl;
  out << setw(12) << ru1.ru_maxrss <<   " | " << setw(12) << ru2.ru_maxrss <<   " | integral max resident set size" << std::endl;
  out << setw(12) << ru1.ru_ixrss <<    " | " << setw(12) << ru2.ru_ixrss <<    " | integral shared text memory size" << std::endl;
  out << setw(12) << ru1.ru_idrss <<    " | " << setw(12) << ru2.ru_idrss <<    " | integral unshared data size" << std::endl;
  out << setw(12) << ru1.ru_isrss <<    " | " << setw(12) << ru2.ru_isrss <<    " | integral unshared stack size" << std::endl;
  out << setw(12) << ru1.ru_minflt <<   " | " << setw(12) << ru2.ru_minflt <<   " | page reclaims" << std::endl;
  out << setw(12) << ru1.ru_majflt <<   " | " << setw(12) << ru2.ru_majflt <<   " | page faults" << std::endl;
  out << setw(12) << ru1.ru_nswap <<    " | " << setw(12) << ru2.ru_nswap <<    " | swaps" << std::endl;
  out << setw(12) << ru1.ru_inblock <<  " | " << setw(12) << ru2.ru_inblock <<  " | block input operations" << std::endl;
  out << setw(12) << ru1.ru_oublock <<  " | " << setw(12) << ru2.ru_oublock <<  " | block output operations" << std::endl;
  out << setw(12) << ru1.ru_msgsnd <<   " | " << setw(12) << ru2.ru_msgsnd <<   " | messages sent" << std::endl;
  out << setw(12) << ru1.ru_msgrcv <<   " | " << setw(12) << ru2.ru_msgrcv <<   " | messages received" << std::endl;
  out << setw(12) << ru1.ru_nsignals << " | " << setw(12) << ru2.ru_nsignals << " | signals received" << std::endl;
  out << setw(12) << ru1.ru_nvcsw <<    " | " << setw(12) << ru2.ru_nvcsw <<    " | voluntary context switches" << std::endl;
  out << setw(12) << ru1.ru_nivcsw <<   " | " << setw(12) << ru2.ru_nivcsw <<   " | involuntary context switches" << std::endl;
  out << "==============================================================" << std::endl;
}

int main(int argc, char *argv[])
{

#ifdef TIMING_ANALYSIS
	ClockTimer timer;
	ClockTimer total_timer;
	TimingAnalysis ta;	
	total_timer.start();
#endif

	/* Quando se está usando o tal log4cxx */
	InitializePZLOG("log4cxx.cfg");
	
	//int dim = 2;				// 2 dim	
	//int maxlevel = 6;			// 2 dim
	//int sublevel = 4;			// 2 dim
	
	
	int plevel = 1;
	int numthreads_assemble = -1;
	int numthreads_nsubmesh_assemble = -1;
	int numthreads_decompose = -1;
	int numthreads_multiply = -1;
	int nsubstruct = -1;
	 
	
	plevel = atoi(argv[1]);
	numthreads_nsubmesh_assemble = atoi(argv[2]);
	numthreads_assemble = atoi(argv[3]);
	numthreads_decompose = atoi(argv[4]);
	numthreads_multiply = atoi(argv[5]);
	nsubstruct = atoi(argv[6]);
	 
      
	//std::cout << "Enter PLevel\n";
	//std::cin >> plevel;
	//std::cout << "Enter NThreads\n";
	//std::cin >> numthreads;
	//std::cout << "Enter NSubstructs\n";
	//std::cin >> nsubstruct;
	//std::cout << nsubstruct << std::endl;
	
	std::cout << "Read values for plevel " << plevel << "; nthreads assemble " << numthreads_assemble << 
	"; numsubmesh simultaneous assemble " << numthreads_nsubmesh_assemble << "; nthreads_decompose " << numthreads_decompose << 
	"; numthreads for multiplications " << numthreads_multiply <<" nsubs " << nsubstruct << std::endl;
	//std::cout.flush();
	TPZPairStructMatrix::gNumThreads = numthreads_assemble;		
	// tempo.fnumthreads_nsubmesh_assemble = numthreads_nsubmesh_assemble;		// alimenta timeTemp com o numero de submalhas assembladas simultaneamente
	// tempo.fNumthreads_assemble = numthreads_assemble;							// alimenta timeTemp com o numero de threads para assemblagem
	// tempo.fNumthreads_multiply = numthreads_multiply;								// alimenta timeTemp com o numero de threads para multiplicacao
	// tempo.fNumthreads_decompose = numthreads_decompose;						// alimenta timeTemp com o numero de threads para decomposicao
	tempo.fPolyOrder = plevel;														// alimenta timeTemp com a ordem polinomial
	TPZGeoMesh *gmesh = 0;
	{
#ifdef TIMING_ANALYSIS
		TIME_SEC_BEG_LOG(perflog, timer,"Reading building mesh from ... ");
#endif // TIMING_ANALYSIS                                                                                                    
		
		//TPZGenSubStruct sub(dim,maxlevel,sublevel); // 2 dim
		gmesh = MalhaPredio();

#ifdef TIMING_ANALYSIS
		TIME_SEC_END_LOG(perflog, ta, timer,"Reading building mesh from ... ");
#endif // TIMING_ANALYSIS                                                                                                    

		REAL height = Height(gmesh);

		TPZCompEl::SetgOrder(plevel);
		
		//TPZAutoPointer<TPZCompMesh> cmesh = sub.GenerateMesh(); // 2 dim
		
#ifdef TIMING_ANALYSIS
		TIME_SEC_BEG_LOG(perflog, timer,"Build cmesh");
#endif		
		TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
		InsertElasticity(cmesh);
		cmesh->AutoBuild();
#ifdef TIMING_ANALYSIS
		TIME_SEC_END_LOG(perflog, ta, timer,"Build cmesh");
#endif
		
		
		/*TPZfTime start; 		TPZAnalysis analysis(cmesh);
		TPZParFrontStructMatrix<TPZFrontSym> str(cmesh.operator->());
		str.SetNumberOfThreads(8);
		str.SetNumThreads(4);
		TPZStepSolver step;
		step.SetDirect(ECholesky);
		str.SetNumThreads(4);
		analysis.SetSolver(step);
		analysis.SetStructuralMatrix(str);
		analysis.Run();
		 */
		
		//std::cout << "Tempo de resolucao " << start.ReturnTimeDouble() << std::endl;
		//return 0;
		TIME_SEC_BEG_LOG(perflog, timer,"Building the DohrmanStructMatrix");
		
		TPZDohrStructMatrix dohrstruct(cmesh,numthreads_nsubmesh_assemble,numthreads_decompose);
		
		std::cout << "Number of equations " << cmesh->NEquations() << std::endl;
		std::cout <<  "Number of elements " << cmesh->NElements() << std::endl;
		
		tempo.fNumEq = cmesh->NEquations();												// alimenta timeTemp com o numero de equacoes da malha
		tempo.fNumberofElements = cmesh->NElements();								// alimenta timeTemp com o numero de elementos da malha
		
		dohrstruct.IdentifyExternalConnectIndexes();
		TIME_SEC_END_LOG(perflog,ta,timer,"Building the DohrmanStructMatrix");
		
		std::cout << "Substructuring the mesh\n";
		TPZfTime timetosub; // init of timer

		TIME_SEC_BEG_LOG(perflog, timer,"SubStructure: partition the mesh in submeshes");
		dohrstruct.SubStructure(nsubstruct);					
		TIME_SEC_END_LOG(perflog, ta, timer,"SubStructure: partition the mesh in submeshes");

		//nsubstruct = SubStructure(cmesh, height/8);
		tempo.fNumSub = nsubstruct;						// alimenta timeTemp com o numero de substruturas
		tempo.ft0sub = timetosub.ReturnTimeDouble();		// end of timer
		std::cout << tempo.ft0sub << std::endl;
		
//		sub.SubStructure();
#ifdef LOG4CXX
		{
			std::stringstream str;
			cmesh->Print(str);
			LOGPZ_DEBUG(logger,str.str());
		}
#endif
		
		gmesh = cmesh->Reference();
		
		
		//dohrstruct.SetNumThreads(numthreads_assemble);
		
		TPZAutoPointer<TPZGuiInterface> gui;
		TPZFMatrix rhs(cmesh->NEquations(),1,0.); // New SubStruct;
		TIME_SEC_BEG_LOG(perflog, timer,"CreateAssemble");
		TPZAutoPointer<TPZMatrix> dohr = dohrstruct.CreateAssemble(rhs, gui);
		TIME_SEC_END_LOG(perflog,ta,timer,"CreateAssemble");
		TPZAutoPointer<TPZMatrix> precond = dohrstruct.Preconditioner();
		TPZDohrMatrix<TPZDohrSubstructCondense> *dohrptr = dynamic_cast<TPZDohrMatrix<TPZDohrSubstructCondense> *> (dohr.operator->());
		TPZDohrPrecond<TPZDohrSubstructCondense> *precondptr = dynamic_cast<TPZDohrPrecond<TPZDohrSubstructCondense> *> (precond.operator->());
		dohrptr->SetNumThreads(numthreads_multiply);
		precondptr->SetNumThreads(numthreads_multiply);
		
		
		TPZFMatrix diag(dohr->Rows(),1,5.), produto(dohr->Rows(),1);
		std::cout << "Numero de equacoes " << dohr->Rows() << std::endl;
		tempo.fNumEqCoarse = dohr->Rows();											// alimenta timeTemp com o numero de equacoes coarse

		TIME_SEC_BEG_LOG(perflog, timer,"Multiply");
		dohr->Multiply(diag,produto);
		TIME_SEC_END_LOG(perflog, ta,timer,"Multiply");		
		
		if (!dohrptr) {
			DebugStop();
		}
		dohrptr->AdjustResidual(produto);
		
#ifdef LOG4CXX
		{
			std::stringstream sout;
			produto.Print("O valor do produto", sout );
			diag.Print("O valor da diagonal",sout);
			LOGPZ_DEBUG(loggerconverge,sout.str())
		}
#endif
		diag.Zero();
		TPZStepSolver pre(precond);
		pre.SetMultiply();
		TPZStepSolver cg(dohr);
		//  void SetCG(const int numiterations,const TPZMatrixSolver &pre,const REAL tol,const int FromCurrent);
		
		cg.SetCG(1000,pre,1.e-8,0);
		

		TPZfTime timetosolve; // init of timer
		TIME_SEC_BEG_LOG(perflog, timer,"cg.Solve");
		//EBORIN: O outro substruct utiliza produto em vez de produto                                                                                                                                         
		TIME_SEC_BEG_LOG(perflog, timer,"cg.Solve");
		cg.Solve(produto,diag);
		TIME_SEC_END_LOG(perflog, ta,timer,"cg.Solve");
		tempo.ft6iter = timetosolve.ReturnTimeDouble(); // end of timer
		cout << "Total: " << tempo.ft6iter << std::endl;
		
		cout << "Tempos para multiplicacao: " << tempo.fMultiply << std::endl;
		cout << "Tempos para precondicionamento: " << tempo.fPreCond << std::endl;

		
		string FileName;
		FileName = "Times_in_Line.txt";
		ofstream OutputFile;
		
		bool shouldprint = tempo.NeedsHeader(FileName);			// verify the need of a header
		OutputFile.open(FileName.c_str(), ios::app);					// creates the file
		if (shouldprint == true) tempo.PrintHeader(OutputFile);		// prints the header if It is the first time the program is executed
		
			
		tempo.PrintLine(OutputFile);		// print all the information in one line
		
		//TPZDataBase data;
		//data.Read(FileName);				
		
#ifdef LOG4CXX
		{
			std::stringstream sout;
			diag.Print("Resultado do processo antes do ajuste",sout);
			LOGPZ_INFO(loggerconverge,sout.str())
		}
#endif
		
		dohrptr->AddInternalSolution(diag);
#ifdef LOG4CXX
		{
			std::stringstream sout;
			diag.Print("Resultado do processo iterativo",sout);
			LOGPZ_INFO(loggerconverge,sout.str())
		}
#endif		
	}
	TIME_SEC_BEG_LOG(perflog, timer,"Final steps");
	delete gmesh;
	TIME_SEC_END_LOG(perflog, ta, timer,"Final steps");

	total_timer.stop();
	struct rusage self, children;
	getrusage(RUSAGE_SELF, &self);
	getrusage(RUSAGE_CHILDREN, &children);
	ta.share_report(std::cout, total_timer.getUnits());
	printrusage2(std::cout, "Self", self, "Children", children);

	return EXIT_SUCCESS;
}


int main2(int argc, char *argv[])
{
	/* Quando se está usando o tal log4cxx */
	InitializePZLOG("log4cxx.cfg");
	
	/*
	 TPZFMatrix teste(2,2);
	 TPZFMatrix parte;
	 teste(0,0)=1;
	 teste(0,1)=2;
	 teste(1,0)=3;
	 teste(1,1)=4;
	 teste.GetSub(0,0,2,1,parte);
	 cout << parte << endl;
	 cout << "Hello, world!" << endl;
	 */
	
	
	/**
	 TPZDohrSubstruct meuobjeto;
	 TPZDohrMatrix *matriz = new TPZDohrMatrix();
	 TPZDohrPrecond *precond = new TPZDohrPrecond();
	 TPZStepSolver dohrprecond(precond);
	 dohrprecond.SetMultiply();
	 TPZStepSolver cg(matriz);
	 cg.SetCG(10,dohrprecond,1.e-7,1);
	 TPZFMatrix rhs,result;
	 cg.Solve(rhs,result);*/
	//meuobjeto.
	/*  int dim = 2;
	 TPZGenSubStruct sub(dim,6,3);*/
	//	int dim = 2;
	//	int maxlevel = 6;
	//	int sublevel = 3;
	//	int plevel = 3;
	int dim = 2;
	int maxlevel = 4;
	int sublevel = 1;
	int plevel = 2;
	TPZGenSubStruct sub(dim,maxlevel,sublevel);
	int nk = 8;
	//	int ik;
	//	for(ik=1; ik<nk; ik++)
	//	{
	//		sub.fK[ik] = 1.;//50.*ik;
	//	}
	//sub.fMatDist = TPZGenSubStruct::RandomMat;
	
	TPZCompEl::SetgOrder(plevel);
	
	sub.GenerateMesh();
	
	/*
	 TPZAutoPointer<TPZDohrAssembly> dohrassembly = new TPZDohrAssembly;
	 TPZDohrMatrix<TPZDohrSubstruct> *dohrptr = new TPZDohrMatrix<TPZDohrSubstruct>(dohrassembly);
	 TPZAutoPointer<TPZMatrix> dohr(dohrptr);
	 sub.InitializeDohr(dohr,dohrassembly);
	 
	 
	 // loop over the substructures
	 // This is a lengthy process which should run on the remote processor
	 //	void InitializeMatrices(TPZSubCompMesh *sub, TPZAutoPointer<TPZDohrSubstruct> substruct,  TPZDohrAssembly &dohrassembly);
	 
	 
	 dohrptr->Initialize();
	 #ifdef LOG4CXX
	 {
	 std::stringstream sout;
	 dohrptr->Print("DohrMatrix without condensation", sout);
	 LOGPZ_DEBUG(logger,sout.str())
	 }
	 #endif
	 TPZDohrPrecond<TPZDohrSubstruct> *precondptr = new TPZDohrPrecond<TPZDohrSubstruct>(*dohrptr,dohrassembly);
	 precondptr->Initialize();
	 TPZAutoPointer<TPZMatrix> precond(precondptr);
	 
	 */
	
	TPZAutoPointer<TPZDohrAssembly> dohrassembly2 = new TPZDohrAssembly;
	TPZDohrMatrix<TPZDohrSubstructCondense> *dohrptr2 = new TPZDohrMatrix<TPZDohrSubstructCondense>(dohrassembly2);
	dohrptr2->SetNumThreads(4);
	TPZAutoPointer<TPZMatrix> dohr2(dohrptr2);
	sub.InitializeDohrCondense(dohr2,dohrassembly2);
	dohrptr2->Initialize();
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		dohr2->Print("The dohr matrix condensed",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	
	
	
#ifdef LOG4CXX
	std::stringstream sout;
	sout << "Three dimensional substructures, maxlevel " << maxlevel << " level of substructures " << sublevel << std::endl;
	sout << "Number of substructures " << dohrptr2->SubStructures().size() << std::endl;
	sout << "Interpolation order " << plevel;
	LOGPZ_DEBUG(loggerconverge,sout.str());
#endif
	
	
	TPZDohrPrecond<TPZDohrSubstructCondense> *precondptr2 = new TPZDohrPrecond<TPZDohrSubstructCondense>(*dohrptr2,dohrassembly2);
	precondptr2->Initialize();
	
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Printing after creating the preconditioner\n";
		dohrptr2->Print("After creating the preconditioner",sout);
		LOGPZ_DEBUG(loggerconverge,sout.str());
	}
#endif
	
	
	
	TPZAutoPointer<TPZMatrix> precond2(precondptr2);
	
	
	TPZFMatrix diag(dohr2->Rows(),1,5.), produto(dohr2->Rows(),1), produto2(dohr2->Rows(),1);
	precondptr2->Multiply(diag,produto);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		produto.Print("O valor do produto", sout );
		LOGPZ_DEBUG(loggerconverge,sout.str())
	}
#endif
	precondptr2->Multiply(diag,produto2);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		produto2.Print("O valor do produto2", sout );
		LOGPZ_DEBUG(loggerconverge,sout.str())
	}
#endif
	//#define TOTAL
#ifdef TOTAL
	{
		int dim=dohr->Rows();
		
		TPZFMatrix teste1(dim,dim);
		teste1.Zero();
		TPZFMatrix teste2(dim,dim);
		teste2.Zero();
		
		int i,j;
		TPZFMatrix col(dim,1); //column of the identity matrix
		TPZFMatrix resul(dohr->Rows(),1);
		for (i=0;i<dim;i++) {
			col.Zero();
			col(i,0)=1;
			precondptr->MultAdd(col,col,resul,1,0,0,1);
			for (j=0;j<dim;j++) {
				teste1(i,j) = resul(j,0);
				teste2(i,j) = resul(j,0);
			}
		}
		teste1.Transpose();
		teste1 -= teste2;
		std::cout << "Norma da diferenca das matrizes " << Norm(teste1) << std::endl;
	}
#endif
	std::cout << "Numero de equacoes " << dohr2->Rows() << std::endl;
	//  produto.Print("The value of the product is");
#ifndef MAKEINTERNAL
	diag(0,0) = 0.;
#endif
	dohr2->Multiply(diag,produto);
	dohrptr2->AdjustResidual(produto);
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		produto.Print("O valor do produto", sout );
		diag.Print("O valor da diagonal",sout);
		LOGPZ_DEBUG(loggerconverge,sout.str())
	}
#endif
	diag.Zero();
	TPZStepSolver pre(precond2);
	pre.SetMultiply();
	TPZStepSolver cg(dohr2);
	//  void SetCG(const int numiterations,const TPZMatrixSolver &pre,const REAL tol,const int FromCurrent);
	
	cg.SetCG(500,pre,1.e-8,0);
	cg.Solve(produto,diag);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		diag.Print("Resultado do processo antes do ajuste",sout);
		LOGPZ_INFO(loggerconverge,sout.str())
	}
#endif
	
	dohrptr2->AddInternalSolution(diag);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		diag.Print("Resultado do processo iterativo",sout);
		LOGPZ_INFO(loggerconverge,sout.str())
	}
#endif
	//diag.Print("Resultado do solve");
	/* Solve
	 TPZFMatrix *teste = new TPZFMatrix(2,2);
	 (*teste)(0,0)=1;
	 (*teste)(0,1)=2;
	 (*teste)(1,0)=3;
	 (*teste)(1,1)=4;
	 TPZStepSolver coef;
	 coef.SetMatrix(teste);
	 coef.SetDirect(ELU);
	 TPZFMatrix resul(2,2);
	 resul(0,0)=2;
	 resul(0,1)=3;
	 resul(1,0)=4;
	 resul(1,1)=5;
	 TPZFMatrix res(2,2);
	 coef.Solve(resul,res);
	 cout << res << endl;*/
	
	
	return EXIT_SUCCESS;
}

void InsertElasticity(TPZAutoPointer<TPZCompMesh> mesh)
{
	mesh->SetDimModel(3);
	int nummat = 1;
	REAL E = 1.e6;
	REAL poisson = 0.3;
	TPZManVector<REAL> force(3,0.);
	force[1] = 20.;
	TPZElasticity3D *elast = new TPZElasticity3D(nummat,E,poisson,force);
	TPZAutoPointer<TPZMaterial> elastauto(elast);
	TPZFMatrix val1(3,3,0.),val2(3,1,0.);
	TPZBndCond *bc = elast->CreateBC(elastauto, -1, 0, val1, val2);
	TPZAutoPointer<TPZMaterial> bcauto(bc);
	mesh->InsertMaterialObject(elastauto);
	mesh->InsertMaterialObject(bcauto);
}

TPZGeoMesh *MalhaPredio()
{
	//int nBCs = 1;
	int numnodes=-1;
	int numelements=-1;
	
	string FileName;
	FileName = "8andares02.txt";
	
	{
		bool countnodes = false;
		bool countelements = false;
		
		ifstream read (FileName.c_str());
		
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			if(str == "Coordinates") countnodes = true;
			if(str == "end coordinates") countnodes = false;
			if(countnodes) numnodes++;
			
			if(str == "Elements") countelements = true;
			if(str == "end elements") countelements = false;
			if(countelements) numelements++;
		}
	}
	
	TPZGeoMesh * gMesh = new TPZGeoMesh;
	
	gMesh -> NodeVec().Resize(numnodes);
	
	TPZVec <int> TopolTetra(4);
	
	const int Qnodes = numnodes;
	TPZVec <TPZGeoNode> Node(Qnodes);
	
	//setting nodes coords
	int nodeId = 0, elementId = 0, matElId = 1;
	
	ifstream read;
	read.open(FileName.c_str());
	
	double nodecoordX , nodecoordY , nodecoordZ ;
	
	char buf[1024];
	read.getline(buf, 1024);
	read.getline(buf, 1024);
	std::string str(buf);
	int in;
	for(in=0; in<numnodes; in++)
	{ 
		read >> nodeId;
		read >> nodecoordX;
		read >> nodecoordY;
		read >> nodecoordZ;
		Node[nodeId-1].SetNodeId(nodeId);
		Node[nodeId-1].SetCoord(0,nodecoordX);
		Node[nodeId-1].SetCoord(1,nodecoordY);
		Node[nodeId-1].SetCoord(2,nodecoordZ);
		gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
		
		
	}
	
	{
		
		read.close();
		read.open(FileName.c_str());
		
		
		
		int l , m = numnodes+5;
		for(l=0; l<m; l++)
		{
			read.getline(buf, 1024);
		}
		
		
		int el;
		int matBCid = -1;
		//std::set<int> ncoordz; //jeitoCaju
		for(el=0; el<numelements; el++)
		{
			read >> elementId;
			read >> TopolTetra[0]; //node 1
			read >> TopolTetra[1]; //node 2
			read >> TopolTetra[2]; //node 3
			read >> TopolTetra[3]; //node 4
			
			// O GID comeca com 1 na contagem dos nodes, e nao zero como no PZ, assim o node 1 na verdade é o node 0
			TopolTetra[0]--;
			TopolTetra[1]--;
			TopolTetra[2]--;
			TopolTetra[3]--;
			
			int index;
			TPZGeoEl * tetra = gMesh->CreateGeoElement(ETetraedro, TopolTetra, matElId, index);
			
			// Colocando as condicoes de contorno
			TPZVec <TPZGeoNode> Nodefinder(4);
			TPZVec <REAL> nodecoord(3);
			TPZVec<int> ncoordzVec(0); int sizeOfVec = 0;
			//ncoordz.clear(); //jeitoCaju
			for (int i = 0; i < 4; i++) 
			{
				Nodefinder[i] = gMesh->NodeVec()[TopolTetra[i]];
				Nodefinder[i].GetCoordinates(nodecoord);
				if (nodecoord[2] == 0.)
				{
					//ncoordz.insert(TopolTetra[i]); //jeitoCaju
					sizeOfVec++;
					ncoordzVec.Resize(sizeOfVec);
					ncoordzVec[sizeOfVec-1] = TopolTetra[i];
				}
			}
			//if(ncoordz.size() == 3) //jeitoCaju
			if(ncoordzVec.NElements() == 3)
			{
				/*
				 //jeitoCaju
				 for(int s = tetra->NNodes(); s < tetra->NSides(); s++)
				 {
				 TPZGeoElSide tetraSide(tetra, s);
				 if(tetraSide.NSideNodes() != 3)
				 {
				 continue;
				 }
				 
				 bool ok = true;
				 for(int n = 0; n < tetraSide.NSideNodes(); n++)
				 {
				 int node = tetraSide.SideNodeIndex(n);
				 
				 if(ncoordz.find(node) == ncoordz.end())
				 {
				 ok = false;
				 break;
				 }
				 }
				 if(ok)
				 {
				 TPZGeoElBC(tetraSide,matBCid, *gMesh);						
				 break;
				 }
				 }
				 */
				
				int lado = tetra->WhichSide(ncoordzVec);
				TPZGeoElSide tetraSide(tetra, lado);
				TPZGeoElBC(tetraSide,matBCid); //, *gMesh);		
				//std::cout << "BC #" << nBCs << std::endl;
				//nBCs++;
			}
		}
		
		gMesh->BuildConnectivity();
		
	}
	
	
	
	// identificando as superficies que terao cond de contorno. Coord z dos 3 nos = 0
	//	for (int el = 0; el < numnodes-1; el++) 
	//	{
	//		Nodefind[el] = gMesh->NodeVec()[el];
	//
	//	}
	//	Nodefind.Print(std::cout);
	//	std::cout.flush();
	
	//TPZGeoElBC(TPZGeoEl *el,int side,int matid, TPZGeoMesh &mesh);
	//TPZGeoElBC(TPZGeoElSide &elside,int matid, TPZGeoMesh &mesh);
	
	ofstream arg("malhaPZ.txt");
	gMesh->Print(arg);
	
	return gMesh;
	
}

/*
int SubStructureNathan (TPZAutoPointer<TPZCompMesh> cmesh, REAL height)
{
	int nelem = cmesh->NElements();
	TPZManVector<int> subindex(nelem)
 
}
*/

int SubStructure(TPZAutoPointer<TPZCompMesh> cmesh, REAL height)
{
	int nelem = cmesh->NElements();
	TPZManVector<int> subindex(nelem,-1);
	int iel;
	int nsub = 0;
	for (iel=0; iel<nelem; iel++) 
	{
		TPZCompEl *cel = cmesh->ElementVec()[iel];
		if (!cel) {
			continue;
		}
		TPZGeoEl *gel = cel->Reference();
		if (!gel) {
			continue;
		}
		int nsides = gel->NSides();
		TPZManVector<REAL> center(gel->Dimension(),0.), xco(3,0.);
		gel->CenterPoint(nsides-1,center);
		gel->X(center,xco);
		REAL z = xco[2];
		int floor = (int) z/height;
		nsub = (floor+1) > nsub ? (floor+1) : nsub;
		subindex[iel] = floor;
	}
	
#ifdef DEBUG 
	{
		TPZGeoMesh *gmesh = cmesh->Reference();
		int nelgeo = gmesh->NElements();
		TPZVec<int> domaincolor(nelgeo,-999);
		int cel;
		int nel = cmesh->NElements();
		for (cel=0; cel<nel; cel++) {
			TPZCompEl *compel = cmesh->ElementVec()[cel];
			if(!compel) continue;
			TPZGeoEl *gel = compel->Reference();
			if (!gel) {
				continue;
			}
			domaincolor[gel->Index()] = subindex[cel];
		}
		ofstream vtkfile("partition.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, domaincolor);
	}
#endif
	
	int isub;
	TPZManVector<TPZSubCompMesh *> submeshes(nsub,0);
	for (isub=0; isub<nsub; isub++) 
	{
		int index;
		std::cout << '^'; std::cout.flush();
		submeshes[isub] = new TPZSubCompMesh(cmesh,index);
		
		if (index < subindex.NElements()) 
		{
			subindex[index] = -1;
		}
	}
	for (iel=0; iel<nelem; iel++) 
	{
		int domindex = subindex[iel];
		if (domindex >= 0) 
		{
			TPZCompEl *cel = cmesh->ElementVec()[iel];
			if (!cel) 
			{
				continue;
			}
			submeshes[domindex]->TransferElement(cmesh.operator->(),iel);
		}
	}
	cmesh->ComputeNodElCon();
	for (isub=0; isub<nsub; isub++) 
	{
		submeshes[isub]->MakeAllInternal();
		std::cout << '*'; std::cout.flush();
	}
	
	cmesh->ComputeNodElCon();
	cmesh->CleanUpUnconnectedNodes();
	return nsub;
	
}

REAL Height(TPZGeoMesh *gmesh)
{
	TPZAdmChunkVector<TPZGeoNode> &nodevec = gmesh->NodeVec();
	int nnodes = nodevec.NElements();
	int in;
	REAL maxz = 0.;
	for (in=0; in<nnodes; in++) {
		REAL z = nodevec[in].Coord(2);
		maxz = (maxz < z) ? z : maxz;
	}
	return maxz;
}

