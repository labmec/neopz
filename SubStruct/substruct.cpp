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

#include "tpzdohrassembly.h"

#include "pzlog.h"
#include "tpzgensubstruct.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr loggerconverge(Logger::getLogger("pz.converge"));
static LoggerPtr logger(Logger::getLogger("main"));
#endif

using namespace std;

int main(int argc, char *argv[])
{
	/* Quando se está usando o tal log4cxx */
	InitializePZLOG("log4cxx.cfg");
	
	int dim = 2;
	int maxlevel = 8;
	int sublevel = 4;
	int plevel = 3;
	int numthreads = 3;
	TPZGeoMesh *gmesh = 0;
	{
		TPZGenSubStruct sub(dim,maxlevel,sublevel);

		TPZCompEl::SetgOrder(plevel);
		
		TPZAutoPointer<TPZCompMesh> cmesh = sub.GenerateMesh();
		gmesh = cmesh->Reference();
		
		TPZDohrStructMatrix dohrstruct(cmesh);
		dohrstruct.SetNumThreads(numthreads);
		
		TPZAutoPointer<TPZGuiInterface> gui;
		TPZFMatrix rhs;
		TPZAutoPointer<TPZMatrix> dohr = dohrstruct.CreateAssemble(rhs, gui);
		TPZAutoPointer<TPZMatrix> precond = dohrstruct.Preconditioner();
		
		
		TPZFMatrix diag(dohr->Rows(),1,5.), produto(dohr->Rows(),1);
		std::cout << "Numero de equacoes " << dohr->Rows() << std::endl;
		dohr->Multiply(diag,produto);
		
		TPZDohrMatrix<TPZDohrSubstructCondense> *dohrptr = dynamic_cast<TPZDohrMatrix<TPZDohrSubstructCondense> *> (dohr.operator->());
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
		
		cg.SetCG(100,pre,1.e-8,0);
		cg.Solve(produto,diag);
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
	delete gmesh;
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
