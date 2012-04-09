/**
 * @file
 * @brief Contains the implementation of the TPZDohrPrecond methods. 
 */
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
#include "tpzdohrprecond.h"
//#include <iostream>
//#include <cstdlib>
#include "tpzdohrsubstructCondense.h"

#include "pzskylmat.h"

#include "pzvisualmatrix.h"

#include "tpzdohrassemblelist.h"

#include <sstream>
#include "pzlog.h"

#include "TPZfTime.h"
#include "TPZTimeTemp.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("substruct.dohrprecond"));
static LoggerPtr loggerv1v2(Logger::getLogger("substruct.v1v2"));
#endif

using namespace std;

template<class TSubStruct>
TPZDohrPrecond<TSubStruct>::TPZDohrPrecond(TPZDohrMatrix<TSubStruct> &origin, TPZAutoPointer<TPZDohrAssembly> assemble)
: TPZMatrix<REAL>(origin), fGlobal(origin.SubStructures()), fCoarse(0), fNumCoarse(origin.NumCoarse()), fNumThreads(0), fAssemble(assemble)
{
	fNumThreads = origin.NumThreads();
	//  Initialize();
}

template<class TSubStruct>
TPZDohrPrecond<TSubStruct>::TPZDohrPrecond(const TPZDohrPrecond<TSubStruct> &cp) : TPZMatrix<REAL>(cp), fGlobal(cp.fGlobal), fCoarse(0), 
fNumCoarse(cp.fNumCoarse), fNumThreads(cp.fNumThreads), fAssemble(cp.fAssemble) 
{
	if (cp.fCoarse) {
		fCoarse = (TPZStepSolver<REAL> *) cp.fCoarse->Clone();
	}
}
template<class TSubStruct>
TPZDohrPrecond<TSubStruct>::~TPZDohrPrecond()
{
	if (fCoarse) 
	{
		delete fCoarse;
		fCoarse = 0;
	}
}

template<class TSubStruct>
void TPZDohrPrecond<TSubStruct>::MultAdd(const TPZFMatrix<REAL> &x,const TPZFMatrix<REAL> &y, TPZFMatrix<REAL> &z, const REAL alpha,const REAL beta,const int opt,const int stride) const {
	if ((!opt && Cols() != x.Rows()*stride) || Rows() != x.Rows()*stride)
		Error( "Operator* <matrices with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		Error ("TPZFMatrix::MultiplyAdd incompatible dimensions\n");
	}
	TPZfTime precondi; // init of timer
	int rows = Rows();
	int cols = Cols();
	int c;
	PrepareZ(y,z,beta,opt,stride);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		x.Print("x entry vector",sout);
		LOGPZ_DEBUG(logger,sout.str())
		LOGPZ_DEBUG(loggerv1v2,sout.str())
	}
#endif
	TPZFMatrix<REAL> v1(rows,x.Cols(),0.);
	TPZFMatrix<REAL> v2(cols,x.Cols(),0.);
	if(fNumThreads <= 0)
	{
		ComputeV1(x,v1);
		ComputeV2(x,v2);
	}
	else 
	{
		TPZVec<pthread_t> AllThreads(fNumThreads+2);
		TPZDohrPrecondThreadV1Data<TSubStruct> v1threaddata(this,x,v1);
		
		pthread_create(&AllThreads[0], 0, TPZDohrPrecondThreadV1Data<TSubStruct>::ComputeV1, &v1threaddata);
		//		TPZDohrPrecondThreadV1Data<TSubStruct>::ComputeV1(&v1threaddata);
		
		TPZAutoPointer<TPZDohrAssembleList> assemblelist = new TPZDohrAssembleList(fGlobal.size(),v2,this->fAssemble);
		
		
		TPZDohrPrecondV2SubDataList<TSubStruct> v2work(assemblelist);
		typename std::list<TPZAutoPointer<TSubStruct> >::const_iterator it;
		
		int isub=0;
		//Criar tarefa que execute a distribuicao de cada elemento do fGlobal
		for(it= fGlobal.begin(); it != fGlobal.end(); it++,isub++)
		{
			TPZFMatrix<REAL> *Residual_local = new TPZFMatrix<REAL>;
			fAssemble->Extract(isub,x,*Residual_local);
			TPZDohrPrecondV2SubData<TSubStruct> data(isub,*it,Residual_local);
			v2work.AddItem(data);
		}
		
		int i;
		for (i=0; i<fNumThreads; i++) {
			pthread_create(&AllThreads[i+2], 0, TPZDohrPrecondV2SubDataList<TSubStruct>::ThreadWork, &v2work);
		}
		//		v2work.ThreadWork(&v2work);
		
		pthread_create(&AllThreads[1], 0, TPZDohrAssembleList::Assemble, assemblelist.operator->());
		//		assemblelist->Assemble(assemblelist.operator->());
		
		for (i=0; i<fNumThreads+2; i++) {
			void *result;
			pthread_join(AllThreads[i], &result);
		}
		//		ComputeV2(x,v2);
	}
	v2 += v1;
	
#ifndef MAKEINTERNAL	
	isub=0;
	//Criar tarefa que execute a distribuicao de cada elemento do fGlobal
	for(it= fGlobal.begin(); it != fGlobal.end(); it++,isub++)
	{
		TPZFNMatrix<100> v2Expand((*it)->fNEquations,1,0.), v3Expand((*it)->fNEquations,1,0.);
		int neqs = (*it)->fGlobalEqs.NElements();
		TPZFMatrix<REAL> v3_local(neqs,1,0.), v2_local(neqs,1,0.);
		fAssemble->Extract(isub,v2,v2_local);
		int i;
		for (i=0;i<neqs;i++) 
		{
			std::pair<int,int> ind = (*it)->fGlobalEqs[i];
			v2Expand(ind.first,0) += v2_local(i,0);
		}
		
		(*it)->Contribute_v3_local(v2Expand,v3Expand);
		for (i=0;i<neqs;i++) 
		{
			std::pair<int,int> ind = (*it)->fGlobalEqs[i];
			v3_local(i,0) += v3Expand(ind.first,0);
		}
#ifdef LOG4CXX
		{
			std::stringstream sout;
			v2Expand.Print("v1+v2 Expand",sout);
			v3Expand.Print("v3 Expand", sout);
			v2_local.Print("v1+v2 local",sout);
			v3_local.Print("v3 local",sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		fAssemble->Assemble(isub,v3_local,v2);
	}		
#endif
	// wait task para finalizacao da chamada
	// esperar a versao correta do v1
	/* Sum v1+v2+v3 with z */
    int xcols = x.Cols();
    for (int ic=0; ic<xcols; ic++) 
    {
        for (c=0; c<rows; c++) {
            z(c,ic) += v2(c,ic);
        }
    }
	tempo.fPreCond.Push(precondi.ReturnTimeDouble()); // end of timer
}


template<class TSubStruct>
void TPZDohrPrecond<TSubStruct>::Initialize()
{
	//Compute the skyline of the coarse equations
	TPZManVector<int> skyline(fNumCoarse);
	int ic;
	for (ic=0; ic<fNumCoarse; ic++) {
		skyline[ic] = ic;
	}
	int nsub = fAssemble->fCoarseEqs.NElements();
	int isub;
	for (isub=0; isub<nsub; isub++) {
		int nc = fAssemble->fCoarseEqs[isub].NElements();
		int ic;
		int mineq = fAssemble->fCoarseEqs[isub][0];
		for (ic=0; ic<nc; ic++) {
			int eq = fAssemble->fCoarseEqs[isub][ic];
			mineq = mineq > eq ? eq : mineq;
		}
		for (ic=0; ic<nc; ic++) {
			int eq = fAssemble->fCoarseEqs[isub][ic];
			if(skyline[eq] > mineq) skyline[eq] = mineq;
		}
	}
	/* Computing K(c) */
	TPZMatrix<REAL> *coarse = new TPZSkylMatrix<REAL>(fNumCoarse,skyline);
#ifdef DEBUG
	{
		TPZFMatrix<REAL> coarse2(*coarse);
		for (isub=0; isub<nsub; isub++) {
			int nc = fAssemble->fCoarseEqs[isub].NElements();
			int ic;
			for (ic=0; ic<nc; ic++) {
				int ieq = fAssemble->fCoarseEqs[isub][ic];
				int jc;
				for (jc=0; jc<nc; jc++) {
					int jeq = fAssemble->fCoarseEqs[isub][jc];
					coarse2(ieq,jeq) = 1.;
				}
			}
			
		}
		VisualMatrix(coarse2,"CoarseMatrix.vtk");
	}
#endif
	typename std::list<TPZAutoPointer<TSubStruct> >::iterator it;
	int count = 0;
	for(it= fGlobal.begin(); it != fGlobal.end(); it++,count++) 
	{
		(*it)->Contribute_Kc(*coarse,fAssemble->fCoarseEqs[count]);
	}
	fCoarse = new TPZStepSolver<REAL>(coarse);
}

template<class TSubStruct>
void TPZDohrPrecond<TSubStruct>::ComputeV1(const TPZFMatrix<REAL> &x, TPZFMatrix<REAL> &v1) const
{
	/* Computing r(c) */
	TPZFMatrix<REAL> CoarseResidual(fNumCoarse,x.Cols());
	CoarseResidual.Zero();
	typename std::list<TPZAutoPointer<TSubStruct> >::const_iterator it;
	
	int isub = 0;
	for(it= fGlobal.begin(); it != fGlobal.end(); it++, isub++) {
		TPZFMatrix<REAL> xloc, CoarseResidual_local;
		fAssemble->Extract(isub,x,xloc);
		//		(*it)->LoadWeightedResidual(xloc);
		(*it)->Contribute_rc_local(xloc,CoarseResidual_local);
		fAssemble->AssembleCoarse(isub,CoarseResidual_local,CoarseResidual);
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		CoarseResidual.Print("Coarse Residual",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	/* Computing K(c)_inverted*r(c) and stores it in "product" */
	fCoarse->SetDirect(ELDLt);
	//Dado global 
	TPZFMatrix<REAL> CoarseSolution(fNumCoarse,x.Cols());
	fCoarse->Solve(CoarseResidual,CoarseSolution);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		CoarseSolution.Print("CoarseSolution",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	isub=0;
	//Criar tarefa que execute a distribuicao de cada elemento do fGlobal
	for(it= fGlobal.begin(); it != fGlobal.end(); it++,isub++)
	{
		// Gerenciamento Global->Local sobre o product
		//product é administrado pelo DM mas permanece no processador 0
		// tarefa separada, expansao da solucao coarse
		TPZFMatrix<REAL> v1_local,CoarseSolution_local;
		fAssemble->ExtractCoarse(isub,CoarseSolution,CoarseSolution_local);
		(*it)->Contribute_v1_local(v1_local,CoarseSolution_local);
		
		fAssemble->Assemble(isub,v1_local,v1);
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		v1.Print("v1 vector",sout);
		LOGPZ_DEBUG(loggerv1v2,sout.str())
	}
#endif	
}

template<class TSubStruct>
void TPZDohrPrecond<TSubStruct>::ComputeV2(const TPZFMatrix<REAL> &x, TPZFMatrix<REAL> &v2) const
{
	
	typename std::list<TPZAutoPointer<TSubStruct> >::const_iterator it;
	
	int isub=0;
	//Criar tarefa que execute a distribuicao de cada elemento do fGlobal
	for(it= fGlobal.begin(); it != fGlobal.end(); it++,isub++)
	{
		// contribute v2 deve ser uma tarefa inicializada mais cedo
		TPZFNMatrix<100> Residual_local,v2_local;
		fAssemble->Extract(isub,x,Residual_local);
		(*it)->Contribute_v2_local(Residual_local,v2_local);
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Substructure " << isub << std::endl;
			Residual_local.Print("Residual local",sout);
            v2_local.Print("v2_local",sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		//		v2_local += v1_local;
		fAssemble->Assemble(isub,v2_local,v2);
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		v2.Print("v2 vector",sout);
		LOGPZ_DEBUG(loggerv1v2,sout.str())
	}
#endif
}		

template<class TSubStruct>
void *TPZDohrPrecondV2SubDataList<TSubStruct>::ThreadWork(void *voidptr)
{
	TPZDohrPrecondV2SubDataList<TSubStruct>	*myptr = (TPZDohrPrecondV2SubDataList<TSubStruct> *) voidptr;
	TPZDohrPrecondV2SubData<TSubStruct> data = myptr->PopItem();
	while (data.IsValid()) {
		data.fSubStructure->Contribute_v2_local(data.fInput_local,data.fv2_local->fAssembleData);
		myptr->fAssemblyStructure->AddItem(data.fv2_local);
		data = myptr->PopItem();
	}
	return voidptr;
}


template class TPZDohrPrecond<TPZDohrSubstruct>;
template class TPZDohrPrecond<TPZDohrSubstructCondense>;