/**
 * @file
 * @brief Contains the implementation of the TPZDohrMatrix methods. 
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
#include "tpzdohrmatrix.h"
//#include "tpzdohrsubstruct.h"

#include "tpzdohrassembly.h"
#include "pzlog.h"

#include "TPZfTime.h"
#include "TPZTimeTemp.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("substruct.dohrsubstruct"));
#endif


template<class TSubStruct>
TPZDohrMatrix<TSubStruct>::TPZDohrMatrix(TPZAutoPointer<TPZDohrAssembly> assembly)
: TPZMatrix(), fNumThreads(0), fAssembly(assembly)
{
}


template<class TSubStruct>
TPZDohrMatrix<TSubStruct>::~TPZDohrMatrix()
{
}


template<class TSubStruct>
void TPZDohrMatrix<TSubStruct>::MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
										const REAL alpha,const REAL beta,const int opt,const int stride) const
{
	TPZfTime mult;
	if ((!opt && Cols() != x.Rows()*stride) || Rows() != x.Rows()*stride)
		Error( "Operator* <matrixs with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		Error ("TPZFMatrix::MultiplyAdd incompatible dimensions\n");
	}
	PrepareZ(y,z,beta,opt,stride);
	
	typename SubsList::const_iterator iter;
	int isub = 0;
	if (fNumThreads == 0) {
		for (iter=fGlobal.begin();iter!=fGlobal.end();iter++,isub++) {
			TPZFMatrix xlocal,zlocal;
			fAssembly->Extract(isub,x,xlocal);
			zlocal.Redim(xlocal.Rows(),xlocal.Cols());
			(*iter)->ContributeKULocal(alpha,xlocal,zlocal);
			fAssembly->Assemble(isub,zlocal,z);
			//         z.Print("Resultado intermediario");
		}		
	}
	else {
		TPZAutoPointer<TPZDohrAssembleList> assemblelist = new TPZDohrAssembleList(fGlobal.size(),z,this->fAssembly);
		
		TPZDohrThreadMultList<TSubStruct> multwork(x,alpha,fAssembly,assemblelist);
		typename std::list<TPZAutoPointer<TSubStruct> >::const_iterator iter;
		int isub=0;
		for (iter=fGlobal.begin(); iter!=fGlobal.end(); iter++,isub++) {
			TPZDohrThreadMultData<TSubStruct> data(isub,*iter);
			multwork.AddItem(data);
		}
		TPZVec<pthread_t> AllThreads(fNumThreads+1);
		int i;
		for (i=0; i<fNumThreads; i++) {
			pthread_create(&AllThreads[i+1], 0, TPZDohrThreadMultList<TSubStruct>::ThreadWork, &multwork);
		}
		pthread_create(&AllThreads[0], 0, TPZDohrAssembleList::Assemble, assemblelist.operator->());
		
		for (i=0; i<fNumThreads+1; i++) {
			void *result;
			pthread_join(AllThreads[i], &result);
		}
	}
	tempo.fMultiply.Push(mult.ReturnTimeDouble());
}

template<class TSubStruct>
void TPZDohrMatrix<TSubStruct>::Initialize() 
{
	std::cout << "Number of substructures " << fGlobal.size() << std::endl;
	tempo.fNumSub = fGlobal.size();																// alimenta timeTemp com o numero de substruturas
	TPZFMatrix diag(Rows(),1,0.);
	typename SubsList::iterator iter;
	int isub = 0;
	for (iter=fGlobal.begin();iter!=fGlobal.end();iter++,isub++) {
        //Basic initialization for each substructure (compute the matrices)
        //(*iter)->Initialize();
		TPZFMatrix diaglocal;
        (*iter)->ContributeDiagonalLocal(diaglocal);
		LOGPZ_DEBUG(logger,"Before assemble diagonal")
		this->fAssembly->Assemble(isub,diaglocal,diag);
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Substructure " << isub << " ";
			diag.Print("Global Diagonal matrix",sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
        std::cout << '*';
        std::cout.flush();
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		diag.Print("Global Diagonal matrix",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	std::cout << std::endl;
	for (iter=fGlobal.begin(),isub=0;iter!=fGlobal.end();iter++,isub++) {
        //Computes the Weights for each substructure
		TPZFMatrix diaglocal;
		this->fAssembly->Extract(isub,diag,diaglocal);
        (*iter)->ComputeWeightsLocal(diaglocal);
		
	}
}

/**
 * Adjust the residual to zero the residual of the internal connects
 */
template<class TSubStruct>
void TPZDohrMatrix<TSubStruct>::AdjustResidual(TPZFMatrix &res)
{
	typename SubsList::iterator iter;
	for (iter=fGlobal.begin();iter!=fGlobal.end();iter++) {
		(*iter)->AdjustResidual(res);
	}
}

/**
 * Add the solution corresponding to the internal residual
 */
template<class TSubStruct>
void TPZDohrMatrix<TSubStruct>::AddInternalSolution(TPZFMatrix &solution)
{
	typename SubsList::iterator iter;
	for (iter=fGlobal.begin();iter!=fGlobal.end();iter++) {
		(*iter)->AddInternalSolution(solution);
	}
}

template<class TSubStruct>
void *TPZDohrThreadMultList<TSubStruct>::ThreadWork(void *ptr)
{
	TPZDohrThreadMultList<TSubStruct> *myptr = (TPZDohrThreadMultList<TSubStruct> *) ptr;
	TPZDohrThreadMultData<TSubStruct> runner = myptr->PopItem();
	while (runner.IsValid()) {
		TPZFMatrix xlocal;
		myptr->fAssembly->Extract(runner.fisub,*(myptr->fInput),xlocal);
		TPZAutoPointer<TPZDohrAssembleItem> assembleItem = new TPZDohrAssembleItem(runner.fisub,xlocal.Rows(),xlocal.Cols());
		runner.fSub->ContributeKULocal(myptr->fAlpha,xlocal,assembleItem->fAssembleData);
		myptr->fAssemblyStructure->AddItem(assembleItem);
		runner = myptr->PopItem();
	}
	return ptr;
}

template class TPZDohrMatrix<TPZDohrSubstruct>;
template class TPZDohrMatrix<TPZDohrSubstructCondense>;

