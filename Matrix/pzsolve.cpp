/**
 * @file
 * @brief Contains the implementation of the TPZMatrixSolver methods.
 */

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzmatred"));
#endif

#include "pzsolve.h"

#include <stdlib.h>
using namespace std;

/**
 
 
 Destructor
 */
template <class TVar>
TPZSolver<TVar>::~TPZSolver()
{
}

template <class TVar>
TPZMatrixSolver<TVar>::TPZMatrixSolver(TPZAutoPointer<TPZMatrix<TVar> > Refmat) :
fScratch()
{
	fContainer = Refmat;
}

template<class TVar>
TPZMatrixSolver<TVar>::TPZMatrixSolver() :
fScratch()
{
}

//misael
template <class TVar>
TPZMatrixSolver<TVar>::TPZMatrixSolver(const TPZMatrixSolver<TVar> &Source) :
fScratch()
{
	fReferenceMatrix = Source.fReferenceMatrix;
	fContainer = Source.fContainer;
}

// philippe 6/2/97
template <class TVar>
TPZMatrixSolver<TVar>::~TPZMatrixSolver()
{
}

template <class TVar>
void TPZMatrixSolver<TVar>::SetMatrix(TPZAutoPointer<TPZMatrix<TVar> > Refmat)
{
	fContainer = Refmat;
}

template <class TVar>
void TPZMatrixSolver<TVar>::ResetMatrix()
{
	TPZAutoPointer<TPZMatrix<TVar> > reset;
	fContainer = reset;
}

template <class TVar>
void TPZMatrixSolver<TVar>::ShareMatrix(TPZMatrixSolver<TVar> &other)
{
	if (this == &other)
		return;
	fContainer = other.fContainer;
}
template <class TVar>
void TPZMatrixSolver<TVar>::Write(TPZStream &buf, int withclassid)
{
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Entering " << __PRETTY_FUNCTION__;
		LOGPZ_DEBUG(logger,sout.str());
	}
#endif
	TPZSaveable::Write(buf,withclassid);
	if(fContainer)
	{
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "fContainer AutoPointer valid on " << __PRETTY_FUNCTION__;
			LOGPZ_DEBUG(logger,sout.str());
		}
#endif
		
		fContainer->Write(buf, 1);
	}
	else
	{
		int flag = -1;
		buf.Write(&flag, 1);
	}
	if(fReferenceMatrix)
	{
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "fReferenceMatrix AutoPointer valid! It Shouldn't ! Expect Trouble " << __PRETTY_FUNCTION__;
			LOGPZ_WARN(logger,sout.str());
		}
#endif
		fReferenceMatrix->Write(buf, 1);
	}
	else
	{
		int flag = -1;
		buf.Write(&flag, 1);
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Leaving" << __PRETTY_FUNCTION__;
		LOGPZ_DEBUG(logger,sout.str());
	}
#endif
}

template <class TVar>
void TPZMatrixSolver<TVar>::Read(TPZStream &buf, void *context)
{
	TPZSaveable::Read(buf,context);
	fContainer = dynamic_cast<TPZMatrix<TVar> *>(TPZSaveable::Restore(buf, context));
	fReferenceMatrix = dynamic_cast<TPZMatrix<TVar> *>(TPZSaveable::Restore(buf, context));
}

template class TPZMatrixSolver<REAL>;
//template class TPZRestoreClass< TPZMatrixSolver, TPZMATRIXSOLVER_ID>;


