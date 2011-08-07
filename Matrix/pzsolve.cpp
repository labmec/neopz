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
TPZSolver::~TPZSolver()
{
}

TPZMatrixSolver::TPZMatrixSolver(TPZAutoPointer<TPZMatrix> Refmat) :
fScratch()
{
	fContainer = Refmat;
}

TPZMatrixSolver::TPZMatrixSolver() :
fScratch()
{
}

//misael
TPZMatrixSolver::TPZMatrixSolver(const TPZMatrixSolver &Source) :
fScratch()
{
	fReferenceMatrix = Source.fReferenceMatrix;
	fContainer = Source.fContainer;
}

// philippe 6/2/97
TPZMatrixSolver::~TPZMatrixSolver()
{
}

void TPZMatrixSolver::SetMatrix(TPZAutoPointer<TPZMatrix> Refmat)
{
	fContainer = Refmat;
}

void TPZMatrixSolver::ResetMatrix()
{
	TPZAutoPointer<TPZMatrix> reset;
	fContainer = reset;
}

void TPZMatrixSolver::ShareMatrix(TPZMatrixSolver &other)
{
	if (this == &other)
		return;
	fContainer = other.fContainer;
}

void TPZMatrixSolver::Write(TPZStream &buf, int withclassid)
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

void TPZMatrixSolver::Read(TPZStream &buf, void *context)
{
	TPZSaveable::Read(buf,context);
	fContainer = dynamic_cast<TPZMatrix *>(TPZSaveable::Restore(buf, context));
	fReferenceMatrix = dynamic_cast<TPZMatrix *>(TPZSaveable::Restore(buf, context));
}

//template class TPZRestoreClass< TPZMatrixSolver, TPZMATRIXSOLVER_ID>;

