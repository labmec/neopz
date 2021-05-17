/**
 * @file
 * @brief Contains the implementation of the TPZMatrixSolver methods.
 */

#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.tpzmatred");
#endif

#include "TPZMatrixSolver.h"

#include "Hash/TPZHash.h"
#include "TPZPersistenceManager.h"

#include <stdlib.h>
using namespace std;

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

template <class TVar>
TPZMatrixSolver<TVar>::TPZMatrixSolver(const TPZMatrixSolver<TVar> &Source) :
fScratch()
{
	fReferenceMatrix = Source.fReferenceMatrix;
	fContainer = Source.fContainer;
}

template <class TVar>
TPZMatrixSolver<TVar>::~TPZMatrixSolver()
{
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
void TPZMatrixSolver<TVar>::Write(TPZStream &buf, int withclassid) const {
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "Entering " << __PRETTY_FUNCTION__;
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    TPZSolver::Write(buf, withclassid);
    if (fContainer) {
#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "fContainer AutoPointer valid on " << __PRETTY_FUNCTION__;
            LOGPZ_DEBUG(logger, sout.str());
        }
#endif
    }
    TPZPersistenceManager::WritePointer(fContainer.operator->(), &buf);
    
    if (fReferenceMatrix) {
#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "fReferenceMatrix AutoPointer valid! It Shouldn't ! Expect Trouble " << __PRETTY_FUNCTION__;
            LOGPZ_WARN(logger, sout.str());
        }
#endif
    }
    TPZPersistenceManager::WritePointer(fReferenceMatrix.operator->(), &buf);
    
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "Leaving" << __PRETTY_FUNCTION__;
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
}

template <class TVar>
void TPZMatrixSolver<TVar>::Read(TPZStream &buf, void *context)
{
	TPZSolver::Read(buf,context);
	fContainer = TPZAutoPointerDynamicCast<TPZMatrix<TVar>>(TPZPersistenceManager::GetAutoPointer(&buf));
	fReferenceMatrix = TPZAutoPointerDynamicCast<TPZMatrix<TVar>>(TPZPersistenceManager::GetAutoPointer(&buf));
}

template<class TVar>
int TPZMatrixSolver<TVar>::ClassId() const{
    return Hash("TPZMatrixSolver") ^ ClassIdOrHash<TVar>() ^ TPZSolver::ClassId() << 2;
}


template class TPZMatrixSolver<float>;
template class TPZMatrixSolver<std::complex<float> >;

template class TPZMatrixSolver<double>;
template class TPZMatrixSolver<std::complex<double> >;

template class TPZMatrixSolver<long double>;
template class TPZMatrixSolver<std::complex<long double> >;

template class TPZMatrixSolver<Fad<float> >;
template class TPZMatrixSolver<Fad<double> >;
template class TPZMatrixSolver<Fad<long double> >;