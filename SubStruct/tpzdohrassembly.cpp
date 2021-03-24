/**
 * @file
 * @brief Contains the implementation of the TPZDohrAssembly methods. 
 * @author Philippe Devloo
 * @since 04/03/2009
 */

#include "tpzdohrassembly.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("substruct.dohrassembly");
#endif

// sum the values in the local matrix into the global matrix
template<class TVar>
void TPZDohrAssembly<TVar>::Assemble(int isub, const TPZFMatrix<TVar> &local, TPZFMatrix<TVar> &global) const
{
	TPZVec<int> &avec = fFineEqs[isub];
	int neq = avec.NElements();
    int ncols = local.Cols();
	int ieq;
    for (int ic=0; ic<ncols; ic++) 
    {
        for(ieq=0; ieq<neq; ieq++)
        {
            global(avec[ieq],ic) += local.GetVal(ieq,ic);
        }
    }
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << "Assembling destination indices " << avec << std::endl;
		local.Print("Input vector",sout);
		global.Print("Resulting vector",sout);
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str());
		}
	}
#endif
}

// extract the values from the global matrix into the local matrix
template<class TVar>
void TPZDohrAssembly<TVar>::Extract(int isub, const TPZFMatrix<TVar> &global, TPZFMatrix<TVar> &local) const
{
	TPZVec<int> &avec = fFineEqs[isub];
	int neq = avec.NElements();
    int ncols = global.Cols();
	int ieq;
	local.Resize(neq,ncols);
    for (int ic=0; ic<ncols; ic++) 
    {
        for(ieq=0; ieq<neq; ieq++)
        {
            local(ieq,ic) = global.GetVal(avec[ieq],ic);
        }
    }
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << "sub structure " << isub << " Extracting destination indices " << avec << std::endl;
		local.Print("extracted vector",sout);
		global.Print("Global vector",sout);
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str());
		}
	}
#endif
}

// sum the values in the local matrix into the global matrix
template<class TVar>
void TPZDohrAssembly<TVar>::AssembleCoarse(int isub, const TPZFMatrix<TVar> &local, TPZFMatrix<TVar> &global) const
{
	TPZVec<int> &avec = fCoarseEqs[isub];
	int neq = avec.NElements();
	int ieq;
    int ncols = local.Cols();
    for (int ic=0; ic<ncols; ic++) 
    {
        for(ieq=0; ieq<neq; ieq++)
        {
            global(avec[ieq],ic) += local.GetVal(ieq,ic);
        }
    }
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << "Assembling destination indices " << avec << std::endl;
		local.Print("Input vector",sout);
		global.Print("Resulting vector",sout);
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str());
		}
	}
#endif
}

// extract the values from the global matrix into the local matrix
template<class TVar>
void TPZDohrAssembly<TVar>::ExtractCoarse(int isub, const TPZFMatrix<TVar> &global, TPZFMatrix<TVar> &local) const
{
	TPZVec<int> &avec = fCoarseEqs[isub];
	int neq = avec.NElements();
	int ieq;
    int ncols = global.Cols();
	local.Resize(neq,ncols);
    for (int ic=0; ic<ncols; ic++) 
    {
        for(ieq=0; ieq<neq; ieq++)
        {
            local(ieq,ic) = global.GetVal(avec[ieq],ic);
        }
    }
}

template<class TVar>
int TPZDohrAssembly<TVar>::ClassId() const {
    return Hash("TPZDohrAssembly") ^ ClassIdOrHash<TVar>() << 1;
}

/** @brief method for streaming the object to a stream */
template<class TVar>
void TPZDohrAssembly<TVar>::Write(TPZStream &buf, int withclassid) const
{
    int nfine = fFineEqs.size();
    buf.Write(&nfine,1);
    for (int f=0; f<nfine; f++) {
        buf.Write( fFineEqs[f]);
    }
    int ncoarse = fCoarseEqs.size();
    buf.Write(&ncoarse,1);
    for (int nc=0; nc<ncoarse; nc++) {
        buf.Write( fCoarseEqs[nc]);
    }
}

/** @brief method for reading the object for a stream */
template<class TVar>
void TPZDohrAssembly<TVar>::Read(TPZStream &buf, void *context)
{
    int nfine;
    buf.Read(&nfine);
    fFineEqs.resize(nfine);
    for (int f=0; f<nfine; f++) {
        buf.Read( fFineEqs[f]);
    }
    int ncoarse;
    buf.Read(&ncoarse);
    fCoarseEqs.resize(ncoarse);
    for (int nc=0; nc<ncoarse; nc++) {
        buf.Read( fCoarseEqs[nc]);
    }    
}

template class TPZDohrAssembly<float>;
template class TPZDohrAssembly<double>;
template class TPZDohrAssembly<long double>;

template class TPZDohrAssembly<std::complex<float> >;
template class TPZDohrAssembly<std::complex<double> >;
template class TPZDohrAssembly<std::complex<long double> >;
