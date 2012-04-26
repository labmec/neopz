/**
 * @file
 * @brief Contains the implementation of the TPZDohrAssembly methods. 
 */
/*
 *  tpzdohrassembly.cpp
 *  SubStruct
 *
 *  Created by Philippe Devloo on 04/03/09.
 *  Copyright 2009 UNICAMP. All rights reserved.
 *
 */

#include "tpzdohrassembly.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("substruct.dohrassembly"));
#endif

// sum the values in the local matrix into the global matrix
void TPZDohrAssembly::Assemble(int isub, const TPZFMatrix<REAL> &local, TPZFMatrix<REAL> &global)
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
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Assembling destination indices " << avec << std::endl;
		local.Print("Input vector",sout);
		global.Print("Resulting vector",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

// extract the values from the global matrix into the local matrix
void TPZDohrAssembly::Extract(int isub, const TPZFMatrix<REAL> &global, TPZFMatrix<REAL> &local)
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
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "sub structure " << isub << " Extracting destination indices " << avec << std::endl;
		local.Print("extracted vector",sout);
		global.Print("Global vector",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

// sum the values in the local matrix into the global matrix
void TPZDohrAssembly::AssembleCoarse(int isub, const TPZFMatrix<REAL> &local, TPZFMatrix<REAL> &global)
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
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Assembling destination indices " << avec << std::endl;
		local.Print("Input vector",sout);
		global.Print("Resulting vector",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

// extract the values from the global matrix into the local matrix
void TPZDohrAssembly::ExtractCoarse(int isub, const TPZFMatrix<REAL> &global, TPZFMatrix<REAL> &local)
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

/**
 * @brief For each substructure the equation numbering of the substructures
 * 
 * The order of the equations follows the ordering of the connects
 */
//TPZVec< TPZVec< int > > fFineEqs;

/** @brief For each substructure the equation numbering of the coarse equations */
//TPZVec< TPZVec< int > > fCoarseEqs;


/** @brief method for streaming the object to a stream */
void TPZDohrAssembly::Write(TPZStream &out)
{
    int nfine = fFineEqs.size();
    out.Write(&nfine,1);
    for (int f=0; f<nfine; f++) {
        TPZSaveable::WriteObjects(out, fFineEqs[f]);
    }
    int ncoarse = fCoarseEqs.size();
    out.Write(&ncoarse,1);
    for (int nc=0; nc<ncoarse; nc++) {
        TPZSaveable::WriteObjects(out, fCoarseEqs[nc]);
    }
}

/** @brief method for reading the object for a stream */
void TPZDohrAssembly::Read(TPZStream &input)
{
    int nfine;
    input.Read(&nfine);
    fFineEqs.resize(nfine);
    for (int f=0; f<nfine; f++) {
        TPZSaveable::ReadObjects(input, fFineEqs[f]);
    }
    int ncoarse;
    input.Read(&ncoarse);
    fCoarseEqs.resize(ncoarse);
    for (int nc=0; nc<ncoarse; nc++) {
        TPZSaveable::ReadObjects(input, fCoarseEqs[nc]);
    }    
}

