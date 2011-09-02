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
void TPZDohrAssembly::Assemble(int isub, const TPZFMatrix &local, TPZFMatrix &global)
{
	TPZVec<int> &avec = fFineEqs[isub];
	int neq = avec.NElements();
	int ieq;
	for(ieq=0; ieq<neq; ieq++)
	{
		global(avec[ieq],0) += local.GetVal(ieq,0);
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
void TPZDohrAssembly::Extract(int isub, const TPZFMatrix &global, TPZFMatrix &local)
{
	TPZVec<int> &avec = fFineEqs[isub];
	int neq = avec.NElements();
	int ieq;
	local.Resize(neq,1);
	for(ieq=0; ieq<neq; ieq++)
	{
		local(ieq,0) = global.GetVal(avec[ieq],0);
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
void TPZDohrAssembly::AssembleCoarse(int isub, const TPZFMatrix &local, TPZFMatrix &global)
{
	TPZVec<int> &avec = fCoarseEqs[isub];
	int neq = avec.NElements();
	int ieq;
	for(ieq=0; ieq<neq; ieq++)
	{
		global(avec[ieq],0) += local.GetVal(ieq,0);
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
void TPZDohrAssembly::ExtractCoarse(int isub, const TPZFMatrix &global, TPZFMatrix &local)
{
	TPZVec<int> &avec = fCoarseEqs[isub];
	int neq = avec.NElements();
	int ieq;
	local.Resize(neq,1);
	for(ieq=0; ieq<neq; ieq++)
	{
		local(ieq,0) = global.GetVal(avec[ieq],0);
	}
}
