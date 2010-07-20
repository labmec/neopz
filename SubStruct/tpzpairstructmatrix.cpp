/*
 *  tpzpairstructmatrix.cpp
 *  SubStruct
 *
 *  Created by Philippe Devloo on 20/04/09.
 *  Copyright 2009 UNICAMP. All rights reserved.
 *
 */

#include "tpzpairstructmatrix.h"
#include "TPZTimer.h"
#include "pzelmat.h"
#include "pzcompel.h"
#include "pzstrmatrix.h"

using namespace std;

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.tpzpairstructmatrix"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
#endif

void TPZPairStructMatrix::Assemble(int mineq, int maxeq, TPZMatrix *first, TPZMatrix *second, TPZFMatrix &rhs)
{
	int iel;
	TPZCompMesh &mesh = *fMesh;
	int nelem = mesh.NElements();
	TPZElementMatrix ek(&mesh, TPZElementMatrix::EK),ef(&mesh, TPZElementMatrix::EF);
	
	TPZTimer calcstiff("Computing the stiffness matrices");
	TPZTimer assemble("Assembling the stiffness matrices");
	TPZAdmChunkVector<TPZCompEl *> &elementvec = mesh.ElementVec();
	
	for(iel=0; iel < nelem; iel++) {
		TPZCompEl *el = elementvec[iel];
		if(!el) continue;
		calcstiff.start();
		
		el->CalcStiff(ek,ef);
		
		calcstiff.stop();
		assemble.start();
		
		if(!el->HasDependency()) {
			ek.ComputeDestinationIndices();
			if(mineq != -1 && maxeq != -1)
			{
				TPZStructMatrix::FilterEquations(ek.fSourceIndex,ek.fDestinationIndex,mineq,maxeq);
			}
			first->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
			rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
			PermuteScatter(ek.fDestinationIndex);
			second->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
#ifdef LOG4CXX
			if(loggerel->isDebugEnabled())
			{
				std::stringstream sout;
				ek.fMat.Print("Element stiffness matrix",sout);
				ef.fMat.Print("Element right hand side", sout);
				LOGPZ_DEBUG(loggerel,sout.str())
			}
#endif
		} else {
			// the element has dependent nodes
			ek.ApplyConstraints();
			ef.ApplyConstraints();
			ek.ComputeDestinationIndices();
			if(mineq != -1 & maxeq != -1)
			{
				TPZStructMatrix::FilterEquations(ek.fSourceIndex,ek.fDestinationIndex,mineq,maxeq);
			}
			first->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
			rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
			PermuteScatter(ek.fDestinationIndex);
			second->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
		}
		
#ifndef WIN32
		assemble.stop();
#endif		
	}//fim for iel
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Number of equations " << first->Rows() << std::endl;
		sout << calcstiff.processName() << " " << calcstiff << std::endl;
		sout << assemble.processName() << " " << assemble;
		/*    stiffness.Print("Matriz de Rigidez: ",sout);
		 rhs.Print("Vetor de Carga: ",sout);*/
		LOGPZ_DEBUG(logger,sout.str().c_str());
	}
#endif
	
}

void TPZPairStructMatrix::PermuteScatter(TPZVec<int> &index)
{
	int nel = index.NElements();
	int iel;
	for(iel = 0; iel<nel; iel++)
	{
		index[iel] = fPermuteScatter[index[iel]];
	}
}
