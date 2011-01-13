/*
 *  tpbrthermaldisc.cpp
 *  IP3D_v5
 *
 *  Created by Philippe Devloo on 1/9/11.
 *  Copyright 2011 UNICAMP. All rights reserved.
 *
 */


#include "tpbrthermaldisc.h"
#include "pzsbndmat.h"
#include "pzstepsolver.h"

/// Compute the stiffness matrix
void TPBRThermalDisc::ComputeStiffness()
{
	TPZSBMatrix	*bnd = new TPZSBMatrix(fNElements+1,1);
	REAL delx = fDomainSize/fNElements;
	REAL diag = fTimeStep*fK/delx+fCp*delx/2;
	REAL offdiag = -fTimeStep*fK/delx;
	bnd->PutVal(0, 0, diag);
	bnd->PutVal(0, 1, offdiag);
	bnd->PutVal(fNElements, fNElements, diag);
	diag *= 2.;
	for (int i=1; i<fNElements; i++) {
		bnd->PutVal(i, i, diag);
		bnd->PutVal(i, i+1, offdiag);
	}
	TPZStepSolver *step = new TPZStepSolver(bnd);
	step->SetDirect(ECholesky);
	fSolver = step;
	fUnitFluxSolution.Redim(fNElements+1,1);
	TPZFMatrix rhs(fNElements,1,0.);
	rhs(0,0) = fTimeStep;
	fSolver->Solve(rhs, fUnitFluxSolution, 0);
}

/// Compute the next solution
void TPBRThermalDisc::NextSolution(REAL inletTemp, TPZFMatrix &prevSol, TPZFMatrix &nextSol)
{
	TPZFMatrix rhs(fNElements+1,1);
	REAL delx = fDomainSize/fNElements;
	rhs(0,0) = prevSol(0,0)*fCp*delx/2.;
	rhs(fNElements,0) = prevSol(fNElements,0)*delx/2.;
	for (int i=1; i<fNElements; i++) {
		rhs(i,0) = prevSol(i,0)*delx;
	}
	nextSol.SetSize(fNElements+1, 1);
	fSolver->Solve(prevSol,nextSol);
	REAL scale = (inletTemp-nextSol(0,0))/fUnitFluxSolution(0,0);
	nextSol += scale*fUnitFluxSolution;
}

