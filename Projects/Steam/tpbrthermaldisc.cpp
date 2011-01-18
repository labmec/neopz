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
#include "pzskylmat.h"
#include "pzstepsolver.h"

/// Compute the stiffness matrix
void TPBRThermalDisc::ComputeStiffness()
{
	TPZManVector<int> skyline(fNElements+1,0);
	for (int i=1; i<=fNElements; i++) {
		skyline[i] = i-1;
	}
	TPZSkylMatrix	*skyl = new TPZSkylMatrix(fNElements+1,skyline);
	REAL delx = fDomainSize/fNElements;
	REAL diag = fTimeStep*fK/delx+fCp*delx/2;
	REAL offdiag = -fTimeStep*fK/delx;
	skyl->PutVal(0, 0, diag);
	skyl->PutVal(0, 1, offdiag);
	skyl->PutVal(fNElements, fNElements, diag);
	diag *= 2.;
	for (int i=1; i<fNElements; i++) {
		skyl->PutVal(i, i, diag);
		skyl->PutVal(i, i+1, offdiag);
	}
	skyl->Print("Band matrix",std::cout);
	TPZStepSolver *step = new TPZStepSolver(skyl);
	step->SetDirect(ECholesky);
	fSolver = step;
	fUnitFluxSolution.Redim(fNElements+1,1);
	TPZFMatrix rhs(fNElements+1,1,0.);
	rhs(0,0) = fTimeStep;
	fSolver->Solve(rhs, fUnitFluxSolution, 0);
	skyl->Print("Band matrix after decompose",std::cout);
	fUnitFluxSolution.Print("UnitFlux",std::cout);
}

/// Compute the next solution
void TPBRThermalDisc::NextSolution(REAL inletTemp, TPZFMatrix &prevSol, TPZFMatrix &nextSol, REAL &flux)
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
	flux = scale;
	nextSol += scale*fUnitFluxSolution;
}

