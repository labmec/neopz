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
void TPBRThermalDiscretization::ComputeStiffness()
{
	TPZManVector<int> skyline(fNElements+1,0);
	for (int i=1; i<=fNElements; i++) {
		skyline[i] = i-1;
	}
	TPZSkylMatrix	*skyl = new TPZSkylMatrix(fNElements+1,skyline);
	REAL delx = fDomainSize/fNElements;
	REAL diag = fTimeStep*fK/delx+fCp*fDensity*delx/2; // [KJ/(m2 C)]
	REAL offdiag = -fTimeStep*fK/delx; // [KJ/(m2 C)]
	skyl->PutVal(0, 0, diag);
	skyl->PutVal(0, 1, offdiag);
	skyl->PutVal(fNElements, fNElements, diag);
	diag *= 2.;
	for (int i=1; i<fNElements; i++) {
		skyl->PutVal(i, i, diag);
		skyl->PutVal(i, i+1, offdiag);
	}
//	skyl->Print("Band matrix",std::cout);
	TPZStepSolver *step = new TPZStepSolver(skyl);
	step->SetDirect(ECholesky);
	fSolver = step;
	fUnitFluxSolution.Redim(fNElements+1,1);
	TPZFMatrix rhs(fNElements+1,1,0.);
	rhs(0,0) = 1.; // [KJ/m2]
	fSolver->Solve(rhs, fUnitFluxSolution, 0);
//	skyl->Print("Band matrix after decompose",std::cout);
//	fUnitFluxSolution.Print("UnitFlux",std::cout);
//    std::cout << "Energy of unit flux " << Energy(fUnitFluxSolution) << std::endl;
}

/// Compute the next solution
void TPBRThermalDiscretization::NextSolution(REAL inletTemp, TPZFMatrix &prevSol, TPZFMatrix &nextSol, REAL &flux)
{
	TPZFNMatrix<201> rhs(fNElements+1,1);
	REAL delx = fDomainSize/fNElements;
	rhs(0,0) = prevSol(0,0)*fCp*fDensity*delx/2.;
	rhs(fNElements,0) = prevSol(fNElements,0)*fCp*fDensity*delx/2.; // [KJ/m2]
	for (int i=1; i<fNElements; i++) {
		rhs(i,0) = prevSol(i,0)*delx*fCp*fDensity; // [KJ/m2]
	}
	nextSol.SetSize(fNElements+1, 1); //[C]
	fSolver->Solve(rhs,nextSol);
	REAL scale = (inletTemp-nextSol(0,0))/fUnitFluxSolution(0,0); // [KJ/m2]
	flux = scale; // [KJ/m2]
	nextSol += scale*fUnitFluxSolution;
}

