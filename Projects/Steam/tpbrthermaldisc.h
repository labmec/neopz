/*
 *  tpbrthermaldisc.h
 *  IP3D_v5
 *
 *  Created by Philippe Devloo on 1/9/11.
 *  Copyright 2011 UNICAMP. All rights reserved.
 *
 */
#ifndef TPBRTHERMALDISC
#define TPBRTHERMALDISC

#include "pzreal.h"
#include "pzsolve.h"

/// Definition of the properties and mesh size of the confining layer
class TPBRThermalDisc 
{
private:
	
	/// total size of the domain
	REAL fDomainSize;
	
	/// number of elements
	int fNElements;
	
	/// thermal conductivity (unit : )
	REAL fK;
	
	/// thermal capacity
	REAL fCp;
	
	/// initial temperature (C)
	REAL fInitialTemp;
	
	/// matrix solution procedure for computing the next solution
	TPZAutoPointer<TPZSolver> fSolver;
	
	/// timestep
	REAL fTimeStep;
	
	/// variation of the temperature with the thermal flux
	TPZFMatrix fUnitFluxSolution;
	
public:
	
	/// constructor
	TPBRThermalDisc(REAL domainsize, int nelements, REAL cp, REAL K, REAL initialtemp) : 
		fDomainSize(domainsize), fNElements(nelements), fK(K), fCp(cp), fInitialTemp(initialtemp), fTimeStep(-1), fUnitFluxSolution()
	{
	}
	
	/// Set the timestep
	void SetTimeStep(REAL delt)
	{
		fTimeStep = delt;
	}
	
	/// Compute the stiffness matrix
	void ComputeStiffness();
	
	/// Compute the derivative of the heat flux rate with respect to the inlet temperature
	REAL DQDT()
	{
		return 1./fUnitFluxSolution(0,0);
	}
	
	/// Compute the energy associated with the solution
	REAL Energy(TPZFMatrix &solution)
	{
#ifdef DEBUG
		if (solution.Rows() != fNElements+1) {
			DebugStop();
		}
#endif
		REAL delx = fDomainSize/fNElements;
		REAL VarEnergy = (solution(0,0)+solution(fNElements,0))/2;
		int i;
		for (i=1; i<fNElements; i++) {
			VarEnergy += solution(i,0);
		}
		VarEnergy *= delx*fCp;
		return VarEnergy;
	}
	
	/// Compute the next solution
	void NextSolution(REAL inletTemp, TPZFMatrix &prevSol, TPZFMatrix &nextSol, REAL &flux);
};
#endif