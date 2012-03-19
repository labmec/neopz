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
class TPBRThermalDiscretization 
{
private:
	
	/// total size of the domain [m]
	REAL fDomainSize;
	
	/// number of elements
	int fNElements;
	
	/// thermal conductivity [KJ/(m C s]
	REAL fK;
	
	/// thermal capacity [KJ/(kg C]
	REAL fCp;
    
    /// density of the material
    REAL fDensity;
	
	/// initial temperature (C)
	REAL fInitialTemp;
	
	/// matrix solution procedure for computing the next solution
	TPZAutoPointer<TPZSolver<REAL> > fSolver;
	
	/// timestep [s]
	REAL fTimeStep;
	
	/// variation of the temperature with the thermal flux
	TPZFMatrix<REAL> fUnitFluxSolution;
	
public:
	/// create an invalid object
    TPBRThermalDiscretization() :fDomainSize(-1.),fNElements(0),fK(0.),fCp(0.),fDensity(-1.),fInitialTemp(0.)
    {
    }
	/// constructor
    /**
     * domainsize [m]
     * cp thermal capacity [KJ/Kg]
     * K conductivity [KJ/(s m C)]
     * density [Kg/m3]
     * initialtemp [C]
     */
	TPBRThermalDiscretization(REAL domainsize, int nelements, REAL cp, REAL K, REAL density, REAL initialtemp) : 
		fDomainSize(domainsize), fNElements(nelements), fK(K), fCp(cp), fDensity(density), fInitialTemp(initialtemp), fTimeStep(-1), fUnitFluxSolution()
	{
	}
    
    TPBRThermalDiscretization &operator=(const TPBRThermalDiscretization &copy)
    {
        fDomainSize = copy.fDomainSize;
        fNElements = copy.fNElements;
        fK = copy.fK;
        fCp = copy.fCp;
        fDensity = copy.fDensity;
        fInitialTemp = copy.fInitialTemp;
        fSolver = copy.fSolver;
        fTimeStep = copy.fTimeStep;
        fUnitFluxSolution = copy.fUnitFluxSolution;
        return *this;
    }
	
    TPBRThermalDiscretization(const TPBRThermalDiscretization &copy)
    {
        fDomainSize = copy.fDomainSize;
        fNElements = copy.fNElements;
        fK = copy.fK;
        fCp = copy.fCp;
        fDensity = copy.fDensity;
        fInitialTemp = copy.fInitialTemp;
        fSolver = copy.fSolver;
        fTimeStep = copy.fTimeStep;
        fUnitFluxSolution = copy.fUnitFluxSolution;
    }
	
	/// Set the timestep
	void SetTimeStep(REAL delt)
	{
		fTimeStep = delt;
        ComputeStiffness();
	}
    
    void InitializeSolution(TPZFMatrix<REAL> &sol)
    {
        sol.Redim(fNElements+1, 1);
        sol += fInitialTemp;
    }
	
	/// Compute the stiffness matrix
	void ComputeStiffness();
	
	/// Compute the derivative of the heat flux rate with respect to the inlet temperature
	REAL DQDT()
	{
		return 1./fUnitFluxSolution(0,0); // [KJ/(m2 C)]
	}
	
	/// Compute the energy associated with the solution
	REAL Energy(TPZFMatrix<REAL> &solution)
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
		VarEnergy *= delx*fCp*fDensity; // [KJ/m2]
		return VarEnergy;
	}
	
	/// Compute the next solution
	void NextSolution(REAL inletTemp, TPZFMatrix<REAL> &prevSol, TPZFMatrix<REAL> &nextSol, REAL &flux);
};
#endif