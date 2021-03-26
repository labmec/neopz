/*
 *  tpbrsolutionlist.h
 *  IP3D_v5
 *
 *  Created by Philippe Devloo on 1/9/11.
 *  Copyright 2011 UNICAMP. All rights reserved.
 *
 */
#ifndef TPBRSOLUTIONLIST
#define TPBRSOLUTIONLIST

#include <vector>
#include "pzreal.h"
#include "tpbrthermalsolution.h"
#include "tpbrthermaldisc.h"

/// Contains a list of thermal solutions
class TPBRSolutionList 
{
private:
	/// list of solutions
	std::vector<TPBRThermalSolution> fList;
	/// discretization
	TPBRThermalDiscretization fDiscretization;
    /// current time step
    REAL fDelt;

public:
    
    /// initialize the discretization
    void SetDiscretization(TPBRThermalDiscretization &discretization)
    {
        fDiscretization = discretization;
    }
    
	/// total energy of the solution list
	REAL Energy();

	/// Compute the solution for the next timestep
	void AdvanceAllSolutions(REAL delt, REAL inletTemp, REAL &flux, REAL &DQDT, bool storesolution);
	
    /// Set the timestep
    void SetDelt(REAL delt);

    /// Compute the variation of the flux with respect to the inlet temperature
//	REAL DQDT(REAL delt, REAL inletTemp);
    
	/// total energy of the solution list
	REAL Energy(int icell);
    
	/// Compute the solution for the next timestep
	void AdvanceSolution(int icell, REAL inletTemp, REAL &flux, REAL &DQDT, bool storesolution);
	
    /// Compute the variation of the flux with respect to the inlet temperature
    /**
     * inletTemp : temperature at the inlet (input)
     * Flux : flux corresponding to the inlet temperature
     * return : variation of the flux with temperature
     */
	REAL DQDT(int icell, REAL inletTemp, REAL &Flux);
    
    /// Add a solution to the list
    void AddSolution(TPBRThermalSolution &nextsol);
    
    /// Remove all solutions
    void ClearSolutions()
    {
        fList.resize(0);
    }
    
    /// Set the timestep and compute the stiffness matrix
    void SetTimeStep(REAL delt);
    
};
#endif
