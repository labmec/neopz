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

#include <list>
#include "pzreal.h"
#include "tpbrthermalsolution.h"
#include "tpbrthermaldisc.h"

/// Contains a list of thermal solutions
class TPBRSolutionList 
{
private:
	/// list of solutions
	std::list<TPBRThermalSolution> fList;
	/// discretization
	TPBRThermalDiscretization fDiscretization;

public:
    
    /// initialize the discretization
    void SetDiscretization(TPBRThermalDiscretization &discretization)
    {
        fDiscretization = discretization;
    }
    
	/// total energy of the solution list
	REAL Energy();
	/// Compute the solution for the next timestep
	void AdvanceSolution(REAL delt, REAL inletTemp, REAL &flux, REAL &DQDT, bool storesolution);
	
    /// Compute the variation of the flux with respect to the inlet temperature
	REAL DQDT(REAL delt, REAL inletTemp);
    
    /// Add a solution to the list
    void AddSolution(TPBRThermalSolution &nextsol);
    
};
#endif
