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
	TPBRThermalDisc fDiscretization;

public:
	/// total energy of the solution list
	REAL Energy()
	{
		REAL varEnergy = 0.;
		std::list<TPBRThermalSolution>::iterator it = fList.begin();
		while (it != fList.end()) {
			varEnergy += it->Energy();
		}
		return varEnergy;
	}
	/// Compute the solution for the next timestep
	void AdvanceSolution(REAL delt, REAL inletTemp, REAL &flux, REAL &DQDT);
	
	REAL DQDT(REAL delt, REAL inletTemp);
};
#endif
