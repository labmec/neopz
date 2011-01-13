/*
 *  tpbrsolutionlist.cpp
 *  IP3D_v5
 *
 *  Created by Philippe Devloo on 1/9/11.
 *  Copyright 2011 UNICAMP. All rights reserved.
 *
 */

#include "tpbrsolutionlist.h"

/// Compute the solution for the next timestep
void TPBRSolutionList::AdvanceSolution(REAL delt, REAL inletTemp)
{
	fDiscretization.SetTimeStep(delt);
	fDiscretization.ComputeStiffness();
	TPZFMatrix nextsol;
	std::list<TPBRThermalSolution>::iterator it = fList.begin();
	while (it != fList.end()) {
		fDiscretization.NextSolution(inletTemp,it->Solution(),nextsol);
		REAL energy = fDiscretization.Energy(nextsol)*it->Area();
		it->UpdateSolution(nextsol,energy);
	}
}
