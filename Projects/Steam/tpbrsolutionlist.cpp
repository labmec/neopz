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
void TPBRSolutionList::AdvanceSolution(REAL delt, REAL inletTemp, REAL &flux, REAL &DQDT)
{
	fDiscretization.SetTimeStep(delt);
	fDiscretization.ComputeStiffness();
	TPZFMatrix nextsol;
	flux = 0.;
	DQDT = 0.;
	REAL localDQDT = fDiscretization.DQDT();
	std::list<TPBRThermalSolution>::iterator it = fList.begin();
	while (it != fList.end()) {
		REAL localFlux;
		REAL area = it->Area();
		fDiscretization.NextSolution(inletTemp,it->Solution(),nextsol, localFlux);
		REAL energy = fDiscretization.Energy(nextsol)*area;
		it->UpdateSolution(nextsol,energy);
		flux += area*localFlux;
		DQDT += area*localDQDT;
		
	}
}
