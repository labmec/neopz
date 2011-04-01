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
void TPBRSolutionList::AdvanceSolution(REAL delt, REAL inletTemp, REAL &flux, REAL &DQDT, bool storesolution)
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
        if(storesolution)
        {
            REAL energy = fDiscretization.Energy(nextsol)*area;
            it->UpdateSolution(nextsol,energy);
            it++;
        }
		flux += area*localFlux;
		DQDT += area*localDQDT;
		
	}
}

/// Add a solution to the list
void TPBRSolutionList::AddSolution(TPBRThermalSolution &nextsol)
{
    fDiscretization.InitializeSolution(nextsol.Solution());
    fList.push_back(nextsol);
}

/// total energy of the solution list
REAL TPBRSolutionList::Energy()
{
    REAL varEnergy = 0.;
    std::list<TPBRThermalSolution>::iterator it;
    for (it=fList.begin(); it != fList.end(); it++) 
    {
        varEnergy += it->Energy();
    }
    return varEnergy;
}

