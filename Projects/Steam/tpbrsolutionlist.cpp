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
void TPBRSolutionList::AdvanceAllSolutions(REAL delt, REAL inletTemp, REAL &flux, REAL &DQDT, bool storesolution)
{
	fDiscretization.SetTimeStep(delt);
	fDiscretization.ComputeStiffness();
    fDelt = delt;
	TPZFMatrix<REAL> nextsol;
	flux = 0.;
	DQDT = 0.;
	REAL localDQDT = fDiscretization.DQDT();
	std::vector<TPBRThermalSolution>::iterator it = fList.begin();
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

/// Set the timestep
void TPBRSolutionList::SetDelt(REAL delt)
{
    if(fDelt != delt)
    {
        fDiscretization.SetTimeStep(delt);
    }
    fDelt = delt;
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
    std::vector<TPBRThermalSolution>::iterator it;
    for (it=fList.begin(); it != fList.end(); it++) 
    {
        varEnergy += it->Energy();
    }
    return varEnergy;
}

/// total energy of the solution list
REAL TPBRSolutionList::Energy(int icell)
{
    return fList[icell].Energy();
}

/// Compute the solution for the next timestep
void TPBRSolutionList::AdvanceSolution(int icell, REAL inletTemp, REAL &flux, REAL &DQDT, bool storesolution)
{
	TPZFMatrix<REAL> nextsol;
	flux = 0.;
	DQDT = 0.;
	REAL localDQDT = fDiscretization.DQDT();
    REAL localFlux;
    REAL area = fList[icell].Area();
    fDiscretization.NextSolution(inletTemp,fList[icell].Solution(),nextsol, localFlux);
    if(storesolution)
    {
        REAL energy = fDiscretization.Energy(nextsol)*area;
        fList[icell].UpdateSolution(nextsol,energy);
	}
    flux += area*localFlux;
    DQDT += area*localDQDT;		
    
}

/// Compute the variation of the flux with respect to the inlet temperature
REAL TPBRSolutionList::DQDT(int icell, REAL inletTemp, REAL &flux)
{
	TPZFNMatrix<201> nextsol;
	flux = 0.;
	REAL dfluxdt = 0.;
	REAL localDQDT = fDiscretization.DQDT();
    REAL localFlux;
    REAL area = fList[icell].Area();
    fDiscretization.NextSolution(inletTemp,fList[icell].Solution(),nextsol, localFlux);
    flux += area*localFlux;
    dfluxdt += area*localDQDT;
	return dfluxdt;	
    
}


