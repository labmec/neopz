/*
 *  SteamMesh.cpp
 *  PZ
 *
 *  Created by Philippe Devloo on 3/5/11.
 *  Copyright 2011 UNICAMP. All rights reserved.
 *
 */

#include "tpbrsteammesh.h"

#include "tpbrsteamflux.h"
#include "tpbrcellconservation.h"
#include "pzbndmat.h"

#ifdef _AUTODIFF

// ESaturationWater, ESaturationOil, ESaturationSteam, ETemperature, EPressureWater, EPressureSteam, EPressureOil,EPhaseChange

/// initialize the state variables
TPBrSteamMesh::TPBrSteamMesh(int numcells, REAL temperature, REAL pressure, REAL WellRadius, REAL ReservoirRadius, REAL oilsaturation)
{
	fNumCells = numcells;
	fPrevState.Redim(NumEquations(), 1);
	fWellRadius = WellRadius;
	fReservoirRadius = ReservoirRadius;
	SetGeometricProgression(1.5);
	fPrevState(0,0) = pressure;
    fDelt = 1.;
	int ic;
	for (ic=0; ic<numcells; ic++) {
		TPZManVector<int> equation,state;
		CellDestination(ic, equation, state);
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::ESaturationWater],0) = 1.-oilsaturation;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::ESaturationOil],0) = oilsaturation;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::ESaturationSteam],0) = 0.;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::ETemperature],0) = temperature;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::EPressureWater],0) = pressure;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::EPressureOil],0) = pressure;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::EPressureSteam],0) = pressure;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::EPhaseChange],0) = 0.;
	}
	fNextState = fPrevState;
}

int TPBrSteamMesh::NumEquations()
{
	return fNumCells*TPBrCellConservation::NumCellEq+(fNumCells+1)*TPBrSteamFlux::NumFluxEq+TPBrSteamFlux::NumInletVars;
}

/// bandwidth of the global system of equations
int TPBrSteamMesh::Bandwidth()
{
	return (TPBrSteamFlux::NumFluxEq-1)+TPBrCellConservation::NumCellEq;
}


/// extract the state variables of the cell
void TPBrSteamMesh::ExtractCellState(int cell, TPZFMatrix &glob, TPZVec<REAL> &cellstate)
{
	TPZManVector<int> equation,state;
	CellDestination(cell,state,equation);
	int nstate = state.NElements();
    cellstate.Resize(nstate);
	int is;
	for (is=0; is<nstate; is++) {
		cellstate[is] = glob(state[is],0);
	}
}

/// extract the state variables of the interface element
void TPBrSteamMesh::ExtractInterfaceState(int interface, TPZFMatrix &glob, TPZVec<REAL> &interfacestate)
{
	TPZManVector<int> equation,state;
	FluxDestination(interface,state,equation);
	int nstate = state.NElements();
    interfacestate.Resize(nstate);
	int is;
	for (is=0; is<nstate; is++) {
		interfacestate[is] = glob(state[is],0);
	}
}

/// assemble the contribution of the cell
void TPBrSteamMesh::AssembleCell(int cell, TPZFMatrix &ekcell, TPZFMatrix &efcell, TPZMatrix &glob, TPZFMatrix &res)
{
	TPZManVector<int> equation,state;
	CellDestination(cell,equation,state);
	AssembleElement(ekcell,efcell,equation,state,glob,res);
}

/// assemble the contribution of the cell
void TPBrSteamMesh::AssembleInterface(int interface, TPZFMatrix &ekinterface, TPZFMatrix &efinterface, TPZMatrix &glob, TPZFMatrix &res)
{
	TPZManVector<int> equation,state;
	FluxDestination(interface,equation,state);
	AssembleElement(ekinterface,efinterface,equation,state,glob,res);
}

/// abstract assemble procedure
void TPBrSteamMesh::AssembleElement(TPZFMatrix &ek, TPZFMatrix &ef, TPZVec<int> &equation, TPZVec<int> &state, TPZMatrix &glob, TPZFMatrix &res)
{
	int neq = equation.NElements();
	int nstate = state.NElements();
	int ieq;
	int is;
	for (ieq=0; ieq<neq; ieq++) 
	{
		res(equation[ieq],0) += ef(ieq,0);
		for (is=0; is<nstate; is++) {
			glob(equation[ieq],state[is]) += ek(ieq,is);
		}
	}
}

/// equation and state destination indexes for a cell
void TPBrSteamMesh::CellDestination(int cell, TPZVec<int> &equation, TPZVec<int> &state)
{
	int firsteq = TPBrSteamFlux::NumFluxEq*(cell+1)+TPBrSteamFlux::NumInletVars+TPBrCellConservation::NumCellEq*cell;
	int firststate = firsteq-TPBrSteamFlux::NumFluxEq;
	int neq = TPBrCellConservation::NumCellEq;
	int nstate = TPBrSteamFlux::NumFluxEq*2+TPBrCellConservation::NumCellEq;
	if (cell==-1) {
		firsteq = 0;
		firststate = 0;
		neq = TPBrSteamFlux::NumInletVars;
		nstate = TPBrSteamFlux::NumInletVars+TPBrSteamFlux::NumFluxEq+TPBrCellConservation::NumCellEq;
	}
	equation.Resize(neq);
	state.Resize(nstate);
	for (int ieq=0; ieq<neq; ieq++) {
		equation[ieq] = firsteq+ieq;
	}
	for (int is=0; is<nstate; is++) {
		state[is] = firststate+is;
	}
}

/// equation and state destination indexes for an interface
void TPBrSteamMesh::FluxDestination(int interface, TPZVec<int> &equation, TPZVec<int> &state)
{
	int firsteq = TPBrSteamFlux::NumFluxEq*(interface)+TPBrSteamFlux::NumInletVars+TPBrCellConservation::NumCellEq*interface;
	int firststate = firsteq-TPBrCellConservation::NumCellEq;
	int neq = TPBrSteamFlux::NumFluxEq;
	int nstate = TPBrSteamFlux::NumFluxEq+TPBrCellConservation::NumCellEq*2;
	if (interface == 0) {
//		firsteq=0;
		firststate=0;
//		neq = TPBrSteamFlux::NumInletVars+TPBrSteamFlux::NumFluxEq;
		nstate = TPBrSteamFlux::NumInletVars+TPBrSteamFlux::NumFluxEq+TPBrCellConservation::NumCellEq;
	}
	else if(interface == fNumCells) {
		nstate=TPBrSteamFlux::NumFluxEq+TPBrCellConservation::NumCellEq;
	}

	equation.Resize(neq);
	state.Resize(nstate);
	for (int ieq=0; ieq<neq; ieq++) {
		equation[ieq] = firsteq+ieq;
	}
	for (int is=0; is<nstate; is++) {
		state[is] = firststate+is;
	}
}

/// compute the tangent matrix and residual
void TPBrSteamMesh::ComputeTangent(TPZMatrix &tangent, TPZFMatrix &residual)
{
	tangent.Zero();
	residual.Zero();
	TPZManVector<REAL,TPBrSteamFlux::NumInletVars> inlet(TPBrSteamFlux::NumInletVars,0);
	TPZManVector<REAL,TPBrSteamFlux::NumFluxEq> leftface(TPBrSteamFlux::NumFluxEq,0), rightface(TPBrSteamFlux::NumFluxEq,0),face(TPBrSteamFlux::NumFluxEq,0);
	TPZManVector<REAL,TPBrCellConservation::NumCellEq> leftstate(TPBrCellConservation::NumCellEq,0),rightstate(TPBrCellConservation::NumCellEq,0),state(TPBrCellConservation::NumCellEq,0);
	TPZManVector<REAL,TPBrCellConservation::NumCellEq> prevstate(TPBrCellConservation::NumCellEq);
	TPZFMatrix ek, ef;
	REAL delx,area,delt,volume;
	delt = fDelt;
	ExtractCellState(-1, fNextState, inlet);
	ExtractInterfaceState(0, fNextState, face);
	ExtractCellState(0, fNextState, rightstate);
	REAL coord = NodeCoord(0);
	area = 2.*M_PI*coord;
	delx = fFirstCellSize*(1.+1./fGeometricProgression)/2.;
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "before computing the inlet calcstiff\n";
        sout << "inlet state " << inlet << std::endl;
        sout << "face state " << face << std::endl;
        sout << "cell state " << rightstate << std::endl;
    }
#endif
	fInterface.InletCalcStiff(inlet,rightstate,face,delx,area,delt,ek,ef);
	AssembleInterface(0, ek, ef, tangent, residual);
	
	int icell;
	for (icell=0; icell<fNumCells; icell++) {
		ExtractInterfaceState(icell, fNextState, leftface);
		ExtractCellState(icell, fNextState, state);
		ExtractCellState(icell, fPrevState, prevstate);
		ExtractInterfaceState(icell+1, fNextState, rightface);
		REAL leftco = NodeCoord(icell);
		REAL rightco = NodeCoord(icell+1);
		volume = M_PI*(rightco*rightco-leftco*leftco);
		fCell.CalcStiff(leftface, state, rightface, prevstate, volume, delt, ek, ef);
		AssembleCell(icell, ek, ef, tangent, residual);
	}
	int iface;
	for (iface=1; iface<fNumCells; iface++) {
		ExtractCellState(iface, fNextState, leftstate);
		ExtractInterfaceState(iface, fNextState, face);		
		ExtractCellState(iface+1, fNextState, rightstate);
		delx = (CellSize(iface-1)+CellSize(iface))/2.;
		area = M_PI*2.*NodeCoord(iface);
		fInterface.CalcStiff(leftface, rightface, face, delx, area, delt, ek, ef);
		AssembleInterface(iface, ek, ef, tangent, residual);
	}
	
	ExtractCellState(fNumCells-1, fNextState, leftstate);
	ExtractInterfaceState(fNumCells, fNextState, face);
	delx = (CellSize(fNumCells)+CellSize(fNumCells-1))/2.;
	area = M_PI*2.*NodeCoord(fNumCells);
	fInterface.OutletCalcStiff(rightstate,face,delx,area,delt,ek,ef);
	AssembleInterface(fNumCells, ek, ef, tangent, residual);
	
}

/// perform a time step
void TPBrSteamMesh::TimeStep(REAL delt)
{
	fDelt = delt;
	REAL resnorm;
	resnorm = Iterate();
	while (resnorm > 1.e-6) {
		resnorm = Iterate();
	}
	fPrevState = fNextState;
}

/// perform an iteration returning the residual
REAL TPBrSteamMesh::Iterate()
{
	int neq = NumEquations();
	TPZFBMatrix bnd(neq,Bandwidth());
	TPZFMatrix rhs(neq,1);
	ComputeTangent(bnd, rhs);
	REAL residual = Norm(rhs);
	bnd.SolveDirect(rhs,ELU);
	fNextState -= rhs;
	return residual;
}


void TPBrSteamMesh::Print(std::ostream &out)
{
    out << "Mesh for steam flux\n";
    out << "Number of cells " << fNumCells << std::endl;
    out << "Well radius " << fWellRadius << std::endl;
    out << "Reservoir radius " << fReservoirRadius << std::endl;
    out << "Geometric progression " << fGeometricProgression << std::endl;
    out << "First cell size " << fFirstCellSize << std::endl;
    out << "Time step " << fDelt << std::endl;
    out << "Previous state " << fPrevState << std::endl;
    fCell.Print(out);
    fInterface.Print();
    out << "Cell indexing\n";
    int icell;
    for (icell=-1; icell<fNumCells; icell++) {
    	TPZManVector<int> equation,state;
        CellDestination(icell,state,equation);
        out << "Cell " << icell << " Equations " << equation << " state vars " << state << std::endl;
    }
    
    out << "Cell states\n";
    for (icell=-1; icell<fNumCells; icell++) {
        TPZManVector<REAL> cellstate;
        ExtractCellState(icell, fPrevState, cellstate);
        out << "Cell " << icell << " cell state " << cellstate << std::endl;
    }
    out << "Interface indexing\n";
    int iface;
    for(iface = 0; iface <= fNumCells ; iface++)
    {
        TPZManVector<int> equation,state;
        FluxDestination(iface,state,equation);
        out << "Face " << iface << " Equations " << equation << " state vars " << state << std::endl;
    }
    out << "Interface states\n";
    for(iface = 0; iface <= fNumCells ; iface++)
    {
        TPZManVector<REAL> facestate;
        ExtractInterfaceState(iface, fPrevState, facestate);
        out << "Face " << iface << " face state " << facestate << std::endl;
    }
}


#endif