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
#include "PropertiesTable.h"
#include "pzbndmat.h"
#include "pzgengrid.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("br.steammesh"));
#endif

#ifdef _AUTODIFF

//WaterDataInStateOfSaturation waterdata;

// ESaturationWater, ESaturationOil, ESaturationSteam, ETemperature, EPressureWater, EPressureSteam, EPressureOil,EPhaseChange

/// initialize the state variables
TPBrSteamMesh::TPBrSteamMesh(int numcells, REAL temperature, REAL pressure, REAL WellRadius, REAL ReservoirRadius, REAL reservoirheight, REAL oilsaturation) :
    fHeight(reservoirheight)
{
	fNumCells = numcells;
	fPrevState.Redim(NumEquations(), 1);
	fWellRadius = WellRadius;
	fReservoirRadius = ReservoirRadius;
    if(fNumCells > 1)
    {
        REAL progression = TPZGenGrid::GeometricProgression(fWellRadius, ReservoirRadius-fWellRadius, fNumCells);
        SetGeometricProgression(progression);        
    }
    else
    {    
        SetGeometricProgression(1.);
    }
	fPrevState(TPBrSteamFlux::EInletPressure,0) = 0;
    fPrevState(TPBrSteamFlux::EInletTemperature,0) = temperature;
    fDelt = 1.;
	int ic;
	for (ic=0; ic<numcells; ic++) {
		TPZManVector<int> equation,state;
		CellDestination(ic, equation, state);
        LOGPZ_DEBUG(logger, state[TPBrSteamFlux::EDarcyVelocityWater])
        fPrevState(state[TPBrSteamFlux::EDarcyVelocityWater],0) = 1.;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::ESaturationWater],0) = 1.-oilsaturation;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::ESaturationOil],0) = oilsaturation;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::ESaturationSteam],0) = 0.;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::ETemperature],0) = temperature;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::EPressureWater],0) = 0;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::EPressureOil],0) = 0;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::EPressureSteam],0) = 0;
		fPrevState(state[TPBrSteamFlux::NumFluxEq+TPBrCellConservation::EPhaseChange],0) = 0.;
	}
	fNextState = fPrevState;
    TPBrSteamFlux::fFarfieldPressureOil = pressure;
	TPBrSteamFlux::fFarfieldPressureWater = pressure;
	TPBrSteamFlux::fFarfieldPressureSteam = pressure;
	TPBrSteamFlux::fFarfieldTemperature = temperature;
	TPBrSteamFlux::fFarfieldSaturationOil = oilsaturation;
	TPBrSteamFlux::fFarfieldSaturationWater = 1.-oilsaturation;
	TPBrSteamFlux::fFarfieldSaturationSteam = 0.;
    
    TPBrScales::fReferencePressure = pressure;
	
    TPBRThermalDiscretization disc(fLayer.fHeight, fLayer.fNumElements, fLayer.fThermalCapacity, fLayer.fConductivity, fLayer.fDensity, temperature);
    fThermal.SetDiscretization(disc);
    InitializeThermalProblem();
    fThermal.SetDelt(fDelt);

}

int TPBrSteamMesh::NumEquations()
{
    // each cell and flux contribute to an independent equation
	return fNumCells*TPBrCellConservation::NumCellEq+(fNumCells+1)*TPBrSteamFlux::NumFluxEq+TPBrSteamFlux::NumInletVars;
}

/// bandwidth of the global system of equations
int TPBrSteamMesh::Bandwidth()
{
	return TPBrSteamFlux::NumInletVars+ 2*(TPBrSteamFlux::NumFluxEq)+2*TPBrCellConservation::NumCellEq;
}


/// extract the state variables of the cell
void TPBrSteamMesh::ExtractCellState(int cell, TPZFMatrix<REAL> &glob, TPZVec<REAL> &cellstate)
{
	TPZManVector<int> equation,state;
	CellDestination(cell,equation,state);
	int nstate = equation.NElements();
    cellstate.Resize(nstate);
	int is;
	for (is=0; is<nstate; is++) {
		cellstate[is] = glob(equation[is],0);
	}
}

/// extract the state variables of the interface element
void TPBrSteamMesh::ExtractInterfaceState(int interface, TPZFMatrix<REAL> &glob, TPZVec<REAL> &interfacestate)
{
	TPZManVector<int> equation,state;
	FluxDestination(interface,equation,state);
	int nstate = equation.NElements();
    interfacestate.Resize(nstate);
	int is;
	for (is=0; is<nstate; is++) {
		interfacestate[is] = glob(equation[is],0);
	}
}

/// extract the state variables of the interface element
void TPBrSteamMesh::ExtractInletfaceState(TPZFMatrix<REAL> &glob, TPZVec<REAL> &interfacestate)
{
	TPZManVector<int> equation,state;
	FluxDestinationInlet(equation,state);
	int nstate = equation.NElements();
    interfacestate.Resize(nstate);
	int is;
	for (is=0; is<nstate; is++) {
		interfacestate[is] = glob(equation[is],0);
	}
}

/// assemble the contribution of the cell
void TPBrSteamMesh::AssembleCell(int cell, TPZFMatrix<REAL> &ekcell, TPZFMatrix<REAL> &efcell, TPZMatrix<REAL> &glob, TPZFMatrix<REAL> &res,
                                 TPZVec<REAL> &statescales)
{
	TPZManVector<int> equation,state;
	CellDestination(cell,equation,state);
    TPZManVector<REAL> eqscale,cellstatescale;
    TPBrCellConservation::Scales(eqscale, cellstatescale);
    int neq = ekcell.Rows();
    int nstate = ekcell.Cols();
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "eqscale " << eqscale << std::endl;
        sout << "statescale " << cellstatescale << std::endl;
        sout << "equations of cell " << equation << std::endl;
        sout << "state indices " << state << std::endl;
        sout << "globscale ";
        for (int jst=0; jst<nstate; jst++) {
            sout << statescales[state[jst]] << " ";
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    for(int ieq=0; ieq<neq; ieq++)
    {
        for(int jst=0; jst<nstate; jst++)
        {
            ekcell(ieq,jst) /= eqscale[ieq];
            ekcell(ieq,jst) *= statescales[state[jst]];
        }
        efcell(ieq,0) /= eqscale[ieq];
    }
#ifdef LOG4CXX
    {
        std::stringstream sout;
        ekcell.Print("Cell Stiffness",sout);
        efcell.Print("Cell Rhs ", sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
	AssembleElement(ekcell,efcell,equation,state,glob,res);
}

/// assemble the contribution of the cell
void TPBrSteamMesh::AssembleInterface(int interface, TPZFMatrix<REAL> &ekinterface, TPZFMatrix<REAL> &efinterface, TPZMatrix<REAL> &glob, TPZFMatrix<REAL> &res, TPZVec<REAL> &statescales)
{
	TPZManVector<int> equation,state;
    if(interface == 0)
    {
        FluxDestinationInlet(equation, state);
    }
    else
    {
        FluxDestination(interface,equation,state);
    }
    TPZManVector<REAL> eqscales,interfacestatescale;
    if(interface == 0)
    {
        TPBrSteamFlux::InletScales(eqscales, interfacestatescale);
    }
    else
    {
        TPBrSteamFlux::Scales(eqscales, interfacestatescale);
    }
    int neq = ekinterface.Rows();
    int nstate = ekinterface.Cols();
    for(int ieq=0; ieq<neq; ieq++)
    {
        for(int jst=0; jst<nstate; jst++)
        {
            ekinterface(ieq,jst) /= eqscales[ieq];
            ekinterface(ieq,jst) *= statescales[state[jst]];
        }
        efinterface(ieq,0) /= eqscales[ieq];
    }
#ifdef LOG4CXX
    {
        std::stringstream sout;
        ekinterface.Print("Interface Stiffness",sout);
        efinterface.Print("Interface Rhs ", sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
	AssembleElement(ekinterface,efinterface,equation,state,glob,res);
}

/// abstract assemble procedure
void TPBrSteamMesh::AssembleElement(TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZVec<int> &equation, TPZVec<int> &state, TPZMatrix<REAL> &glob, TPZFMatrix<REAL> &res)
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
	int firsteq = FirstCellEquation(cell);
	int firststate = FirstInterfaceEquation(cell);
	int neq = TPBrCellConservation::NumCellEq;
	int nstate = TPBrSteamFlux::NumFluxEq*2+TPBrCellConservation::NumCellEq;
	if (cell==-1) {
        DebugStop();
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
    int eq = FirstInterfaceEquation(cell);
    int counter = 0;
    int is;
	for (is=0; is<TPBrSteamFlux::NumFluxEq; is++,counter++) {
		state[is] = eq+counter;
	}
    eq = FirstCellEquation(cell);
    counter = 0;
    for (; is<TPBrCellConservation::NumCellEq+TPBrSteamFlux::NumFluxEq; is++,counter++) {
        state[is] = eq+counter;
    }
    eq = FirstInterfaceEquation(cell+1);
    counter = 0;
    for (; is<TPBrCellConservation::NumCellEq+2*TPBrSteamFlux::NumFluxEq; is++,counter++) {
        state[is] = eq+counter;
    }
}

/// equation and state destination indexes for an interface
void TPBrSteamMesh::FluxDestination(int interface, TPZVec<int> &equation, TPZVec<int> &state)
{
	int firsteq = FirstInterfaceEquation(interface);
	int neq = TPBrSteamFlux::NumFluxEq;
	int nstate = TPBrSteamFlux::NumFluxEq+TPBrCellConservation::NumCellEq*2;
    if(interface == 0)
    {
        nstate = 0;
    }
    else if(interface == fNumCells) {
		nstate=TPBrSteamFlux::NumFluxEq+TPBrCellConservation::NumCellEq;
	}

	equation.Resize(neq);
	state.Resize(nstate);
	for (int ieq=0; ieq<neq; ieq++) {
		equation[ieq] = firsteq+ieq;
	}
    if(interface == 0) return;
    int eq = FirstCellEquation(interface-1);
    int counter = 0;
    int is;
	for (is=0; is<TPBrCellConservation::NumCellEq; is++,counter++) {
		state[is] = eq+counter;
	}
    eq = FirstInterfaceEquation(interface);
    counter = 0;
    for (; is<TPBrCellConservation::NumCellEq+TPBrSteamFlux::NumFluxEq; is++,counter++) {
        state[is] = eq+counter;
    }
    if(interface != fNumCells)
    {
        eq = FirstCellEquation(interface);
        counter = 0;
        for (; is<2*TPBrCellConservation::NumCellEq+TPBrSteamFlux::NumFluxEq; is++,counter++) {
            state[is] = eq+counter;
        }
    }
}

/// equation and state destination indexes for an interface
void TPBrSteamMesh::FluxDestinationInlet(TPZVec<int> &equation, TPZVec<int> &state)
{
	int firsteq = 0;
	int firststate = 0;
	int neq = TPBrSteamFlux::NumInletVars+TPBrSteamFlux::NumFluxEq;
    int nstate = TPBrSteamFlux::NumInletVars+TPBrSteamFlux::NumFluxEq+TPBrCellConservation::NumCellEq;
    
	equation.Resize(neq);
	state.Resize(nstate);
	for (int ieq=0; ieq<neq; ieq++) {
		equation[ieq] = firsteq+ieq;
	}
    int is;
	for (is=0; is<TPBrSteamFlux::NumInletVars+TPBrSteamFlux::NumFluxEq; is++) {
		state[is] = firststate+is;
	}
    int eq = FirstCellEquation(0);
    int counter = 0;
	for (; is<nstate; is++,counter++) {
		state[is] = eq+counter;
	}

}

/// compute the tangent matrix and residual
void TPBrSteamMesh::ComputeTangent(TPZMatrix<REAL> &tangent, TPZFMatrix<REAL> &residual, TPZVec<REAL> &statescales)
{
	tangent.Zero();
	residual.Zero();
	TPZManVector<REAL,TPBrSteamFlux::NumInletVars> inlet(TPBrSteamFlux::NumInletVars,0);
	TPZManVector<REAL,TPBrSteamFlux::NumFluxEq+TPBrSteamFlux::NumInletVars> leftface(TPBrSteamFlux::NumFluxEq,0), rightface(TPBrSteamFlux::NumFluxEq,0),face(TPBrSteamFlux::NumFluxEq+TPBrSteamFlux::NumInletVars,0);
	TPZManVector<REAL,TPBrCellConservation::NumCellEq> leftstate(TPBrCellConservation::NumCellEq,0),rightstate(TPBrCellConservation::NumCellEq,0),state(TPBrCellConservation::NumCellEq,0);
	TPZManVector<REAL,TPBrCellConservation::NumCellEq> prevstate(TPBrCellConservation::NumCellEq);
	TPZFMatrix<REAL> ek, ef;
	REAL delx,area,delt,volume;
	delt = fDelt;
//	ExtractCellState(-1, fNextState, inlet);
	ExtractInletfaceState(fNextState, face);
	ExtractCellState(0, fNextState, rightstate);
	REAL coord = NodeCoord(0);
	area = 2.*M_PI*coord*fHeight;
	delx = fFirstCellSize*(1.+1./fGeometricProgression)/2.;
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "before computing the inlet calcstiff\n";
//        sout << "inlet state " << inlet << std::endl;
        sout << "face state " << face << std::endl;
        sout << "cell state " << rightstate << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	fInterface.InletCalcStiff(rightstate,face,delx,area,delt,ek,ef);
	AssembleInterface(0, ek, ef, tangent, residual,statescales);
	
	int icell;
	for (icell=0; icell<fNumCells; icell++) {
		ExtractInterfaceState(icell, fNextState, leftface);
		ExtractCellState(icell, fNextState, state);
		ExtractCellState(icell, fPrevState, prevstate);
		ExtractInterfaceState(icell+1, fNextState, rightface);
		REAL leftco = NodeCoord(icell);
		REAL rightco = NodeCoord(icell+1);
		volume = M_PI*(rightco*rightco-leftco*leftco)*fHeight;
		fCell.CalcStiff(leftface, state, rightface, prevstate, volume, delt, ek, ef);

        // contribute the flux of the impermeous layers
        REAL inletTemp = state[TPBrCellConservation::ETemperature];
        REAL Flux;
        REAL dfluxdT = fThermal.DQDT(icell, inletTemp, Flux);
        ef(TPBrCellConservation::EEnergyCons,0) -= 2.*Flux;
        ek(TPBrCellConservation::EEnergyCons,TPBrSteamFlux::NumFluxEq+TPBrCellConservation::ETemperature) -= 2.*dfluxdT;

        
        
		AssembleCell(icell, ek, ef, tangent, residual,statescales);
	}
	int iface;
	for (iface=1; iface<fNumCells; iface++) {
		ExtractCellState(iface-1, fNextState, leftstate);
		ExtractInterfaceState(iface, fNextState, face);		
		ExtractCellState(iface, fNextState, rightstate);
		delx = (CellSize(iface-1)+CellSize(iface))/2.;
		area = M_PI*2.*NodeCoord(iface)*fHeight;
		fInterface.CalcStiff(leftstate, rightstate, face, delx, area, delt, ek, ef);
		AssembleInterface(iface, ek, ef, tangent, residual,statescales);
	}
	
	ExtractCellState(fNumCells-1, fNextState, leftstate);
	ExtractInterfaceState(fNumCells, fNextState, face);
	delx = (CellSize(fNumCells)+CellSize(fNumCells-1))/2.;
	area = M_PI*2.*NodeCoord(fNumCells)*fHeight;
	fInterface.OutletCalcStiff(rightstate,face,delx,area,delt,ek,ef);
	AssembleInterface(fNumCells, ek, ef, tangent, residual,statescales);
	
}

/// compute the tangent matrix and residual
void TPBrSteamMesh::ComputeTangent2(TPZMatrix<REAL> &tangent, TPZFMatrix<REAL> &residual, TPZVec<REAL> &statescales)
{
	tangent.Zero();
	residual.Zero();
	TPZManVector<REAL,TPBrSteamFlux::NumInletVars> inlet(TPBrSteamFlux::NumInletVars,0);
	TPZManVector<REAL,TPBrSteamFlux::NumFluxEq+TPBrSteamFlux::NumInletVars> leftface(TPBrSteamFlux::NumFluxEq,0), rightface(TPBrSteamFlux::NumFluxEq,0),face(TPBrSteamFlux::NumFluxEq+TPBrSteamFlux::NumInletVars,0);
	TPZManVector<REAL,TPBrCellConservation::NumCellEq> leftstate(TPBrCellConservation::NumCellEq,0),rightstate(TPBrCellConservation::NumCellEq,0),state(TPBrCellConservation::NumCellEq,0);
	TPZManVector<REAL,TPBrCellConservation::NumCellEq> prevstate(TPBrCellConservation::NumCellEq);
	TPZFMatrix<REAL> ek, ef;
	REAL delx,area,delt,volume;
	delt = fDelt;
    //	ExtractCellState(-1, fNextState, inlet);
	ExtractInletfaceState(fNextState, face);
	ExtractCellState(0, fNextState, rightstate);
	REAL coord = NodeCoord(0);
	area = 2.*M_PI*coord*fHeight;
	delx = fFirstCellSize*(1.+1./fGeometricProgression)/2.;
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "before computing the inlet calcstiff\n";
        //        sout << "inlet state " << inlet << std::endl;
        sout << "face state " << face << std::endl;
        sout << "cell state " << rightstate << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
	fInterface.InletCalcStiff(rightstate,face,delx,area,delt,ek,ef);
	AssembleInterface(0, ek, ef, tangent, residual,statescales);
	
    
	int icell;
	for (icell=0; icell<fNumCells; icell++) {
        // ESTE ESTA DANDO ERRADO!!
		ExtractInterfaceState(icell, fNextState, leftface);
		ExtractCellState(icell, fNextState, state);
		ExtractCellState(icell, fPrevState, prevstate);
		ExtractInterfaceState(icell+1, fNextState, rightface);
		REAL leftco = NodeCoord(icell);
		REAL rightco = NodeCoord(icell+1);
		volume = M_PI*(rightco*rightco-leftco*leftco)*fHeight;
		fCell.CalcStiff(leftface, state, rightface, prevstate, volume, delt, ek, ef);
        
//        ek.Zero();
//        ef.Zero();

        // contribute the flux of the impermeous layers
        REAL inletTemp = state[TPBrCellConservation::ETemperature];
        REAL Flux;
        REAL dfluxdT = fThermal.DQDT(icell, inletTemp, Flux);
        ef(TPBrCellConservation::EEnergyCons,0) -= 2.*Flux;
        ek(TPBrCellConservation::EEnergyCons,TPBrSteamFlux::NumFluxEq+TPBrCellConservation::ETemperature) -= 2.*dfluxdT;
		
        AssembleCell(icell, ek, ef, tangent, residual,statescales);
	}
	int iface;
	for (iface=1; iface<fNumCells; iface++) {
		ExtractCellState(iface-1, fNextState, leftstate);
		ExtractInterfaceState(iface, fNextState, face);		
		ExtractCellState(iface, fNextState, rightstate);
		delx = (CellSize(iface-1)+CellSize(iface))/2.;
		area = M_PI*2.*NodeCoord(iface)*fHeight;
		fInterface.CalcStiff(leftstate, rightstate, face, delx, area, delt, ek, ef);
		AssembleInterface(iface, ek, ef, tangent, residual,statescales);
	}
	
	ExtractCellState(fNumCells-1, fNextState, leftstate);
	ExtractInterfaceState(fNumCells, fNextState, face);
	delx = (CellSize(fNumCells)+CellSize(fNumCells-1))/2.;
	area = M_PI*2.*NodeCoord(fNumCells)*fHeight;
	fInterface.OutletCalcStiff(rightstate,face,delx,area,delt,ek,ef);
	AssembleInterface(fNumCells, ek, ef, tangent, residual,statescales);
}

/// perform a time step
void TPBrSteamMesh::TimeStep(REAL delt)
{
	fDelt = delt;
    fThermal.SetDelt(fDelt);
	REAL resnorm;
	resnorm = Iterate();
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Residual norm " << resnorm;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	while (resnorm > 1.e-5) {
		resnorm = Iterate();
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "Residual norm " << resnorm;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
	}
    REAL leftin,rightout,energybalance;
    REAL energyconfinementlayerbefore = fThermal.Energy();
    UpdateConfinementLayer();
    REAL energyconfinementlayerafter = fThermal.Energy();
    
    EnergyBalance(leftin, rightout, energybalance );
    std::cout << "Energy input " << leftin << " Energy out " << rightout << " Energy balance " << energybalance <<
    " saldo " << leftin-rightout-energybalance << std::endl;
    std::cout << "Energy transferred to the confinement layer " << 2.*(energyconfinementlayerafter-energyconfinementlayerbefore) << std::endl;
    
	fPrevState = fNextState;
}

/// perform an iteration returning the residual
REAL TPBrSteamMesh::Iterate()
{
	int neq = NumEquations();
	TPZFBMatrix<REAL> bnd(neq,Bandwidth());
	TPZFMatrix<REAL> rhs(neq,1);
    TPZStack<REAL> scales;
    StateScales(scales);
	ComputeTangent(bnd, rhs,scales);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        rhs.Print("rhs after ComputeTangent",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	REAL residual = Norm(rhs);
    std::list<int> singular;
	bnd.SolveDirect(rhs,ELU,singular);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        rhs.Print("solution correction before scaling",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int ieq;
    for (ieq=0; ieq<neq; ieq++) {
        rhs(ieq,0) *= -scales[ieq];
    }
#ifdef LOG4CXX
    {
        std::stringstream sout;
        rhs.Print("solution correction",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    REAL scale = LimitRange(fNextState, rhs);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "scale factor applied to the solution correction = " << scale;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    rhs *= scale;
    if(fDelt > 0.) 
    {
//       VerifyTangent(rhs); 
    }
    
	fNextState += rhs;
#ifdef LOG4CXX
    {
        std::stringstream sout;
        fNextState.Print("solution after correction",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    ProjectSolution(fNextState);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        fNextState.Print("solution after projection",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    // put the pressure of the inlet larger than the pressure in the cell
    REAL inletpress = fNextState(TPBrSteamFlux::EInletPressure,0);
    TPZManVector<int> eqs,sts;
    CellDestination(0, eqs, sts);
    int cellpressindex = eqs[TPBrCellConservation::EPressureWater];
    REAL cellpressure = fNextState(cellpressindex,0);
    if(inletpress < cellpressure)
    {
        inletpress = cellpressure + 10000.;
        fNextState(TPBrSteamFlux::EInletPressure,0) = inletpress;
    }
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
    for (icell=0; icell<fNumCells; icell++) {
    	TPZManVector<int> equation,state;
        CellDestination(icell,equation,state);
        out << "Cell " << icell << " Equations " << equation << "\nstate vars " << state << std::endl;
    }
    
    out << "Cell states\n";
    for (icell=0; icell<fNumCells; icell++) {
        TPZManVector<REAL> cellstate;
        ExtractCellState(icell, fPrevState, cellstate);
        out << "Cell " << icell << " cell state " << cellstate << std::endl;
    }
    out << "Interface indexing\n";
    int iface;
    for(iface = 0; iface <= fNumCells ; iface++)
    {
        TPZManVector<int> equation,state;
        FluxDestination(iface,equation,state);
        out << "Face " << iface << " Equations " << equation << "\nstate vars " << state << std::endl;
    }
    out << "Interface states\n";
    for(iface = 0; iface <= fNumCells ; iface++)
    {
        TPZManVector<REAL> facestate;
        if(iface == 0) 
        {
            ExtractInletfaceState(fPrevState, facestate);
        }
        else
        {
            ExtractInterfaceState(iface, fPrevState, facestate);            
        }
        out << "Face " << iface << " face state " << facestate << std::endl;
    }
}

/// compute the scales associated with the state variables
void TPBrSteamMesh::StateScales(TPZStack<REAL> &scales)
{
    scales.Resize(NumEquations(), 0);
    TPZManVector<REAL> eqscales,statescale;
    TPBrSteamFlux::InletScales(eqscales, statescale);
    int nst = statescale.NElements();
    int ist;
    for (ist=0; ist<nst; ist++) {
        scales[ist] = statescale[ist];
    }
    int ic;
    for (ic=0; ic<fNumCells; ic++) {
        TPBrCellConservation::Scales(eqscales, statescale);
        TPZManVector<int> equations,statevar;
        CellDestination(ic, equations, statevar);
        int ist;
        int nst = equations.NElements();
        for (ist=0; ist<nst; ist++) {
            scales[equations[ist]] = statescale[ist];
        }
        TPBrSteamFlux::Scales(eqscales, statescale);
        nst = statescale.NElements();
        FluxDestination(ic+1, equations, statevar);
        for (ist=0; ist<nst; ist++) {
            scales[equations[ist]] = statescale[ist];
        }
    }
}

/// initialize the mass fluxes considering
void TPBrSteamMesh::SetMassFlux(REAL waterflux, REAL oilflux)
{
    int face;
    for(face=0; face<=fNumCells; face++)
    {
        TPZManVector<int> equations, state;
        FluxDestination(face, equations, state);
        int eq0 = equations[0];
        int eqwaterflux = eq0+TPBrSteamFlux::EMassFluxWater;
        fPrevState(eqwaterflux,0) = waterflux;
        fNextState(eqwaterflux,0) = waterflux;
        int eqoilflux = eq0+TPBrSteamFlux::EMassFluxOil;
        fPrevState(eqoilflux,0) = oilflux;
        fNextState(eqoilflux,0) = oilflux;
    }
}

/// initialize the cell pressure
void TPBrSteamMesh::SetCellPressure(int cell, REAL pressure)
{
    TPZManVector<int> equations,state;
    CellDestination(cell, equations, state);
    int eq0 = equations[0];
    int eqpress = eq0+TPBrCellConservation::EPressureWater;
    fPrevState(eqpress,0) = pressure;
    fNextState(eqpress,0) = pressure;
    eqpress = eq0+TPBrCellConservation::EPressureOil;
    fPrevState(eqpress,0) = pressure;
    fNextState(eqpress,0) = pressure;
    eqpress = eq0+TPBrCellConservation::EPressureSteam;
    fPrevState(eqpress,0) = pressure;
    fNextState(eqpress,0) = pressure;
}

/// initialize the temperature
void TPBrSteamMesh::SetTemperature(REAL temperature)
{
    int cell;
    for (cell = 0; cell < fNumCells; cell++) {
        TPZManVector<int> equations,state;
        CellDestination(cell, equations, state);
        int eq0 = equations[0];
        int eqtemp = eq0+TPBrCellConservation::ETemperature;
        fPrevState(eqtemp,0) = temperature;
        fNextState(eqtemp,0) = temperature;

    }
    
}

/// initialize the saturation
void TPBrSteamMesh::SetWaterSaturation(REAL saturation)
{
    int cell;
    for (cell = 0; cell < fNumCells; cell++) {
        TPZManVector<int> equations,state;
        CellDestination(cell, equations, state);
        int eq0 = equations[0];
        int eqsat = eq0+TPBrCellConservation::ESaturationWater;
        fPrevState(eqsat,0) = saturation;
        fNextState(eqsat,0) = saturation;
        eqsat = eq0+TPBrCellConservation::ESaturationOil;
        fPrevState(eqsat,0) = 1.-saturation;
        fNextState(eqsat,0) = 1.-saturation;
        eqsat = eq0+TPBrCellConservation::ESaturationSteam;
        fPrevState(eqsat,0) = 0.;
        fNextState(eqsat,0) = 0.;
    }    
}


/// initialize the injection data
void TPBrSteamMesh::SetWaterInjection(REAL massflux, REAL temperature)
{
    TPBrSteamFlux::fInletMassFlux = massflux;
    SetMassFlux(massflux, 0.);
    SetTemperature(temperature);
    SetWaterSaturation(1.);
    REAL energyflux = TPBrSteamFlux::EnthalpyWater(temperature)*massflux;
    TPBrSteamFlux::fInletEnergyFlux = energyflux;
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Enthalpy water " << TPBrSteamFlux::EnthalpyWater(temperature) << std::endl;
        sout << "mass flux " << massflux << std::endl;
        sout << "energy flux " << energyflux;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int iface;
    for(iface=0; iface<=fNumCells; iface++)
    {
        TPZManVector<int> equations, state;
        FluxDestination(iface, equations, state);
        int eq0 = equations[0];
        int eqenergyflux = eq0+TPBrSteamFlux::EEnergyFlux;
        fPrevState(eqenergyflux,0) = energyflux;
        fNextState(eqenergyflux,0) = energyflux;
    }
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        Print(sout);
        fNextState.Print("Next State",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
	TPZManVector<REAL,TPBrSteamFlux::NumInletVars> inlet(TPBrSteamFlux::NumInletVars,0);
	TPZManVector<REAL,TPBrSteamFlux::NumFluxEq> leftface(TPBrSteamFlux::NumFluxEq,0), rightface(TPBrSteamFlux::NumFluxEq,0),face(TPBrSteamFlux::NumFluxEq,0);
	TPZManVector<REAL,TPBrCellConservation::NumCellEq> leftstate(TPBrCellConservation::NumCellEq,0),rightstate(TPBrCellConservation::NumCellEq,0),state(TPBrCellConservation::NumCellEq,0);
	TPZManVector<REAL,TPBrCellConservation::NumCellEq> prevstate(TPBrCellConservation::NumCellEq);
	TPZFNMatrix<800> ek, ef;
	REAL delx,area,delt;
	delt = fDelt;
    
	
	ExtractCellState(fNumCells-1, fNextState, leftstate);
	ExtractInterfaceState(fNumCells, fNextState, face);
	delx = (CellSize(fNumCells)+CellSize(fNumCells-1))/2.;
	area = M_PI*2.*NodeCoord(fNumCells)*fHeight;
	fInterface.OutletCalcStiff(leftstate,face,delx,area,delt,ek,ef);
    
    int eqmassconswater = TPBrSteamFlux::EMassFluxWaterEq;
    int statepresswater = TPBrCellConservation::EPressureWater;
    REAL delp = ef(eqmassconswater)/ek(eqmassconswater,statepresswater);
    TPZManVector<int> statevars,equations;
    CellDestination(fNumCells-1, equations, statevars);
    int eqpress = equations[TPBrCellConservation::EPressureWater];
    REAL pressure = fNextState(eqpress,0);
    pressure -= delp;
    SetCellPressure(fNumCells-1, pressure);

    for (iface=fNumCells-1; iface>0; iface--) {
		ExtractCellState(iface-1, fNextState, leftstate);
		ExtractInterfaceState(iface, fNextState, face);		
		ExtractCellState(iface, fNextState, rightstate);
		delx = (CellSize(iface-1)+CellSize(iface))/2.;
		area = M_PI*2.*NodeCoord(iface)*fHeight;
		fInterface.CalcStiff(leftstate, rightstate, face, delx, area, delt, ek, ef);
        REAL efval = ef(eqmassconswater,0);
        REAL ekval = ek(eqmassconswater,statepresswater);
        REAL delp = efval/ekval;
        CellDestination(iface-1, equations, statevars);
        int eqpress = equations[TPBrCellConservation::EPressureWater];
        pressure = fNextState(eqpress,0);
        pressure -= delp;
        SetCellPressure(iface-1, pressure);
	}
    
    //	ExtractCellState(-1, fNextState, inlet);
    eqmassconswater = TPBrSteamFlux::NumInletVars+TPBrSteamFlux::EMassFluxWaterEq;
    statepresswater = TPBrSteamFlux::EInletPressure;
    eqpress = statepresswater;
    fNextState(eqpress,0) = pressure+100.;
	ExtractInletfaceState(fNextState, face);
	ExtractCellState(0, fNextState, rightstate);
	REAL coord = NodeCoord(0);
	area = 2.*M_PI*coord*fHeight;
	delx = fFirstCellSize*(1.+1./fGeometricProgression)/2.;
	fInterface.InletCalcStiff(rightstate,face,delx,area,delt,ek,ef);
    REAL efval = ef(eqmassconswater,0);
    REAL ekval = ek(eqmassconswater,statepresswater);
    delp = efval/ekval;
    pressure = fNextState(eqpress,0);
    pressure -= delp;
    fNextState(eqpress,0) = pressure;
    fPrevState(eqpress,0) = pressure;
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        fNextState.Print("Initialized state",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

/// initialize the steam quality for the inlet pressure
void TPBrSteamMesh::SetSteamQuality(REAL quality, REAL referencepressure)
{
    REAL massflux = TPBrSteamFlux::fInletMassFlux;
    REAL temperature = TPBrCellConservation::TemperatureSaturation(referencepressure);
    REAL enthalpywater  = waterdata.getSaturationStateSpecificEnthalpyToLiquidWater(temperature);
    REAL enthalpysteam = waterdata.getSaturationStateSpecificEnthalpyToSteam(temperature);
    REAL masswater = (1.-quality)*massflux;
    REAL masssteam = quality*massflux;
    REAL energyflux = enthalpywater*masswater+enthalpysteam*masssteam;
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Inlet pressure " << referencepressure;
        sout << "\nMass flux " << massflux;
        sout << "\ntemperature " << temperature;
        sout << "\nenthalpywater " << enthalpywater;
        sout << "\nenthalpysteam " << enthalpysteam;
        sout << "\nmass flux water " << masswater;
        sout << "\nmass flux steam " << masssteam;
        sout << "\nenergy flux " << energyflux;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPBrSteamFlux::fInletEnergyFlux = energyflux;
}

/// Limit correction range
REAL TPBrSteamMesh::LimitRange(TPZFMatrix<REAL> &prevsol, TPZFMatrix<REAL> &correction)
{
    TPZManVector<REAL> facestate,cellstate,inletstate,facecorrection,cellcorrection,inletcorrection;
    REAL scale = 1.;
    ExtractInletfaceState(prevsol, inletstate);
    ExtractInletfaceState(correction, inletcorrection);
    ExtractCellState(0, prevsol, cellstate);
    ExtractCellState(0, correction, cellcorrection);
    scale = TPBrSteamFlux::LimitRangeInlet(scale,inletstate,cellstate,inletcorrection,cellcorrection);
    int iface;
    for (iface=1; iface<fNumCells; iface++) {
        ExtractInterfaceState(0,prevsol, facestate);
        ExtractInterfaceState(0,correction, facecorrection);
        scale = TPBrSteamFlux::LimitRange(scale,facestate,facecorrection);
    }
    for (int icell=0; icell<fNumCells; icell++) {
        ExtractCellState(icell, prevsol, cellstate);
        ExtractCellState(icell, correction, cellcorrection);
        scale = TPBrCellConservation::LimitRange(scale,cellstate,cellcorrection);
    }
    return scale;
}

/// Project the solution to an allowable state
void TPBrSteamMesh::ProjectSolution(TPZFMatrix<REAL> &solution)
{
    TPZManVector<REAL> facestate,cellstate,inletstate,facecorrection,cellcorrection,inletcorrection;
    ExtractInletfaceState(solution, facestate);
    ExtractCellState(0, solution, cellstate);
    if(facestate[TPBrSteamFlux::EInletPressure] - cellstate[TPBrCellConservation::EPressureWater] < 100.)
    {
        solution(TPBrSteamFlux::EInletPressure,0) = cellstate[TPBrCellConservation::EPressureWater] + 100.;
    }
    if(facestate[TPBrSteamFlux::EInletSteamSaturation] > 0.01)
    {
        REAL press = facestate[TPBrSteamFlux::EInletPressure]+TPBrScales::fReferencePressure;
        REAL sattemp = TPBrCellConservation::TemperatureSaturation(press);
        solution(TPBrSteamFlux::EInletTemperature,0) = sattemp;
    }
    else if(facestate[TPBrSteamFlux::EInletSteamSaturation] < 0.)
    {
        solution(TPBrSteamFlux::EInletSteamSaturation,0) = 0.;
    }
    if(facestate[TPBrSteamFlux::EInletSteamSaturation] > 1.)
    {
        solution(TPBrSteamFlux::EInletSteamSaturation,0) = 1.;
    }
    for (int icell=0; icell<fNumCells; icell++) {
        ExtractCellState(icell, solution, cellstate);
        if(cellstate[TPBrCellConservation::ESaturationSteam] > 0.01)
        {
            REAL press = cellstate[TPBrCellConservation::EPressureWater]+TPBrScales::fReferencePressure;
            REAL temp = TPBrCellConservation::TemperatureSaturation(press);
            TPZManVector<int> states,eqs;
            CellDestination(icell, eqs, states);
            solution(eqs[TPBrCellConservation::ETemperature],0) = temp;
        }
        else if(cellstate[TPBrCellConservation::ESaturationSteam] < 0.)
        {
            REAL watersat = cellstate[TPBrCellConservation::ESaturationWater];
            REAL oilsat = cellstate[TPBrCellConservation::ESaturationOil];
            REAL sumsat = watersat+oilsat;
            TPZManVector<int> states,eqs;
            CellDestination(icell, eqs, states);
            solution(eqs[TPBrCellConservation::ESaturationWater],0) = watersat/sumsat;
            solution(eqs[TPBrCellConservation::ESaturationOil],0) = oilsat/sumsat;
            solution(eqs[TPBrCellConservation::ESaturationSteam],0) = 0.;
        }
    }    
}

/// verify the consistency of the tangent matrix
void TPBrSteamMesh::VerifyTangent(TPZFMatrix<REAL> &direction)
{
	int neq = NumEquations();
	TPZFBMatrix<REAL> bnd(neq,Bandwidth()),bnd2(neq,Bandwidth());
	TPZFNMatrix<100> rhsref(neq,1),rhs(neq,1),result(neq,1);
	TPZManVector<REAL> scales(neq,1.);
	ComputeTangent2(bnd, rhsref,scales);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        bnd.Print("VerifyTangent matrix ",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	TPZFNMatrix<100> refsol = fNextState;
	TPZManVector<REAL,10> norms(10);
	TPZFNMatrix<100> dirkeep(direction);
	REAL alfa = 0.1;
	int count = 0;
	for (alfa = 0.1; alfa < 1.; alfa += 0.1, count++) {
		TPZFNMatrix<100> alfadir = alfa*direction;
		fNextState = refsol+alfadir;
		ComputeTangent2(bnd2, rhs, scales);
		TPZFNMatrix<100> diff;
		diff = rhs-rhsref;
		bnd.Multiply(alfadir,result);
		diff -= result;
		norms[count] = Norm(diff);
	}
	std::cout << "errors of the tangent prediction " << norms << std::endl;
	fNextState = refsol;
}

/// verify the conservation of energy between prevstate and nextstate
void TPBrSteamMesh::EnergyBalance(REAL &leftin, REAL &rightout, REAL &energychange)
{
    leftin = fNextState(TPBrSteamFlux::NumInletVars+TPBrSteamFlux::EEnergyFlux,0)*fDelt;
    TPZManVector<REAL> facestate,cellstate,inletstate,facecorrection,cellcorrection,inletcorrection;
    ExtractInterfaceState(fNumCells,fNextState, facestate);
    rightout = facestate[TPBrSteamFlux::EEnergyFlux]*fDelt;
    REAL energybefore = 0.;
    REAL energyafter = 0.;
    for (int icell=0; icell<fNumCells; icell++) {
        REAL leftco = NodeCoord(icell);
		REAL rightco = NodeCoord(icell+1);
		REAL volume = M_PI*(rightco*rightco-leftco*leftco)*fHeight;

        ExtractCellState(icell, fPrevState, cellstate);
        energybefore += TPBrCellConservation::Energy(cellstate,volume);
        ExtractCellState(icell, fNextState, cellstate);
        energyafter += TPBrCellConservation::Energy(cellstate,volume);
    }
    energychange = energyafter-energybefore;
}

/// write the steam saturation to a file
void TPBrSteamMesh::PrintSaturations(std::ostream &out)
{
    TPZManVector<REAL> facestate,cellstate,inletstate,facecorrection,cellcorrection,inletcorrection;
    for (int icell=0; icell<fNumCells; icell++) {
        REAL leftco = NodeCoord(icell);
		REAL rightco = NodeCoord(icell+1);
        ExtractCellState(icell, fNextState, cellstate);
        out << (leftco+rightco)/2. << " " << cellstate[TPBrCellConservation::ESaturationSteam] << std::endl;
    }
}

/// write the temperature to a file
void TPBrSteamMesh::PrintTemperature(std::ostream &out)
{
    TPZManVector<REAL> facestate,cellstate,inletstate,facecorrection,cellcorrection,inletcorrection;
    for (int icell=0; icell<fNumCells; icell++) {
        REAL leftco = NodeCoord(icell);
		REAL rightco = NodeCoord(icell+1);
        ExtractCellState(icell, fNextState, cellstate);
        out << (leftco+rightco)/2. << " " << cellstate[TPBrCellConservation::ETemperature] << std::endl;
    }
    
}

/// write the pressure to a file
void TPBrSteamMesh::PrintPressure(std::ostream &out)
{
    TPZManVector<REAL> facestate,cellstate,inletstate,facecorrection,cellcorrection,inletcorrection;
    for (int icell=0; icell<fNumCells; icell++) {
        REAL leftco = NodeCoord(icell);
		REAL rightco = NodeCoord(icell+1);
        ExtractCellState(icell, fNextState, cellstate);
        out << (leftco+rightco)/2. << " " << cellstate[TPBrCellConservation::EPressureWater] << std::endl;
    }
    
}

/// write all solutions for a given time
void TPBrSteamMesh::PrintAll(REAL time)
{
    std::stringstream satname,tempname,pressurename;
    satname << "saturation t = " << time << ".txt";
    tempname << "temperature t = " << time << ".txt";
    pressurename << "pressure t = " << time << ".txt";
    ofstream satfile(satname.str().c_str());
    ofstream tempfile(tempname.str().c_str());
    ofstream presfile(pressurename.str().c_str());
    PrintPressure(presfile);
    PrintSaturations(satfile);
    PrintTemperature(tempfile);

}

/// initialize the data structure of the thermal problem
void TPBrSteamMesh::InitializeThermalProblem()
{
    fThermal.ClearSolutions();
    int icell;
    for (icell=0; icell<fNumCells; icell++) {
        REAL node0 = NodeCoord(icell);
        REAL node1 = NodeCoord(icell+1);
        REAL area = M_PI*(node1*node1-node0*node0);
        TPBRThermalSolution tsol(area);
        fThermal.AddSolution(tsol);
    }
}

/// Update the solution of the confinement layer
void TPBrSteamMesh::UpdateConfinementLayer()
{
    TPZManVector<REAL> facestate,cellstate,inletstate,facecorrection,cellcorrection,inletcorrection;
    for (int icell=0; icell<fNumCells; icell++) {
        ExtractCellState(icell, fNextState, cellstate);
        REAL temp = cellstate[TPBrCellConservation::ETemperature];
        REAL flux,dqdt;
        fThermal.AdvanceSolution(icell, temp, flux, dqdt, true);
    }

    
}

#endif
